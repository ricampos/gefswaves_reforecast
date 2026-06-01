#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
wcpval.py

VERSION AND LAST UPDATE:
 v1.0  10/30/2023
 v1.1  11/03/2023
 v1.2  05/21/2024
 v1.3  01/10/2025
 v1.4  05/01/2025
 v1.5  07/08/2025

PURPOSE:
 Validation of long-term probabilistic wave forecasts using 
 fuzzy verification. This can be used to assess and optimize
 the parameters of operational probability maps.
 The probability maps are based on GEFSv12 and have been validated against
  NDBC buoys, GDAS, and GEFS hindcasts by appending 24-hour segments 
  from consecutive cycles.

USAGE:
 These are auxiliar functions used by the main script fuzzy_verification_GEFS.py

DEPENDENCIES:
 See the imports below.

AUTHOR and DATE:
 10/30/2023: Ricardo M. Campos, first version, named wprob.py
 11/03/2023: Ricardo M. Campos, fuzzy verification including categorical 
  and probabilistic assessments.
 05/21/2024: Ricardo M. Campos, renamed to wcpval.py, new reliability curve plot
 01/10/2025: Ricardo M. Campos, climatology and persistence added to reliability diagrams and ROC
  plots. New updates added to the reliability diagram plots.
 05/01/2025: Ricardo M. Campos, edit in read_data (dependence on NDBC buoy data has been removed)
 07/08/2025: Ricardo M. Campos, bias correction applied to the hindcast used as ground truth

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
# matplotlib.use('Agg')
from matplotlib.dates import DateFormatter
import netCDF4 as nc
import xarray as xr
import numpy as np
import pandas as pd
import pickle
import sys
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sklearn.metrics import roc_curve, auc
from sklearn.calibration import calibration_curve
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import properscoring as ps
from datetime import datetime
import yaml
import mvalstats
import pvalstats
from pvalstats import ModelObsPlot
from lmoments3 import lmom_ratios
import warnings; warnings.filterwarnings("ignore")

sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

pcolors = np.array(['navy', 'firebrick', 'dimgrey', 'darkgreen', 'fuchsia', 'gold', 'blue', 'salmon', 'lime', 'darkviolet', 'yellow',
    'cornflowerblue', 'red', 'green', 'violet', 'orange'])

lstl  = np.array(['-','--','-.'])

# READ DATA
def read_data(wlist,bid,ltime1,ltime2,wvar):

    lstw=int(len(wlist))
    f=nc.Dataset(wlist[0])    
    # stations = np.array(f.variables['ID'][bid]).astype('str')
    stations = np.array([f.variables['ID'][i] for i in bid]).astype('str')
    latm = f.variables['lat'][bid,:]; lonm = f.variables['lon'][bid,:]
    ensm = f.variables['ensemble_member'][:]
    indlat = int(np.floor(len(latm[0,:])/2))
    indlon = int(np.floor(len(lonm[0,:])/2))
    auxt = np.array(f.variables['time'][:])
    # week 2 forecast
    axt = np.array(auxt-auxt.min())/3600.
    indt = np.where( (axt>=ltime1*24) & (axt<=ltime2*24) )[0]
    f.close(); del f, auxt, axt

    # read data
    for i in range(0,lstw):

        f=nc.Dataset(wlist[i])    

        if i==0:
            ctime = np.zeros((len(wlist)),'double')*np.nan
            ftime = np.zeros((len(wlist),len(indt)),'double')*np.nan
            fshape = f.variables[wvar+"_gefs_forecast"].shape
            gefs_hindcast = np.zeros((len(wlist),len(bid),len(indt),fshape[2],fshape[3],fshape[4]),'f')*np.nan
            gefs_forecast = np.zeros((len(wlist),len(bid),len(indt),fshape[2],fshape[3],fshape[4]),'f')*np.nan
            gefs_forecast_p = np.zeros((len(wlist),len(bid),len(indt),fshape[2],fshape[3],fshape[4]),'f')*np.nan
            del fshape

        ftime[i,:] = np.array(f.variables['time'][indt]).astype('double')
        ctime[i] = np.double(f.variables['time'][0])
        gefs_hindcast[i,:,:,:,:,:] = np.array(f.variables[wvar+"_gefs_hindcast"][bid,:,:,:,:][:,indt,:,:,:])
        gefs_forecast[i,:,:,:,:,:] = np.array(f.variables[wvar+"_gefs_forecast"][bid,:,:,:,:][:,indt,:,:,:])
        aux = np.array(f.variables[wvar+"_gefs_forecast"][bid,:,:,:,:][:,0,:,:,:])
        for j in range(0,len(indt)):
            gefs_forecast_p[i,:,j,:,:,:] = aux

        f.close(); del f, aux
        print(" read file "+repr(i)+" of "+repr(lstw-1))

    cdate = [datetime.utcfromtimestamp(ts) for ts in ctime]

    gefsdata = {'cdate': cdate,'ctime': ctime,'ftime': ftime, 'stations': stations, 'lft': len(indt), 'ensm': ensm,
        'latm': latm, 'lonm': lonm, 'indlat': indlat, 'indlon': indlon, 
        'indlat': indlat, 'indlon': indlon, 'indt': indt,
        'gefs_hindcast': gefs_hindcast,'gefs_forecast': gefs_forecast,'gefs_forecast_p': gefs_forecast_p}

    return gefsdata

# READ OBSERVATIONS 
def read_obs(fname,wvar,altselm):
    '''
    wvar: hs u10
    altselm: mean, p80, p90, p95, p95, max
    '''
    ds = xr.open_dataset(fname)
    f = nc.Dataset(fname)

    result = {'bdate': ds.time.values,'btime': f.variables['time'][:], 'bid': f.variables['ID'][:],
        'lat': f.variables['lat'][:], 'lon': f.variables['lon'][:], 
        'mlat': f.variables['mlat'][:], 'mlon': f.variables['mlon'][:], 
        'buoy': f.variables["buoy_"+wvar][:], 'alt': f.variables["alt_"+wvar+"_"+altselm][:]}

    return result


# --- Bias correction using altimeter obs ---

# Quantile Mapping bias-correction
def qm_train(model=None,obs=None,prob=None):
    """ 
    Univariate linear regression calibration using the Quantile Mapping Method
    Fit module
    Input:
    - model data
    - observations ("truth")
    - probability array (optional) to define the percentiles 
    Output:
    - slope
    - intercept
    """

    if (np.size(model)>2)==False or (np.size(obs)>2)==False:
        raise ValueError("Model and Obs arrays must be informed.")

    if np.size(model) != np.size(obs):
        raise ValueError("Model and Obs arrays must have the same sizes.")

    # Probability array
    if prob==None:
        prob=np.append(np.array(np.arange(1.,50.+0.1,1.),'f'),np.array(np.arange(50.5,99.5+0.1,0.5),'f'))
        prob=np.append(prob,np.array([99.7,99.8,99.9,99.95,99.99]))

    slope,intercept = np.polyfit(np.nanpercentile(model,prob),np.nanpercentile(obs,prob),1)
    # model_cal=np.array((model*slope)+intercept)
    print(" QMM linear regression. Slope: "+repr(slope)+"  Intercept: "+repr(intercept))

    return float(slope),float(intercept)

def qmcal(model=None,slope=1.,intercept=0.,pprint='yes'):
    """ 
    Univariate linear regression calibration using the Quantile Mapping Method
    Calibration module based on previously trained qm_train
    Input:
    - model data
    - slope
    - intercept
    Output:
    - calibrated model data
    """

    if (np.size(model)>2)==False:
        raise ValueError("Model array must be informed.")

    model_cal=np.array((model*slope)+intercept)
    if pprint=='yes':
        print(" QMM linear regression. Slope: "+repr(slope)+"  Intercept: "+repr(intercept))

    return np.array(model_cal).astype('float')


def bias_correction(gdata,obs,spctl,wvar,opath):
    '''
    Bias correction of GEFS hindcast data. 
    Using altimeter data for calibration, and buoy data to validate the results.
    '''

    gefs_hindcast = gdata['gefs_hindcast']

    gefs_hindcast_bc = np.copy(gefs_hindcast)
    pctlarr = np.arange(spctl,99.9,0.1)

    merr=np.zeros((len(gdata['stations']),len(gdata['ensm'])),'f')*np.nan
    merr_cal=np.zeros((len(gdata['stations']),len(gdata['ensm'])),'f')*np.nan

    for i in range(0,len(gdata['stations'])):
        if gdata['stations'][i] in obs['bid']:
            inds=np.where(gdata['stations'][i]==obs['bid'])[0][0]
            fmodel=np.array([]); fobs=np.array([])
            c=0
            for j in range(0,gdata['ftime'].shape[0]):
                for k in range(0,gdata['ftime'].shape[1]):
                    indt = np.where( np.abs(gdata['ftime'][j,k]-obs['btime'])<1800. )
                    if np.size(indt)>0:
                        indt=indt[0][0]
                        if np.any( np.isnan(obs['alt'][inds,:,:][indt,:]) == False ):

                            for p in range(0,obs['mlat'].shape[1]):
                                if obs['alt'][inds,:,:][indt,:][p] > 0. :

                                    indplat=np.where( np.abs(gdata['latm'][i,:]-obs['mlat'][inds,p]) == np.min( np.abs(gdata['latm'][i,:]-obs['mlat'][inds,p]) ) )[0][0]
                                    indplon=np.where( np.abs(gdata['lonm'][i,:]-obs['mlon'][inds,p]) == np.min( np.abs(gdata['lonm'][i,:]-obs['mlon'][inds,p]) ) )[0][0]

                                    # statistical comparison between neighboring points
                                    stresult = np.round(np.abs(mvalstats.metrics(gefs_hindcast[j,i,:,0,gdata['indlat'],gdata['indlon']],gefs_hindcast[j,i,:,0,indplat,indplon])[[2,5,7]]),3)
                                    L1, L2, t3_p, t4_p = lmom_ratios(gefs_hindcast[j,i,:,0,gdata['indlat'],gdata['indlon']], nmom=4); l1_p = L2 / L1
                                    L1, L2, t3_n, t4_n = lmom_ratios(gefs_hindcast[j,i,:,0,indplat,indplon], nmom=4); l1_n = L2 / L1

                                    # if stresult[0]<0.1 and stresult[1]<0.1 and stresult[2]>0.8 and np.abs(l1_p-l1_n)<0.05 and np.abs(t3_p-t3_n)<0.05 and np.abs(t4_p-t4_n)<0.05:
                                    if stresult[0]<0.15 and stresult[1]<0.15 and stresult[2]>0.75 and np.abs(l1_p-l1_n)<0.07 and np.abs(t3_p-t3_n)<0.07 and np.abs(t4_p-t4_n)<0.07:

                                        if c==0:
                                            fmodel=np.array([gefs_hindcast[j,i,k,:,indplat,indplon]])
                                            fobs=np.array([np.full( len(gdata['ensm']),obs['alt'][inds,:,:][indt,:][p])])
                                        else:
                                            fmodel=np.append(fmodel, np.array([gefs_hindcast[j,i,k,:,indplat,indplon]]),axis=0)
                                            fobs=np.append(fobs, np.array([np.full( len(gdata['ensm']),obs['alt'][inds,:,:][indt,:][p])]),axis=0)

                                        del indplat,indplon
                                        c=c+1
                                        # print(" Ok, "+gdata['stations'][i]+". "+repr(j)+", "+repr(k))

                                    else:
                                        print(" Neighbour point is statistically different, "+gdata['stations'][i]+". "+repr(j)+", "+repr(k))

            for j in range(0,len(gdata['ensm'])):
                if np.size(np.where(fobs[:,j]>0))>100.:
                    # qmm_slope,qmm_intercept = qm_train(model=fmodel[:,j],obs=fobs[:,j])
                    qmm_slope=1.
                    qmm_intercept = np.nanmean( np.nanpercentile(fobs[:,j],pctlarr) - np.nanpercentile(fmodel[:,j],pctlarr) )
                    # merr = np.array(mvalstats.metrics(fmodel[:,j],fobs[:,j],pctlerr='yes'))[9]

                    gefs_hindcast_bc[:,i,:,:,:,:][:,:,j,:,:] = qmcal(model=gefs_hindcast[:,i,:,:,:,:][:,:,j,:,:],slope=qmm_slope,intercept=qmm_intercept)
                    del qmm_slope,qmm_intercept

                    print(" QMM station "+gdata['stations'][i])

            # Quick assessment against buoy (when available)
            bobs=[]
            c=0
            if np.any(obs['buoy'][inds,:] > 0.)==True:
                for j in range(0,gdata['ftime'].shape[0]):
                    for k in range(0,gdata['ftime'].shape[1]):
                        indt = np.where( np.abs(gdata['ftime'][j,k]-obs['btime'])<1800. )
                        if np.size(indt)>0:
                            if c==0:
                                ghcal = np.array([gefs_hindcast_bc[j,i,k,:,gdata['indlat'],gdata['indlon']]])
                                gh = np.array([gefs_hindcast[j,i,k,:,gdata['indlat'],gdata['indlon']]])
                            else:
                                ghcal = np.append(ghcal,np.array([gefs_hindcast_bc[j,i,k,:,gdata['indlat'],gdata['indlon']]]),axis=0)
                                gh = np.append(gh,np.array([gefs_hindcast[j,i,k,:,gdata['indlat'],gdata['indlon']]]),axis=0)

                            bobs=np.append(bobs,obs['buoy'][inds,indt])
                            c=c+1

                for j in range(0,len(gdata['ensm'])):
                    if np.any(gh[:,j]>0.)==True and np.any(ghcal[:,j]>0.)==True and np.any(bobs>0.)==True :
                        # RMSE pctl95
                        merr[i,j] = mvalstats.metrics(gh[:,j],bobs,pctlerr='yes')[9]
                        merr_cal[i,j] = mvalstats.metrics(ghcal[:,j],bobs,pctlerr='yes')[9]
                        # Plots
                        if j==0:
                            mop=ModelObsPlot(model=[gh[:,j],ghcal[:,j]],obs=bobs,mlabels=["GEFSv12","GEFSv12_BC"],ftag=opath+"/Eval_"+wvar+"_"+gdata['stations'][i]+"_",vaxisname=wvar)
                            mop.qqplot(); mop.scatterplot(); mop.taylordiagram()

                del ghcal, gh

            del bobs

            if np.any(merr[i,:]>-999.)==True and np.any(merr_cal[i,:]>-999.)==True :
                fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
                ax.plot(merr[i,:],'k'); ax.plot(merr_cal[i,:],'r')
                plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
                ax.set_xlabel("Ens Members"); ax.set_ylabel("Bias (Pctl95)")
                plt.tight_layout()
                plt.savefig(opath+"/Eval_"+wvar+"_"+gdata['stations'][i]+"_BiasPctl95.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
                plt.close(fig1); del fig1, ax


    return gefs_hindcast_bc


# ---- STATISTICAL MODELLING -----

# Function to check np.sort using multidimensional array
def tsort(model):

    fmodel=np.copy(model)
    if len(fmodel.shape) == 4:
        for i in range(0,model.shape[1]):
            for j in range(0,model.shape[2]):
                for k in range(0,model.shape[3]):
                    fmodel[:,i,j,k] = np.sort(model[:,i,j,k])

    elif len(fmodel.shape) == 3:
        for i in range(0,model.shape[1]):
            for j in range(0,model.shape[2]):
                fmodel[:,i,j] = np.sort(model[:,i,j])

    elif len(fmodel.shape) == 2:
        for i in range(0,model.shape[1]):
            fmodel[:,i] = np.sort(model[:,i])

    elif len(fmodel.shape) == 2:
        fmodel = np.sort(model[:])

    return fmodel


def nmaxsel(wdata,nmax):

    lstw = int(wdata.shape[0])
    wdata_tmax = np.mean(wdata,axis=1)
    for i in range(0,lstw):
        ind = np.where(wdata[i,:]>-999.)
        if np.size(ind)>int(np.floor(wdata.shape[1]/2)):
            ind = ind[0]
            wdata_tmax[i] = np.nanmean(np.sort(wdata[i,ind])[-nmax::])
        else:
            print(" nmaxsel(): forecast time series incomplete "+repr(i))
            wdata_tmax[i] = np.nan

        del ind

    return wdata_tmax


def probforecast(nmax,gspws,spctl,mod_forecast,qlev):

    lstw=int(mod_forecast.shape[0])
    nensm=int(mod_forecast.shape[2])
    indlat = int(np.floor(mod_forecast.shape[3]/2)); indlon = int(np.floor(mod_forecast.shape[4]/2))

    prob_forecast = np.zeros((lstw,len(qlev)),'f')*np.nan
    for i in range(0,lstw):
        fmod=np.array(mod_forecast[i,:,:,(indlat-gspws):(indlat+gspws+1),:][:,:,:,(indlon-gspws):(indlon+gspws+1)])
        fmod=tsort(fmod)
        fmod=np.copy(fmod[-nmax::,:,:,:])
        fmod=np.array(fmod.reshape(nmax*nensm,fmod.shape[2],fmod.shape[3]))
        fmod=np.array(fmod.reshape(nmax*nensm,fmod.shape[1]*fmod.shape[2]))
        fmod=np.sort(fmod,axis=1)
        ind=np.where(np.mean(fmod,axis=0)>=0.)
        if np.size(ind)>0:
            fmod=np.array(fmod[:,ind[0]])            
            fmod=np.array(fmod[:,int(np.floor(fmod.shape[1]*(spctl/100)))::])
            fmod=fmod.reshape(fmod.shape[0]*fmod.shape[1])
            for j in range(0,len(qlev)):
                prob_forecast[i,j] = np.size(fmod[fmod>qlev[j]]) / np.size(fmod[fmod>0.])

        if i==0:
            fmod_result = np.array([fmod])
        else:
            fmod_result = np.append(fmod_result,np.array([fmod]),axis=0)

    return prob_forecast, fmod_result


def prob_validation(nmax,spws,gspws,spctl,prob_u10_gefs_forecast,fmod_result_u10,prob_hs_gefs_forecast,fmod_result_hs,u10_obs_tmax,hs_obs_tmax,qlev_wnd,qlev_hs,plevels,ftag):

    # Binary Categorical array
    # U10
    fprob_u10_gefs_forecast = probforecast_binary(prob_u10_gefs_forecast,qlev_wnd,plevels)
    # Hs
    fprob_hs_gefs_forecast = probforecast_binary(prob_hs_gefs_forecast,qlev_hs,plevels)

    print(" Prob Forecast and Binary array, OK")

    lplev=int(len(plevels)-1)

    # ----------- Error metrics: TP, FP, TN, FN, CSI;  Brier; ROC_AUC ---------------
    ceval_gefs_u10 = np.zeros((len(qlev_wnd),lplev,3),'f')*np.nan
    for i in range(0,len(qlev_wnd)):
        for j in range(0,lplev):
            bresult = categorical_bin_eval(fprob_u10_gefs_forecast[:,i,j],u10_obs_tmax,qlev_wnd[i])
            ceval_gefs_u10[i,j,0] = float(bresult['POD'])
            ceval_gefs_u10[i,j,1] = float(bresult['FAR'])
            ceval_gefs_u10[i,j,2] = float(bresult['CSI'])
            del bresult

    ceval_gefs_hs = np.zeros((len(qlev_hs),lplev,3),'f')*np.nan
    for i in range(0,len(qlev_hs)):
        for j in range(0,lplev):
            bresult = categorical_bin_eval(fprob_hs_gefs_forecast[:,i,j],hs_obs_tmax,qlev_hs[i])
            ceval_gefs_hs[i,j,0] = float(bresult['POD'])
            ceval_gefs_hs[i,j,1] = float(bresult['FAR'])
            ceval_gefs_hs[i,j,2] = float(bresult['CSI'])
            del bresult

    # --- ROC Curve ---
    # U10
    froc_u10=[]
    for i in range(0,len(qlev_wnd)):
        fftag=ftag+"ProbEvents_U10_"+repr(np.round(qlev_wnd[i],2))
        true_binary = (u10_obs_tmax > qlev_wnd[i]).astype(int)
        roc = roc_plot(true_binary,prob_u10_gefs_forecast[:,i],fftag)
        froc_u10 = np.append(froc_u10,roc); del true_binary, fftag, roc

    # Hs
    froc_hs=[]
    for i in range(0,len(qlev_hs)):
        fftag=ftag+"ProbEvents_Hs_"+repr(np.round(qlev_hs[i],2))
        true_binary = (hs_obs_tmax > qlev_hs[i]).astype(int)
        roc = roc_plot(true_binary,prob_hs_gefs_forecast[:,i],fftag)
        froc_hs = np.append(froc_hs,roc); del true_binary, fftag, roc

    print(" Validation: Categorical OK")
    # --- Brier Score ---
    # U10
    fbriers_u10=[]
    for i in range(0,len(qlev_wnd)):    
        fftag=ftag+"ProbEvents_U10_"+repr(np.round(qlev_wnd[i],2))
        briers = brier_score(prob_u10_gefs_forecast[:,i],u10_obs_tmax,qlev_wnd[i],fftag)
        fbriers_u10 = np.append(fbriers_u10,briers)
        del briers, fftag

    # Hs
    fbriers_hs=[]
    for i in range(0,len(qlev_hs)):    
        fftag=ftag+"ProbEvents_Hs_"+repr(np.round(qlev_hs[i],2))
        briers = brier_score(prob_hs_gefs_forecast[:,i],hs_obs_tmax,qlev_hs[i],fftag)
        fbriers_hs = np.append(fbriers_hs,briers)
        del briers, fftag

    print(" Validation: Brier Score OK")

    # --- CRPS ---
    # U10 
    fftag=ftag+"ProbEvents_U10"
    crps_u10, mean_crps_u10 = crps(fmod_result_u10, u10_obs_tmax, fftag)
    del fftag
    # Hs
    fftag=ftag+"ProbEvents_hs"
    crps_hs, mean_crps_hs = crps(fmod_result_hs, hs_obs_tmax, fftag)
    del fftag

    plt.close('all')
    print(" Validation: CRPS OK")

    # ================================================

    # Save results
    stat_results = {
    'ceval_gefs_hs': ceval_gefs_hs,
    'ceval_gefs_u10': ceval_gefs_u10,
    'froc_hs': froc_hs,
    'froc_u10': froc_u10,
    'fbriers_u10': fbriers_u10,
    'fbriers_hs': fbriers_hs,
    'crps_u10': crps_u10,
    'mean_crps_u10': mean_crps_u10,
    'crps_hs': crps_hs,
    'mean_crps_hs': mean_crps_hs
    }

    # Save the dictionary to a pickle file
    with open(ftag+'STAT.RESULTS.pkl', 'wb') as f:
        pickle.dump(stat_results, f)

    return stat_results
    print(" Results saved in the pickle file STAT.RESULTS.pkl")

    del u10_obs_tmax, hs_obs_tmax, prob_u10_gefs_forecast, prob_hs_gefs_forecast, fprob_u10_gefs_forecast, fprob_hs_gefs_forecast
    del stat_results, ceval_gefs_u10, ceval_gefs_hs, froc_u10, froc_hs, fbriers_u10, fbriers_hs, crps_u10, mean_crps_u10, crps_hs, mean_crps_hs

    print(" - Done Optmz_nmax"+str(int(nmax))+"_spws"+str(int(spws*100)).zfill(3)+"_spctl"+str(int(spctl)))


# Binary Categorical array
def probforecast_binary(prob_forcast,qlev,plevels):

    lstw = int(prob_forcast.shape[0])
    lplev = int(len(plevels)-1)
    fprob_forecast = np.zeros((lstw,len(qlev),lplev),'f')*np.nan

    for i in range(0,lstw):
        for j in range(0,len(qlev)):
            for k in range(0,lplev):
                if prob_forcast[i,j] >= plevels[k]:
                    fprob_forecast[i,j,k] = int(1)
                else:
                    fprob_forecast[i,j,k] = int(0)

    return fprob_forecast


def categorical_bin_eval(forecasted,obs,threshold):

    ind = np.where((forecasted>-999.)&(obs>-999.))
    if np.size(ind)>1:
        forecasted=forecasted[ind[0]]
        obs=obs[ind[0]]

        SIZE = int(len(obs))
        EVENTS = int(np.size(np.where((obs >= threshold)==True)))

        # Calculate True Positives (TP), False Positives (FP), False Negatives (FN), and True Negatives (TN)
        if np.any((obs >= threshold)==True):
            TP = np.sum((forecasted == int(1)) & (obs >= threshold))
            FP = np.sum((forecasted == int(1)) & (obs < threshold))
            FN = np.sum((forecasted == int(0)) & (obs >= threshold))
            TN = np.sum((forecasted == int(0)) & (obs < threshold))

            # POD (Probability of Detection) also name as True Positive Rate (TPR) or Sensitivity
            if (TP + FN) == 0:
                POD = np.nan
            else:
                POD = np.round((TP / (TP + FN)),3)

            # FAR (False Alarm Ratio)
            if (FP + TP)==0:
                FAR = 0
            else:
                FAR = np.round((FP / (FP + TP)),3)

            # Critical Success Index (CSI)
            if (TP + FN + FP)==0:
                CSI = np.nan
            else:
                CSI = np.round((TP / (TP + FN + FP)),3)

            # ETS (Equitable Threat Score)
            ETS = np.round(((TP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN)) / (TP + FN + FP - (TP + FN) * (TP + FP) / (TP + FN + FP + TN))) ,3)
        else:
            TP=np.nan; FP=np.nan; FN=np.nan; TN=np.nan
            POD=np.nan; FAR=np.nan; CSI=np.nan; ETS=np.nan

        result = {'SIZE':SIZE, 'EVENTS':EVENTS,
            'TP':TP, 'FP':FP, 'FN':FN, 'TN':TN,
            'POD':POD, 'FAR':FAR, 'CSI':CSI, 'ETS': ETS}

        return result

def brier_score(prob_forecast,obs,threshold):

    # Brier Score
    true_binary = (obs > threshold).astype(int)
    briers=0
    for i in range(0,len(true_binary)):
        briers = briers+(true_binary[i]-prob_forecast[i])**2

    briers = briers/len(true_binary)

    return float(briers)


def crps(fmod,obs):

    # CRPS
    batch_size = 1000
    if fmod.shape[0]>batch_size:
        crps=[]
        # Process in batches
        for i in range(0, fmod.shape[0], batch_size):
            fmod_batch = fmod[i:i+batch_size]
            obs_batch = obs[i:i+batch_size]
            # Compute CRPS for this batch
            crps_batch = ps.crps_ensemble(obs_batch, fmod_batch)
            crps.append(crps_batch)

        crps = np.array(np.concatenate(crps))
    else:
        crps = ps.crps_ensemble(obs,fmod)

    fcrps_mean = np.nanmean(crps)

    return np.array(crps).astype('float'), float(fcrps_mean)
    del crps, fcrps_mean


def roc_plot(true_binary,predicted_probabilities,ftag):

    # Calculate ROC curve
    fpr, tpr, _ = roc_curve(true_binary, predicted_probabilities)
    # Calculate AUC (Area Under the ROC Curve)
    roc_auc = auc(fpr, tpr)

    if len(ftag)>0:
        # Plot the ROC curve
        fig1, ax = plt.subplots(figsize=(5,4.5))
        ax.plot([0, 1], [0, 1], color='dimgray', lw=2, linestyle='--', alpha=0.8, zorder=1)
        ax.plot(fpr, tpr, color='firebrick', lw=1, marker='.', alpha=0.8, zorder=2)
        ax.plot(fpr, tpr, color='firebrick', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})', zorder=3)
        ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01) 
        ax.set_xlabel('False Positive Rate')
        ax.set_ylabel('True Positive Rate')
        ax.legend(loc='lower right',fontsize=sl-3)
        plt.title('ROC Curve')
        plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
        plt.savefig(ftag+"_ROC.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
            format='png', bbox_inches='tight', pad_inches=0.1)

        plt.close(fig1)

    return float(roc_auc)


def roc_plot2(true_binary,predicted_probabilities,mlabels,plevels,ftag):

    from sklearn.metrics import confusion_matrix, accuracy_score

    predicted_probabilities = np.atleast_2d(predicted_probabilities)
    froc_auc = []

    fig1, ax = plt.subplots(figsize=(5,4.5))
    ax.plot([0, 1], [0, 1], color='dimgray', lw=2, linestyle='--', alpha=0.8, zorder=1)

    for i in range(0,predicted_probabilities.shape[0]):

        # Calculate ROC curve
        fpr, tpr, _ = roc_curve(true_binary, predicted_probabilities[i,:])

        ## Calculate confusion matrix
        for j in range(0,len(plevels)):
            y_pred = (predicted_probabilities[i,:] >= plevels[j]).astype(int)
            tn, fp, fn, tp = confusion_matrix(true_binary,y_pred).ravel()
            pod = tp / (tp + fn)
            far = fp / (fp + tp)
            csi = tp / (tp + fp + fn)
            print("CSI prob lev "+str(plevels[j])+" = "+repr(np.round(csi,3)))

        print("--------"+repr(i))

        # Calculate AUC (Area Under the ROC Curve)
        roc_auc = auc(fpr, tpr)
        froc_auc = np.append(froc_auc,roc_auc)

        # Plot the ROC curve
        if i==0:
            ax.plot(fpr, tpr, color=pcolors[i], lw=1, marker='.', alpha=0.8, zorder=2)      
            ax.plot(fpr, tpr, color=pcolors[i], lw=2, linestyle=lstl[i], label=f'{mlabels[i]} (AUC = {roc_auc:.2f})', zorder=3)
        else:
            ax.plot(fpr, tpr, color=pcolors[i], lw=1, linestyle=lstl[i], label=f'{mlabels[i]} (AUC = {roc_auc:.2f})',alpha=0.7 , zorder=3)

    ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.01) 
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc='lower right',fontsize=sl-3)
    # plt.title('ROC Curve')
    plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
    plt.savefig(ftag+"_ROC.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)

    return np.array(froc_auc).astype('float')


def reliability_curve(true_binary,predicted_probabilities,cfreq,mlabels,nbins,pmax,ftag):

    predicted_probabilities = np.atleast_2d(predicted_probabilities)

    xd = yd = [0,1]
    clim_x = clim_y = [cfreq,cfreq]
    sk_line = [cfreq/2,(1-cfreq)/2+cfreq]

    fig1, ax = plt.subplots(figsize=(5.5,5.))
    ax.axhline(cfreq,color='grey',linestyle='-.',linewidth=1, zorder=2)
    ax.plot(clim_x,yd,color='grey',linestyle=':',linewidth=0.5)
    ax.text(0., cfreq, 'Clim', color='dimgrey', fontsize=6)
    ax.plot(xd,sk_line,color='grey',linestyle='--',linewidth=0.5)
    ax.fill_between(xd,xd,sk_line,facecolor='peachpuff',alpha=0.2)
    ax.fill_betweenx(yd,yd,clim_x,facecolor='peachpuff',alpha=0.2)

    for i in range(0,predicted_probabilities.shape[0]):
        ry, rx = calibration_curve(true_binary, predicted_probabilities[i,:], n_bins=nbins)
        # Max prob 
        indp = np.where(rx<pmax)[0]

        if i==0:
            ax.plot(rx[indp], ry[indp],color=pcolors[i],linestyle=lstl[i], label=mlabels[i], linewidth=2., zorder=4)
            ax.plot(rx[indp], ry[indp], marker='s', linewidth=1, color=pcolors[i], zorder=4)
        else:
            ax.plot(rx[indp], ry[indp],color=pcolors[i],linestyle=lstl[i], label=mlabels[i], linewidth=1.,alpha=0.7, zorder=3)
            # ax.plot(rx[indp], ry[indp], marker='.', linewidth=1, color=pcolors[i], zorder=3)

    line = mlines.Line2D([0, 1], [0, 1], color='dimgrey',linestyle='--', zorder=2)
    transform = ax.transAxes; line.set_transform(transform); ax.add_line(line)
    ax.set_xlabel('Forecast probability')
    ax.set_ylabel('Observed relative frequency')
    ax.legend(loc='lower right',fontsize=sl-3)
    plt.grid(c='silver', ls='--', alpha=0.2, zorder=1)
    plt.ylim(ymax = 1.005, ymin = -0.005)
    plt.xlim(xmax = 1.005, xmin = -0.005)
    # ax_inset = inset_axes(ax, width="30%", height="30%", loc="upper left")
    ax_inset = fig1.add_axes([0.16, 0.75, 0.25, 0.19])
    # Plot the histogram in the inset
    ax_inset.hist(predicted_probabilities[0,:], bins=np.arange(0, 1.1, 0.1), color='gray' , edgecolor='black', alpha=0.5, zorder=2)
    ax_inset.set_xticks(np.arange(0., 1.1, 0.1))
    ax_inset.set_xticklabels([f'{x:.1f}' for x in np.arange(0., 1.1, 0.1)])
    ax_inset.tick_params(axis='x', which='major', labelsize=5)
    ax_inset.set_yticks([]); ax_inset.set_yticklabels([])
    ax_inset.set_xlim(0,1)
    plt.tight_layout()
    # plt.show()
    # ax.set_title("Reliability diagram" )
    plt.savefig(ftag+"_ReliabilityDiagram.png", dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
        format='png', bbox_inches='tight', pad_inches=0.1)

    plt.close(fig1)

