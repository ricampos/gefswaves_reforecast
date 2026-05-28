#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
combine_buoy_alt.py

VERSION AND LAST UPDATE:
 v1.0  06/25/2025
 v1.1  05/27/2026

PURPOSE:
 Part of the fuzzy verification process and dataset. This script includes NDBC buoy data to points (when available)
  in Altimeter.PointExtract*.nc, to build a complete set of observations.

USAGE:
 Edit paths and file names below.

DEPENDENCIES:
 See the imports below.

AUTHOR and DATE:
 06/25/2025: Ricardo M. Campos, first version, named fuzzy_verification_ProbMaps.py
 05/27/2026: Ricardo M. Campos, optimize the code to run faster

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import netCDF4 as nc
import xarray as xr
import numpy as np
import warnings; warnings.filterwarnings("ignore")
# netcdf format
fnetcdf="NETCDF4"

# Observations NDBC, netcdf format
def tseriesnc_ndbc(fname=None,anh=None):
    '''
    Observations NDBC, time series/table, netcdf format
    Input: file name (example: 46047h2016.nc), and anemometer height (optional)
    Output: dictionary containing the arrays: time(seconds since 1970),time(datetime64),lat,lon, 
      and arrays sst,mslp,dwp,tmp,gst(10-m height),wsp(10-m height),wdir,hs,tm,tp,dm
    '''
    if fname==None:
        raise ValueError("NDBC file name must be informed.")

    try:
        ds = xr.open_dataset(fname); f=nc.Dataset(fname)
    except:
        sys.exit(" Cannot open "+fname)
    else:
        btp = f.variables['dominant_wpd'][:,0,0]
        btime = np.array(f.variables['time'][:]).astype('double')
        f.close(); del f 

        if 'wind_spd' in ds.keys():
            bwsp = ds['wind_spd'].values[:,0,0]
            if anh==None:
                anh=4.1 # NDBC (3M buoys)

            # convert wind speed to 10 meters (DNVGL C-205 Table 2-1, confirmed by https://onlinelibrary.wiley.com/doi/pdf/10.1002/er.6382)
            bwsp =  np.copy(((10./anh)**(0.12)) * bwsp)

        bhs = ds['wave_height'].values[:,0,0]

        # Automatic and basic Quality Control
        bwsp[(bwsp<0)|(bwsp>100)]=np.nan
        bhs[(bhs<0)|(bhs>20)]=np.nan
        btp[(btp<0)|(btp>30)]=np.nan

        result={'latitude':np.array(ds['latitude'].values[:]),'longitude':np.array(ds['longitude'].values[:]),
        'time':btime,'wind_spd':bwsp,'hs':bhs,'tp':btp}

        return result
        ds.close()
        del ds,btime,bwsp,bhs,btp


if __name__ == "__main__":

    # Altimeter data, from extract_altimeter.py
    fname='/work/noaa/marine/ricardo.campos/work/analysis/3assessments/fuzzy_verification/data/Atlantic/Altimeter.PointExtract.Atlantic_20201001to20260101.nc'
    # Paths
    ndbcp="/work/noaa/marine/ricardo.campos/data/buoys/NDBC/ncformat/wparam"

    ds = xr.open_dataset(fname)
    mdate = ds['time'].values[:]
    f=nc.Dataset(fname)
    pid = f.variables['ID'][:]
    mtime = f.variables['time'][:]

    blat = f.variables['lat'][:]; blon = f.variables['lon'][:]
    mlat = f.variables['mlat'][:]; mlon = f.variables['mlon'][:]
    u10m = f.variables['u10_mean'][:]; hsm = f.variables['hs_mean'][:]
    u10p80 = f.variables['u10_p80'][:]; hsp80 = f.variables['hs_p80'][:]
    u10p90 = f.variables['u10_p90'][:]; hsp90 = f.variables['hs_p90'][:]
    u10p95 = f.variables['u10_p95'][:]; hsp95 = f.variables['hs_p95'][:]
    u10max = f.variables['u10_max'][:]; hsmax = f.variables['hs_max'][:]

    bhs=np.zeros((len(pid),len(mtime)),'f')*np.nan
    btp=np.zeros((len(pid),len(mtime)),'f')*np.nan
    bu10=np.zeros((len(pid),len(mtime)),'f')*np.nan

    for i in range(0,len(pid)):
        print(" Start "+pid[i])

        # Read NDBC, if available
        years = np.unique(mdate.astype('datetime64[Y]').astype(int) + 1970)
        ahs=[]; atp=[]; au10=[]; btime=[]
        for j in range(0,len(years)):

            try:
                res = tseriesnc_ndbc(fname=ndbcp+"/"+pid[i]+"h"+str(years[j])+".nc",anh=None)
                ahs = np.append(ahs,res['hs'])
                atp = np.append(atp,res['tp'])
                au10 = np.append(au10,res['wind_spd'])
                btime = np.append(btime,res['time'])
                del res
            except:
                print(" Cannot open "+pid[i]+"h"+str(years[j]))
            else:
                print(" Ok "+pid[i]+"h"+str(years[j]))

        if np.any(ahs):
            c=0
            for t in range(0,len(mtime)):
                indt=np.where(np.abs(btime-mtime[t])<1800.)
                if np.size(indt)>0:
                    bhs[i,t] = float(np.nanmean(ahs[indt[0]]))
                    bu10[i,t] = float(np.nanmean(au10[indt[0]]))
                    btp[i,t] = float(np.nanmean(atp[indt[0]]))
                    c=c+1
                    del indt

            print(" Ok "+pid[i]+", with "+repr(c)+" buoy records")

        print(" Done "+pid[i]); print(" ")


    print(' ')
    # Save netcdf output file 
    ncfile = nc.Dataset("Altimeter.Buoy.PointExtract_20201001to20260101.nc", "w", format=fnetcdf)
    ncfile.history="AODN Altimeter data extracted for fixed point outputs, and NDBC buoy data for those points (when available)"
    # create  dimensions
    ncfile.createDimension('points', pid.shape[0] )
    ncfile.createDimension('time', len(mtime) )
    ncfile.createDimension('grid', int(mlat.shape[1]) )
    # create variables.
    vt = ncfile.createVariable('time',np.dtype('float64').char,('time'))
    vbid = ncfile.createVariable('ID',np.dtype('a25'),('points'))
    vblat = ncfile.createVariable('lat',np.dtype('float32').char,('points'))
    vblon = ncfile.createVariable('lon',np.dtype('float32').char,('points'))
    vmlat = ncfile.createVariable('mlat',np.dtype('float32').char,('points','grid'))
    vmlon = ncfile.createVariable('mlon',np.dtype('float32').char,('points','grid'))
    # results
    vu10m = ncfile.createVariable('alt_u10_mean',np.dtype('float32').char,('points','time','grid'))
    vhsm = ncfile.createVariable('alt_hs_mean',np.dtype('float32').char,('points','time','grid'))
    vu10p80 = ncfile.createVariable('alt_u10_p80',np.dtype('float32').char,('points','time','grid'))
    vhsp80 = ncfile.createVariable('alt_hs_p80',np.dtype('float32').char,('points','time','grid'))
    vu10p90 = ncfile.createVariable('alt_u10_p90',np.dtype('float32').char,('points','time','grid'))
    vhsp90 = ncfile.createVariable('alt_hs_p90',np.dtype('float32').char,('points','time','grid'))
    vu10p95 = ncfile.createVariable('alt_u10_p95',np.dtype('float32').char,('points','time','grid'))
    vhsp95 = ncfile.createVariable('alt_hs_p95',np.dtype('float32').char,('points','time','grid'))
    vu10max = ncfile.createVariable('alt_u10_max',np.dtype('float32').char,('points','time','grid'))
    vhsmax = ncfile.createVariable('alt_hs_max',np.dtype('float32').char,('points','time','grid'))
    #
    vbhs = ncfile.createVariable('buoy_hs',np.dtype('float32').char,('points','time'))
    vbu10 = ncfile.createVariable('buoy_u10',np.dtype('float32').char,('points','time'))
    vbtp = ncfile.createVariable('buoy_tp',np.dtype('float32').char,('points','time'))
    # Assign units
    vt.units = 'seconds since 1970-01-01T00:00:00+00:00'
    vblat.units = 'degrees_north' ; vblon.units = 'degrees_east'
    vu10m.units = 'm/s'; vhsm.units = 'm'
    vu10p80.units = 'm/s'; vhsp80.units = 'm'
    vu10p90.units = 'm/s'; vhsp90.units = 'm'
    vu10p95.units = 'm/s'; vhsp95.units = 'm'
    vu10max.units = 'm/s'; vhsmax.units = 'm'
    vbhs.units = 'm'; vbu10.units = 'm/s'; vbtp.units = 's'
    # Allocate Data
    vt[:]=mtime[:]; vbid[:]=pid[:]
    vblat[:]=np.array(blat[:]).astype('float'); vblon[:]=np.array(blon[:]).astype('float')
    vmlat[:,:]=np.array(mlat[:,:]).astype('float'); vmlon[:,:]=np.array(mlon[:,:]).astype('float')
    vu10m[:,:,:]=u10m[:,:,:]; vhsm[:,:,:]=hsm[:,:,:]
    vu10p80[:,:,:]=u10p80[:,:,:]; vhsp80[:,:,:]=hsp80[:,:,:]
    vu10p90[:,:,:]=u10p90[:,:,:]; vhsp90[:,:,:]=hsp90[:,:,:]
    vu10p95[:,:,:]=u10p95[:,:,:]; vhsp95[:,:,:]=hsp95[:,:,:]
    vu10max[:,:,:]=u10max[:,:,:]; vhsmax[:,:,:]=hsmax[:,:,:]
    # 
    vbhs[:,:]=bhs[:,:]; vbu10[:,:]=bu10[:,:]; vbtp[:,:]=btp[:,:]
    # 
    ncfile.close()
    print("Done. Netcdf ok. New file saved: Altimeter.Buoy.PointExtract_20201001to20260101.nc")


