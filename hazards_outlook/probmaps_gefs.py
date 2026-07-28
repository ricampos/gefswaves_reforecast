#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
probmaps_gefs.py

VERSION AND LAST UPDATE:
 v1.0  06/09/2023
 v1.1  01/06/2024
 v1.2  02/24/2024
 v1.3  07/31/2024
 v2.0  09/11/2024
 v2.1  12/17/2025
 v2.2  01/06/2025
 v2.3  07/22/2026

PURPOSE:
 This program makes the probability maps based on the NOAA Global Ensemble
  Forecast System (GEFS), marine forecast. It requires the grib2 files to
  have been previously downloaded (see download_GEFSwaves.sh)
 The global hazards outlook (probability map) shows the probability of 
  certain variable (wind speed, significant wave height, and peak period)
  to have at least one value above certain pre-defined level within the given
  time range.
 For example, selecting the week2 time range (from day 7 to day 14), it
  calculates the probabilities of encontering at least one instant above
  qlev (defined in probmaps_gefs.yaml) during this week (for each grid point).

USAGE:
 The information is passed to the script through 5 input arguments:
  1) .yaml configuration file probmaps_gefs.yaml;
  2) forecast cycle (YYYYMMDDHH) that will be used to compute the probabilities;
  3) initial day to define the time interval to calculate the statistics;
  4) final day to define the time interval to calculate the statistics;
  5) variable (select only one) to be processed: WS10 or Hs;
 The configuration probmaps_gefs.yaml contains the fixed information.
 This script must be run for each variable (WS10, Hs) separately.
 See the probmaps_gefs.yaml for more specific information and how to calibrate
  the probabilities and customize the plots.

 Example (from linux terminal command line):
  python3 probmaps_gefs.py probmaps_gefs_internal.yaml 2023060900 7 14 Hs
  nohup python3 probmaps_gefs.py probmaps_gefs_internal.yaml 2023060900 7 14 Hs >> nohup_probmaps_gefs.out 2>&1 &

OUTPUT:
 png and kmz figures (probability maps) saved in the outpath directory informed in
 the configuration file probmaps_gefs.yaml

DEPENDENCIES:
 See the imports below.

AUTHOR and DATE:
 06/09/2023: Ricardo M. Campos, first version.
 01/06/2024: Ricardo M. Campos, kmz format added and variable Tp removed. The date
  has been removed from the output file name (figure) so it can be replaced automatically
  aften each new run.
 02/24/2024: Ricardo M. Campos, include 5% prob level for hurricane force 
  winds and waves, in green.
 07/31/2024: Ricardo M. Campos, mask up / de-emphasize the tropics between 30S to 30N.
 09/11/2024: Ricardo M. Campos, mask the tropics only during hurricane season. Generalize
  the code to be working with any domain (lat/lon intervals).
 12/17/2025: Ricardo M. Campos, inclusion of KML with polygons and contours, able to open in 
  Google Earth.
 01/06/2025: Ricardo M. Campos, improved KML output.
 07/22/2026: Ricardo M. Campos, fixed contour_to_kml for matplotlib >=3.10, which
  removed QuadContourSet.collections. Each level band is now contoured with its
  own contourf() call instead of indexing into a combined ContourSet, so the
  code works across matplotlib versions.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import pygrib
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from matplotlib.colors import ListedColormap
from scipy.ndimage import gaussian_filter
from shapely.geometry import Polygon
import simplekml
from matplotlib import colors
from shapely.ops import unary_union
import cartopy
import cartopy.crs as ccrs
# Palette and colors for plotting the figures
import sys
import warnings; warnings.filterwarnings("ignore")
# --------------------------------------------------------------------------
sl=13 # plot style configuration
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


def contour_to_kml(lon, lat, Z, levels, lname, tstart, tend, pcolors, filename):

    def mpl_color_to_kml(color, alpha):
        r, g, b = colors.to_rgb(color)
        a = int(alpha * 255)
        return f"{a:02x}{int(b*255):02x}{int(g*255):02x}{int(r*255):02x}"

    def ensure_counter_clockwise(coords):
        if len(coords) < 4:
            return coords

        area = sum((coords[i+1][0] - coords[i][0]) * (coords[i+1][1] + coords[i][1]) 
                  for i in range(len(coords)-1))
        return coords if area < 0 else coords[::-1]

    def ensure_clockwise(coords):
        if len(coords) < 4:
            return coords

        area = sum((coords[i+1][0] - coords[i][0]) * (coords[i+1][1] + coords[i][1]) 
                  for i in range(len(coords)-1))
        return coords if area > 0 else coords[::-1]

    kml = simplekml.Kml()
    kml.document.name = f"Probability {lname}, {tstart} to {tend}"

    # --- 1. DEFINE SHARED STYLES ---
    styles_lookup = {}
    for i, lvl in enumerate(levels):
        shared_style = simplekml.Style()
        shared_style.linestyle.color = mpl_color_to_kml(pcolors[i], 1.0)
        shared_style.linestyle.width = 3
        shared_style.polystyle.color = mpl_color_to_kml(pcolors[i], 0.3)
        shared_style.polystyle.fill = 0
        shared_style.polystyle.outline = 1
        shared_style.polystyle.tessellate = 0
        shared_style.balloonstyle.text = (
            f"<![CDATA[<b>Probability:</b> P ≥ {lvl:.2f}<br>"
            f"<b>Date Range:</b> {tstart} to {tend}]]>")
        styles_lookup[lvl] = shared_style

    # --- 2. DOMAIN BACKGROUND ---
    bg_style = simplekml.Style()
    bg_style.polystyle.color = '4c696969'
    bg_style.polystyle.outline = 0
    lon_min, lon_max = lon.min(), lon.max()
    lat_min, lat_max = lat.min(), lat.max()
    bg = kml.newpolygon(
        name="Model domain",
        outerboundaryis=[(round(lon_min,2), round(lat_min,2), 0),(round(lon_max,2), round(lat_min,2), 0),
                         (round(lon_max,2), round(lat_max,2), 0),(round(lon_min,2), round(lat_max,2), 0),(round(lon_min,2), round(lat_min,2), 0)])
    bg.style = bg_style
    bg.altitudemode = simplekml.AltitudeMode.absolute

    # --- 3. CONTOURS ---
    # matplotlib >=3.10 removed QuadContourSet.collections (previously one
    # contourf() call returned one PathCollection per level band, indexed
    # via csf.collections[i]). To stay compatible across matplotlib
    # versions, each level band is now contoured with its own separate
    # contourf() call, so there's no need to index into a combined,
    # multi-level ContourSet at all.
    for i in range(len(levels) - 1):
        lvl_lower = levels[i]
        lvl_upper = levels[i + 1]
        folder = kml.newfolder(name=f"P {lvl_lower:.2f}–{lvl_upper:.2f}")

        csf = plt.contourf(lon, lat, Z, levels=[lvl_lower, lvl_upper])
        plt.close()

        # matplotlib < 3.8: ContourSet.collections[0].get_paths()
        # matplotlib >= 3.10: ContourSet.get_paths() directly (ContourSet
        # itself is now a single Collection)
        if hasattr(csf, 'collections'):
            paths = csf.collections[0].get_paths()
        else:
            paths = csf.get_paths()

        polys = []
        for path in paths:
            for poly_coords in path.to_polygons():
                if len(poly_coords) < 3:
                    continue
                polys.append(Polygon(poly_coords))

        if not polys:
            continue

        merged = unary_union(polys)
        geoms = [merged] if merged.geom_type == "Polygon" else list(merged.geoms)

        for k, g in enumerate(geoms):
            # Fix winding order
            exterior_coords = [(round(x,2), round(y,2), 0) for x,y in g.exterior.coords]
            exterior_coords = ensure_counter_clockwise(exterior_coords)
            
            interiors_coords = []
            for ring in g.interiors:
                inner = [(round(x,2), round(y,2), 0) for x,y in ring.coords]
                interiors_coords.append(ensure_clockwise(inner))

            pol = folder.newpolygon( name=f"{lname} | P {lvl_lower:.2f}–{lvl_upper:.2f} | #{k+1}",
                outerboundaryis=exterior_coords)

            for inner in interiors_coords:
                pol.innerboundaryis.append(inner)

            pol.style = styles_lookup[lvl_lower]
            pol.style.polystyle.fill = 0
            pol.draworder = int(1+lvl_lower * 100)
            pol.altitudemode = simplekml.AltitudeMode.absolute
            pol.extrude = 0

    # --- 4. SAVE KML ---
    kml.document.timespan.begin = tstart
    kml.document.timespan.end = tend
    kml.save(filename)


if __name__ == "__main__":
    # Input Arguments -----
    # .yaml configuration file
    fconfig=str(sys.argv[1])
    # Forecast Cycle
    fcycle=str(sys.argv[2])
    fcdate=str(fcycle[0:8]); fchour=str(fcycle[8:10])
    # Forecast Lead Time (Day) and intervall
    ltime1=int(sys.argv[3])
    ltime2=int(sys.argv[4])
    # forecast variable: WS10, Hs
    fvarname=str(sys.argv[5])

    # Fixed configuration variables, read yaml file -----------
    print(" "); print(" Reading yaml configuration file ...")
    with open(fconfig, 'r') as file:
        wconfig = yaml.safe_load(file)

    ftag=str(wconfig['ftag'])

    # plot mode (internal/external)
    mode=str(wconfig['mode'])
    # number of ensemble members
    nenm=wconfig['nenm']
    # time resolution
    tres=wconfig['tres']
    # n-max expansion
    nmax=wconfig['nmax']
    # spatial window size (diameter), square (degrees)
    spws=wconfig['spws']
    # Plotting Area
    slonmin=wconfig['lonmin']-spws; slonmax=wconfig['lonmax']+spws
    slatmin=wconfig['latmin']-spws; slatmax=wconfig['latmax']+spws
    # Spatial percentile (grid points)
    spctl=wconfig['spctl']
    # gaussian filter spatial plot smoothing
    gft=wconfig['gft']
    # percentiles for the initial plots
    pctls=np.array(wconfig['pctls']).astype('int')
    # Probability levels
    plevels = np.array(wconfig['plevels']).astype('float')
    hplevels = np.array(wconfig['hplevels']).astype('float')
    # Colors for the probability maps
    pcolors = np.array(wconfig['pcolors']).astype('str')
    hpcolors = np.array(wconfig['hpcolors']).astype('str')

    # Path of GEFSv12 grib2 files
    gefspath=wconfig['gefspath']
    if gefspath[-1] != '/':
        gefspath=gefspath+"/"

    # output path
    outpath=str(wconfig['outpath'])
    if outpath[-1] != '/':
        outpath=outpath+"/"

    umf=1. # unit conversion, when necessary

    # maximum value allowed (quick quality control), variable name (grib2), and levels for the probability plot
    if fvarname.upper() == "WS10" or fvarname.upper() == "WND" or fvarname.upper() == "U10":
        qqvmax=wconfig['qqvmax_wnd']
        mvar=wconfig['mvar_wnd']
        qlev=np.array(wconfig['qlev_wnd']).astype('float')
        qlev_fnames=np.array(wconfig['qlev_wnd_fnames']).astype('str')
        vtickd=int(wconfig['vtickd_wnd'])
        funits=str('knots')
        umf=1.94 # m/s to knots
    elif fvarname.upper() == "HS":
        qqvmax=wconfig['qqvmax_hs']
        mvar=wconfig['mvar_hs']
        funits=str('m')
        qlev=np.array(wconfig['qlev_hs']).astype('float')
        qlev_fnames=np.array(wconfig['qlev_hs_fnames']).astype('str')
        vtickd=int(wconfig['vtickd_hs'])
    else:
        sys.exit(" Input variable "+fvarname+" not included in the list. Please select only one: WS10, Hs.")

    print(" Reading yaml configuration file, OK."); print(" ")

    # Time range forecast intervall string, for the plots
    if (ltime2-ltime1)<0:
        aux=np.copy(ltime1); ltime1=np.copy(ltime2)
        ltime2=np.copy(aux); del aux

    if ltime1==1 and (ltime2-ltime1)==7:
        trfi=str("Week 1 - ")
    elif ltime1==7 and (ltime2-ltime1)==7:
        trfi=str("Week 2 - ")
    elif ltime1==14 and (ltime2-ltime1)==7:
        trfi=str("Week 3 - ")
    elif ltime1==21 and (ltime2-ltime1)==7:
        trfi=str("Week 4 - ")
    elif ltime1==28 and (ltime2-ltime1)==7:
        trfi=str("Week 5 - ")
    elif (ltime2-ltime1)>0:
        trfi=str("Days "+str(ltime1)+"-"+str(ltime2)+" , ")
    elif ltime2==ltime1:
        trfi=str("Day "+str(ltime1)+" , ")
    else:
        trfi=''

    print(" "); print(" 1. Reading Forecast Data ...")

    auxltime = np.arange((ltime1-1)*24,((ltime2)*24)+1,tres)
    auxltime[auxltime>384]=384; auxltime[auxltime<0]=0

    # READ WW3 Ensemble Forecast files. Appending forecast days (time intervall).
    c=0
    for t in range(0,auxltime.shape[0]):
        for enm in range(0,nenm):

            # if enm==0:
            #    fname=gefspath+"gefs."+fcdate+"/"+fchour+"/wave/gridded/gefs.wave.t"+fchour+"z.c"+str(enm).zfill(2)+".global.0p25.f"+str(auxltime[t]).zfill(3)+".grib2"
            # else:
            #    fname=gefspath+"gefs."+fcdate+"/"+fchour+"/wave/gridded/gefs.wave.t"+fchour+"z.p"+str(enm).zfill(2)+".global.0p25.f"+str(auxltime[t]).zfill(3)+".grib2"
            fname=gefspath+"GEFSv12Waves_"+fcdate+fchour+"/gefs.wave."+fcdate+"."+str(enm).zfill(2)+".global.0p25.f"+str(auxltime[t]).zfill(3)+".grib2"

            if c==0:
                ds = xr.open_dataset(fname, engine='cfgrib')
                wtime = np.atleast_1d(np.array(ds.time.values))
                lat = np.array(ds.latitude.values); lat = np.sort(lat); lon = np.array(ds.longitude.values)
                fmod=np.zeros((auxltime.shape[0],nenm,lat.shape[0],lon.shape[0]),'f')*np.nan
                ds.close(); del ds

            # opening and reading with pygrib (faster than xarray cfgrib)
            grbs = pygrib.open(fname); glist=list(grbs.select())
            ind=np.nan
            for i in range(0,len(glist)):
                if mvar in str(glist[i]):
                    ind=int(i)
            if ind>=0:
                grb = grbs.select()[ind]
                fmod[t,enm,:,:]=np.flip(grb.values,axis=0)
            else:
                sys.exit(" Enviromental variable "+mvar+" not found in the grib2 file "+fname)

            c=c+1; del ind

        print(repr(t))

    # Quick simple quality control
    fmod[fmod>=qqvmax]=np.nan; fmod[fmod<0.]=np.nan
    # Unit conversion
    fmod=fmod*umf
    # Select domain of interest
    indlat=np.where((lat>=(slatmin-spws))&(lat<=(slatmax+spws)))
    indlon=np.where((lon>=(slonmin-spws))&(lon<=(slonmax+spws)))
    if np.size(indlon)>0 and np.size(indlat)>0:
        lat=np.copy(lat[indlat[0]]); lon=np.copy(lon[indlon[0]])
        fmod=np.copy(fmod[:,:,indlat[0],:][:,:,:,indlon[0]])
    else:
        sys.exit(" Min/Max lat and lon incorrect when applied to the forecast file. Check longitude standards.")

    del indlat, indlon
    print(" 1. Forecast Data ... OK"); print(" ")

    print(" 2. Space-Time Cells and Probabilities ...")

    # n-max expansion and reshape
    fmod=np.sort(fmod,axis=0)
    fmod=np.copy(fmod[-nmax::,:,:,:])
    fmod=np.array(fmod.reshape(nmax*nenm,lat.shape[0],lon.shape[0]))

    gspws=int(np.floor(spws/np.diff(lat).mean())/2)

    probecdf=np.zeros((qlev.shape[0],lat.shape[0],lon.shape[0]),'f')
    for i in range(0,qlev.shape[0]):
        for j in range(0,lat.shape[0]):
            for k in range(0,lon.shape[0]):
                if np.any(fmod[:,j,k]>0.0):
                    if (j>=gspws) and (j<=lat.shape[0]-gspws) and (k>=gspws) and (k<=lon.shape[0]-gspws):
                        aux=np.array(fmod[:,(j-gspws):(j+gspws+1),:][:,:,(k-gspws):(k+gspws+1)])
                        aux=aux.reshape(nmax*nenm,aux.shape[1]*aux.shape[2])
                        aux=np.sort(aux,axis=1)
                        ind=np.where(np.mean(aux,axis=0)>=0.)
                        if np.size(ind)>0:
                            aux=np.array(aux[:,ind[0]])            
                            aux=np.array(aux[:,int(np.floor(aux.shape[1]*(spctl/100)))::])
                            aux=aux.reshape(aux.shape[0]*aux.shape[1])
                            probecdf[i,j,k] = np.size(aux[aux>qlev[i]]) / np.size(aux[aux>0.])

                        del aux

    print(" 2. Space-Time Cells and Probabilities ... OK"); print(" ")

    # PLOTS
    print(" 3. Probability Maps ...")

    # Probability levels and labels for the plots
    clabels=[];clevels=[]
    for j in range(0,np.size(plevels)-1):
        clabels=np.append(clabels,">"+str(int(plevels[j]*100)).zfill(2)+"%")
        clevels=np.append(clevels,(plevels[j]+plevels[j+1])/2)

    cmap = ListedColormap(pcolors)

    # If hurricane season (May 15th to Nov 31st), it masks (or cuts) the tropics for external plots
    mmdd = wtime[0].astype('datetime64[D]').astype(object).month * 100 + wtime[0].astype('datetime64[D]').astype(object).day
    if mode=='external' and (mmdd>=515 and mmdd<=1130):
        slatmin=wconfig['latminhs']-spws # min latitude for the hurricane season

    # Loop through the intensity levels
    for i in range(0,qlev.shape[0]):

        # Extra percentage level associated with the most extreme case
        if i == int(qlev.shape[0]-1):

            plevels = hplevels
            pcolors = hpcolors

            clabels=[];clevels=[]
            for j in range(0,np.size(plevels)-1):
                clabels=np.append(clabels,">"+str(int(plevels[j]*100))+"%")
                clevels=np.append(clevels,(plevels[j]+plevels[j+1])/2)

            cmap = ListedColormap(pcolors)

        # Figure
        plt.figure(figsize=(9,5.5))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
        ax.set_extent([slonmin+spws,slonmax-spws,slatmin+spws,slatmax-spws], crs=ccrs.PlateCarree())
        # ax.set_extent([slonmin,slonmax,slatmin,slatmax], crs=ccrs.PlateCarree())
        gl = ax.gridlines(crs=ccrs.PlateCarree(),  xlocs=range(-180,180, 20), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
        gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
        ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
        ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
        ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
        ax.coastlines(resolution='50m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)
        title = "Prob "+fvarname+">"+str(qlev[i]).zfill(1)+funits+", Cycle "+fcycle[0:8]+" "+fcycle[8:10]+"Z \n"
        title += r"$\bf{"+trfi+"Valid: "+pd.to_datetime(wtime[0]+np.timedelta64(ltime1,'D')).strftime('%B %d, %Y')+" - "
        title += pd.to_datetime(wtime[0]+np.timedelta64(ltime2,'D')).strftime('%B %d, %Y')+"}$"

        # If external and in the hurricane season, it masks (or cuts) the tropics.
        if mode=='external' and (mmdd>=515 and mmdd<=1130) and np.any((np.linspace(wconfig['latminhs'],slatmax,10)>-25.) & (np.linspace(wconfig['latminhs'],slatmax,10)<25.)):
            contour_data = gaussian_filter(probecdf[i, :, :], gft)
            tropical_mask = (lat >= -30) & (lat <= 30)
            # Expand tropical_mask to match the shape of contour_data
            tropical_mask_expanded = np.tile(tropical_mask[:, np.newaxis], (1, lon.shape[0]))
            # Add transparency to tropical areas
            cs = ax.contourf(lon, lat, np.ma.masked_where(tropical_mask_expanded, contour_data), 
                                          levels=plevels, cmap=cmap, alpha=0.7, zorder=1, transform=ccrs.PlateCarree())

            ax.add_patch(plt.Rectangle((slonmin, -30), 360, 60, transform=ccrs.PlateCarree(), color='whitesmoke', alpha=0.4, zorder=2, hatch='//', edgecolor='none'))
            ax.add_patch(plt.Rectangle((slonmin, -30), 360, 60, transform=ccrs.PlateCarree(), color='whitesmoke', alpha=0.4, zorder=2, hatch='\\', edgecolor='none'))

            plt.text(slonmin+spws-270+1, 0., 'This product is not intended for tropical cyclones', color='k', fontsize=8, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'), zorder=4)
            plt.text(slonmin+spws-270+1, -5., 'Tropical latitudes are masked during hurricane season.', color='k', fontsize=7, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),zorder=4)

        else:
            cs=ax.contourf(lon,lat,gaussian_filter(probecdf[i,:,:],gft),levels=plevels,alpha=0.7,cmap=cmap,zorder=1,transform = ccrs.PlateCarree())

        ax.set_title(title); del title
        plt.tight_layout()
        ax = plt.gca(); pos = ax.get_position(); l, b, w, h = pos.bounds; cax = plt.axes([l+0.06, b-0.07, w-0.12, 0.03]) # setup colorbar axes.
        cbar = plt.colorbar(cs,cax=cax, orientation='horizontal',ticks=clevels, format='%g')
        cbar.ax.set_xticklabels(clabels)
        cbar.ax.tick_params(length=0)
        for label in cbar.ax.get_xticklabels():
            label.set_weight('bold')

        plt.axes(ax); plt.tight_layout()
        plt.text(-90., 76., 'Experimental', color='k', fontsize=13, fontweight='bold')
        figname = outpath+"ProbMap_"+fvarname+"_"+str(qlev[i]).zfill(1)+"_fcst"+str(ltime1).zfill(2)+"to"+str(ltime2).zfill(2)+"_"+mode+"_"+ftag
        plt.savefig(figname+".png", dpi=200, facecolor='w', edgecolor='w',
                orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

        plt.close('all'); del ax

        # KML output file ----------------
        tstart = pd.to_datetime(wtime[0] + np.timedelta64(ltime1, 'D')).strftime('%Y/%m/%d')
        tend = pd.to_datetime(wtime[0] + np.timedelta64(ltime2, 'D')).strftime('%Y/%m/%d')

        contour_to_kml(lon, lat, gaussian_filter(probecdf[i, :, :], gft),levels=plevels,
            lname=qlev_fnames[i],tstart=tstart,tend=tend,pcolors=np.append(pcolors, pcolors[-1]),
            filename=figname + ".kml")

        del figname
        # --------------

        print("   Plot ... qlev "+repr(qlev[i]))

    print(" 3. Probability Plots ... OK"); print(" ")
