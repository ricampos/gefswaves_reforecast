#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
plot_positions.py

VERSION AND LAST UPDATE:
 v1.0  05/13/2025

PURPOSE:
 Small script to plot the positions used for validation
  selected.

"""

import matplotlib
# matplotlib.use('Agg')
import pandas as pd
from matplotlib import ticker
from matplotlib.patches import Circle
from pylab import *
import cartopy
import cartopy.crs as ccrs
import warnings; warnings.filterwarnings("ignore")

sl=14
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

if __name__ == "__main__":

    fname = 'points_sel_Pacific.txt'
    ds_points = pd.read_csv(fname,comment='#',delimiter=r"\s+")

    # ==== Plot point positions ====
    plt.close('all')
    plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([125,300,-50,65.], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("lightgrey"), edgecolor='grey',linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='50m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)

    for i in range(0,len(ds_points)):
        bid=ds_points['ID'][i]; blat=float(ds_points['Lat'][i]); blon=float(ds_points['Lon'][i])
        ax.scatter(blon,blat,s=9, c='b', marker='o',transform=ccrs.PlateCarree(),zorder=3)

    ind=np.array(ds_points[ds_points['ID'].str.startswith('4')].index).astype('int')
    for i in range(0,len(ind)):
        bid=ds_points['ID'][ind[i]]; blat=float(ds_points['Lat'][ind[i]]); blon=float(ds_points['Lon'][ind[i]])
        ax.scatter(blon,blat,s=9, c='g', marker='o',transform=ccrs.PlateCarree(),zorder=3)

    ind=np.array(ds_points[ds_points['ID'].str.startswith('5')].index).astype('int')
    for i in range(0,len(ind)):
        bid=ds_points['ID'][ind[i]]; blat=float(ds_points['Lat'][ind[i]]); blon=float(ds_points['Lon'][ind[i]])
        ax.scatter(blon,blat,s=9, c='g', marker='o',transform=ccrs.PlateCarree(),zorder=3)

    plt.savefig('Point_Positions.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close()



    # --------------------------------------------------
    # Ilustration of grid points and altimeter selection
    alat=np.arange(27.,28.+0.001,0.25)
    alon=np.arange(-84.5,-83.5+0.001,0.25)
    plat=np.arange(27.,28.+0.001,0.5)
    plon=np.arange(-84.5,-83.5+0.001,0.5)
    plat2=np.array([27.25,27.75])
    plon2=np.array([-84.25,-83.75])
    plt.close('all')
    plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([-85.,-82.2,26.2,28.8], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("antiquewhite"), edgecolor='grey', alpha=0.5,linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='10m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)

    for i in range(0,len(alat)):
        for j in range(0,len(alon)):
            ax.scatter(alon[j],alat[i],s=20, c='k', marker='.',transform=ccrs.PlateCarree(),zorder=3)

    ax.scatter(alon[2],alat[2],s=80,linewidths=3, c='k', marker='x',transform=ccrs.PlateCarree(),zorder=5)
    plt.savefig('Illustration_1.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close()


    # 
    plt.close('all')
    plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([-85.,-82.2,26.2,28.8], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("antiquewhite"), edgecolor='grey', alpha=0.5,linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='10m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)

    for i in range(0,len(alat)):
        for j in range(0,len(alon)):
            ax.scatter(alon[j],alat[i],s=10, c='dimgrey', marker='.',transform=ccrs.PlateCarree(),zorder=3)

    for i in range(0,len(plat)):
        for j in range(0,len(plon)):
            ax.scatter(plon[j],plat[i],s=20, c='red', marker='o',transform=ccrs.PlateCarree(),zorder=3)

    for i in range(0,len(plat2)):
        for j in range(0,len(plon2)):
            ax.scatter(plon2[j],plat2[i],s=20, c='red', marker='o',transform=ccrs.PlateCarree(),zorder=3)

    plt.savefig('Illustration2.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close()


    # 
    plt.close('all')
    plt.figure(figsize=(7,5))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
    ax.set_extent([-85.,-82.2,26.2,28.8], crs=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
    gl.xlabel_style = {'size': 9, 'color': 'k','rotation':0}; gl.ylabel_style = {'size': 9, 'color': 'k','rotation':0}
    ax.add_feature(cartopy.feature.OCEAN,facecolor=("white"))
    ax.add_feature(cartopy.feature.LAND,facecolor=("antiquewhite"), edgecolor='grey', alpha=0.5,linewidth=0.5, zorder=2)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='grey', linestyle='-',linewidth=0.5, alpha=1, zorder=3)
    ax.coastlines(resolution='10m', color='dimgrey',linewidth=0.5, linestyle='-', alpha=1, zorder=4)

    for i in range(0,len(alat)):
        for j in range(0,len(alon)):
            ax.scatter(alon[j],alat[i],s=10, c='dimgrey', marker='.',transform=ccrs.PlateCarree(),zorder=3)

    for i in range(0,len(plat)):
        for j in range(0,len(plon)):
            ax.scatter(plon[j],plat[i],s=20, c='red', marker='o',transform=ccrs.PlateCarree(),zorder=3)

            circle = Circle(
                (plon[j], plat[i]),            # (lon, lat)
                radius=0.25,                   # ~0.25 degrees (~25 km)
                facecolor='none',              # or 'blue' if you want filled
                edgecolor='blue',
                linewidth=1.,
                transform=ccrs.PlateCarree(),  # Important!
                zorder=3
            )
            ax.add_patch(circle)

            circle = Circle(
                (plon[j], plat[i]),           # (lon, lat)
                radius=0.25,                  # Approx 25 km in degrees
                facecolor='lightblue',             # Fill color
                edgecolor='none',             # No border (or set to 'blue')
                alpha=0.3,                    # Transparency level (0=transparent, 1=opaque)
                transform=ccrs.PlateCarree(), # Respect map projection
                zorder=2                      # Plot beneath the points
            )
            ax.add_patch(circle)

    for i in range(0,len(plat2)):
        for j in range(0,len(plon2)):
            ax.scatter(plon2[j],plat2[i],s=20, c='red', marker='o',transform=ccrs.PlateCarree(),zorder=3)

            circle = Circle(
                (plon2[j], plat2[i]),            # (lon, lat)
                radius=0.25,                   # ~0.25 degrees (~25 km)
                facecolor='none',              # or 'blue' if you want filled
                edgecolor='blue',
                linewidth=1.,
                transform=ccrs.PlateCarree(),  # Important!
                zorder=3
            )
            ax.add_patch(circle)

            circle = Circle(
                (plon2[j], plat2[i]),           # (lon, lat)
                radius=0.25,                  # Approx 25 km in degrees
                facecolor='lightblue',             # Fill color
                edgecolor='none',             # No border (or set to 'blue')
                alpha=0.3,                    # Transparency level (0=transparent, 1=opaque)
                transform=ccrs.PlateCarree(), # Respect map projection
                zorder=2                      # Plot beneath the points
            )
            ax.add_patch(circle)


    plt.savefig('Illustration3.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

    ax.scatter(alon[2],alat[2],s=80,linewidths=3, c='k', marker='x',transform=ccrs.PlateCarree(),zorder=5)

    plt.savefig('Illustration4.png', dpi=300, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)


