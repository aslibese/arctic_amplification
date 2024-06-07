#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon May 20 2024

This script creates global maps of the warming contributions (K) of individual feedback mechanisms, ERF, OHT and AHT 
during 1979-2023.

@author: aslibese
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from matplotlib.colors import BoundaryNorm
import warnings
warnings.filterwarnings("ignore")
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import os

# load the warming contributions dataset
ds_path = 'data/CMIP6/CanESM5_r1i1p1f1_warming_contr.nc'
warming_contr_ds = xr.open_dataset(ds_path)

# extract model and ensemble name from the path
model = os.path.basename(ds_path).split('_')[0]
ensemble = os.path.basename(ds_path).split('_')[1]

# extract the variable names 
variables = [warming_contr_ds[var].name for var in warming_contr_ds.data_vars]
long_names = ['a) Lapse Rate', 'b) Planck', 'c) Water Vapour$_{LW}$', 'd) Water Vapour$_{SW}$', 'e) Albedo', 'f) Cloud$_{LW}$', 
              'g) Cloud$_{SW}$', 'h) ERF$_{LW}$', 'i) ERF$_{SW}$', 'j) OHU', 'k) AHT']

# set up plot with 11 subplots
fig, axes = plt.subplots(4, 3, figsize=(20, 20), subplot_kw={'projection': ccrs.PlateCarree()})
axes = axes.flatten() # convert 2D array to 1D

# set up colour range and levels for the countour plot (min, max values of the temperature anomaly trends to display)
cmin, cmax, incr = -2.0, 2.0, 0.25 
levels = np.arange(cmin, cmax + incr, incr) 
norm = BoundaryNorm(levels, ncolors=256) # distinct intervals 
cmap = 'RdYlBu_r'

for i, var in enumerate(variables):
    ax = axes[i]
    ax.set_global()
    ax.coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='lightgray', linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12}
    gl.ylabel_style = {'size': 12}

    # plot data
    # contour = ax.contour(warming_contr_ds['lon'], warming_contr_ds['lat'], warming_contr_ds[var], levels=levels, cmap=cmap, extend='both')
    pc = ax.pcolormesh(warming_contr_ds['lon'], warming_contr_ds['lat'], warming_contr_ds[var], cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax.set_title(long_names[i], fontsize=18)  
    # vmin=cmin, vmax=cmax, 

# remove the last subplot as it is extra
fig.delaxes(axes[-1])

plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.1, wspace=0.1, hspace=0.07)

# add colourbar
cbar_ax = fig.add_axes([0.25, 0.05, 0.5, 0.02]) # [left, bottom, width, height]
cbar = fig.colorbar(pc, cax=cbar_ax, orientation='horizontal')
cbar.ax.tick_params(labelsize=24)
cbar.set_label('Warming Contribution (K)', fontsize=18)
labels = np.arange(cmin, cmax + incr, incr * 2)
cbar.set_ticks(labels)
cbar.set_ticklabels([f"{label:.1f}" for label in labels])  

fig.suptitle(f'Global Map of Warming Contributions for 1979-2023 in {model},{ensemble}', fontsize=26, y=0.97)

plt.savefig(f'figures/{model}_{ensemble}_map_warming_contr.png', dpi=300)
