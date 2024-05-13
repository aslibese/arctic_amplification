#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat May 11 2024

This driver script imports the plot_background function from PolarMapMultiple.py
and contains the functions serving to plot temperature trend and amplification as polar maps using a single dataset.

Code adapted from mikarant's GitHub repository: https://github.com/mikarant/arctic-amplification/tree/main
Original code: https://github.com/mikarant/arctic-amplification/blob/main/fig1_manuscript.py

@author: aslibese
"""

import numpy as np   
import matplotlib.pyplot as plt  
from matplotlib import cm
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs # crs: coordinate reference system
from cartopy.util import add_cyclic_point

from PolarMapMultiple import plot_background

# function to plot the Arctic temperature trend as a polar map
def plot_average_ArcticTrend(ds, lower_boundary, add_colorbar=True, ax=None):
	# set up colour range and levels for the countour plot 
	cmin = -1.5 # min countour level (min value of the temperature anomaly trends to display)
	cmax = 1.5 # max countour level
	incr = 0.25 # increment between countour levels
	c_levels = np.arange(cmin, cmax + incr, incr) 
	trend_cmap = 'RdYlBu_r'

	# call the plot_background() function to create a background plot 
	plot_background(ax, lower_boundary)

	# use add_cyclic_point to avoid a gap at the longtitude seam
	# ds['slope']*10 to convert per-year trend to per-decade 
	data_cyclic_point, cyclic_lon = add_cyclic_point(ds['slope']*10, coord=ds['lon'])

	# create filled contour plot 
	filled_contourf = ax.contourf(cyclic_lon, ds['lat'], data_cyclic_point, levels=c_levels, 
		zorder=2, extend='both', # add arrow extentions on both ends of the colourbar to handle values outside of the specificied range
		cmap=trend_cmap, transform=ccrs.PlateCarree())

	cbar = plt.colorbar(filled_contourf, orientation='horizontal', pad=0.05, fraction=0.05)
	cbar.ax.tick_params(labelsize=14)
	cbar.set_label(label='Temperature trend [°C decade⁻¹]', fontsize=16)
	labels = np.arange(cmin, cmax + incr, incr * 3)
	cbar.set_ticks(labels)
	
   
# function to plot the Arctic Amplification (AA) as a polar map
def plot_average_AA(ds, lower_boundary, ax=None):
	# generate base colourmaps
	Reds = cm.get_cmap('Reds', 12) # make 12 discrete shades of red 
	Blues = cm.get_cmap('Blues', 10) # make 10 discrete shades of blue
	
	# modify colourmap colours
	newcolours = Reds(np.linspace(0,1,12)) # array created from 'Reds'
	newcolours[0,:] = Blues(0.7) # replace the first index with dark blue shade: AA = 0-0.5
	newcolours[1,:] = Blues(0.2) # replace the second index with lighter blue shade: AA = 0.5-1

	# loop to create a stepped effect in the colour transition for AA values larger than 1
	for i in range(2,12,2):
		newcolours[i] = newcolours[i+1, :]
	
	newcolours = np.vstack((newcolours, newcolours[-1], newcolours[-1])) # duplicate the last colour twice (highest AA values)
	newcolours[-2:, :3] = newcolours[-2:, :3] * 0.4 # set the last two colours as distinct colours by reducing the brigthness

	new_colourmap = ListedColormap(newcolours) # create a new colourmap from 'newcolours'

	# set out-of-range values below and above the colourmap's range using shades from the 'Blues' and 'Greys' colourmaps
	new_colourmap.set_under(cm.get_cmap('Blues')(0.99), 1.0) 
	new_colourmap.set_over(cm.get_cmap('Greys')(0.75), 1.0)

	# display AA values from 0 to 7 with 0.5 increments
	cmin = 0
	cmax = 7
	incr = 0.25
	c_levels = np.arange(cmin, cmax + incr, incr)

	plot_background(ax, lower_boundary)

	# use add_cyclic_point to avoid a gap at the longtitude seam
	data_cyclic_point, cyclic_lon = add_cyclic_point(ds['amplification'], coord=ds['lon'])

	# create filled contourmap
	filled_contourf = ax.contourf(cyclic_lon, ds['lat'], data_cyclic_point, levels=c_levels, zorder=2, 
		extend = 'both', cmap=new_colourmap, transform=ccrs.PlateCarree())

	cbar = plt.colorbar(filled_contourf, orientation='horizontal', pad=0.05, fraction=0.05)
	cbar.ax.tick_params(labelsize=14)
	cbar.set_label(label='Local amplification', fontsize=16)
	labels = np.arange(cmin, cmax + incr, 1)
	cbar.set_ticks(labels)