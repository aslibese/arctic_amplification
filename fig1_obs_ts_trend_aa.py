#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat May 11 2024

This driver script imports functions from PolarMapMultiple.py and ArcticWarmingUtilities.py,
reads in the four observational datasets as netCDF files.
It generates a single plots with three sub-plots and saves it in the /figures directory:
a) time series showing the Arctic and global temperature anomaly trend during 1950-2023
b) polar map of decadal temperature trend   
c) polar map of local amplification

@author: aslibese
"""

from PolarMapSingle import plot_average_ArcticTrend, plot_average_AA
from ArcticWarmingUtilities import selectPeriod, find_variable_by_dims, annualMean, selectRegion, performLinearRegression, weightedAverage

import numpy as np  
import xarray as xr 
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.lines as mlines
import cartopy.crs as ccrs # crs: coordinate reference system

# set up the figure using GridSpec for a grid layout
fig = plt.figure(figsize=(10, 10), constrained_layout=False)
gs = GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.3, wspace=0.3)

# a) TOP: TIME SERIES
# open the individual observational datasets
BerkeleyEarth = xr.open_dataset('data/BerkeleyEarth_anom.nc')
HadCRUT5 = xr.open_dataset('data/HadCRUT5_anom.nc')
Gistemp = xr.open_dataset('data/gistemp_anom.nc')
ERA5 = xr.open_dataset('data/new_era5_anom.nc')

obs_ds = [BerkeleyEarth, HadCRUT5, Gistemp, ERA5]

# define acceptable dimensions
acceptable_dims = [
	{'time', 'lat', 'lon'},
	{'time', 'latitude', 'longitude'}
]

# define the Arctic lower boundary
arctic_lower = 66.5

labels = ['Berkeley Earth', 'HadCRUT5', 'Gistemp', 'ERA5']
colours = ['b', 'r', 'y', 'g']

ax1 = fig.add_subplot(gs[0, :])
print("Plotting subplot (a)...\n")

for i in range(4): 
    print("Analyzing and plotting: " + labels[i] + "\n")

    # select the dataset for the period interest
    ds = selectPeriod(obs_ds[i], 1950, 2023)

    # use the function to find the variable for temperature anomaly
    var_name = find_variable_by_dims(obs_ds[i], acceptable_dims)

    # take the annual mean of temperature anomalies
    ds_with_annualMean = annualMean(ds, var_name)

    # compute the yearly Global weighted average
    global_weighted_ave = weightedAverage(ds_with_annualMean)

    # compute the Global warming trend for 1979-2023 and fitted y values
    global_selected_per = global_weighted_ave.sel(year = slice(str(1979), str(2023)))
    global_trend, global_intercept, _ = performLinearRegression(global_selected_per.values)
    y_fit_global = np.arange(len(global_selected_per)) * global_trend + global_intercept

    # select the Arctic region based on the lower boundary
    ds_arctic = selectRegion(ds_with_annualMean, arctic_lower) 

    # compute the yearly Arctic weighed average 
    arctic_weighted_ave = weightedAverage(ds_arctic)
    
    # compute the Arctic warming trend for 1979-2023 and fitted y values
    arctic_selected_per = arctic_weighted_ave.sel(year = slice(str(1979), str(2023)))
    arctic_trend, arctic_intercept, _ = performLinearRegression(arctic_selected_per.values)
    y_fit_arctic = np.arange(len(arctic_selected_per)) * arctic_trend + arctic_intercept
	
    # plot the Global temperature anomaly time series and the trend line for each dataset
    ax1.plot(global_weighted_ave.year, global_weighted_ave, color=colours[i], alpha=0.2)
    ax1.plot(global_selected_per.year, y_fit_global, color=colours[i], alpha=0.2)

    # plot the Arctic temperature anomaly time series and the trend line for each dataset
    ax1.plot(arctic_weighted_ave.year, arctic_weighted_ave, color=colours[i], alpha=1)
    ax1.plot(arctic_selected_per.year, y_fit_arctic, color=colours[i], alpha=1)

# add legend and labels
legend_lines = [mlines.Line2D([0], [0], color=color, linestyle='solid') for color in colours] #  dummy lines for legend purposes
ax1.legend(handles=legend_lines, labels=labels, loc='upper left',fontsize=16)
ax1.grid(True, linestyle='-') 
ax1.set_ylabel('Temperature anomaly [°C]', fontsize=16) 
ax1.set_xticks(np.arange(1950, 2021, 10))
ax1.tick_params(axis='both', labelsize=16)
ax1.set_xlim(1950, 2025)

# b) BOTTOM LEFT PLOT: ARCTIC TREND
print("Plotting subplot (b)...\n")
obs_ave_trend_masked = xr.open_dataset('data/obs_ave_trendMasked.nc')

ax2 = fig.add_subplot(gs[1, 0], projection=ccrs.NorthPolarStereo())
# plot Arctic Temperature Trend using the average of four observational datasets
plot_average_ArcticTrend(obs_ave_trend_masked, arctic_lower, ax=ax2)

# c) BOTTOM RIGHT: ARCTIC AMPLIFICATION
print("Plotting subplot (c)...\n")
obs_ave_aa = xr.open_dataset('data/obs_ave_aa.nc')

ax3 = fig.add_subplot(gs[1, 1], projection=ccrs.NorthPolarStereo())
# plot Arctic Amplification using the average of four observational datasets
plot_average_AA(obs_ave_aa, arctic_lower, ax=ax3)

# annotations
ax1.text(0, 1.1, 'a)', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax2.text(0, 1.1, 'b)', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')
ax3.text(0, 1.1, 'c)', transform=ax3.transAxes, fontsize=16, fontweight='bold', va='top', ha='right')

plt.savefig('figures/obs_combined_arcticTrend_aa.png', dpi=300)
