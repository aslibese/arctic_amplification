#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat May 11 2024

This driver script imports functions from PolarMapMultiple.py, 
reads in the netCDF files of the four observational datasets with temperature trend and amplification.
It generates two plots and saves them in the /figures directory:
1) multi-panel plots of decadal temperature trend 
2) multi-panel plots of local amplification  

@author: aslibese
"""

from PolarMapMultiple import plot_multipanel_ArcticTrend, plot_multipanel_AA
  
import xarray as xr 
import matplotlib.pyplot as plt

# open the datasets
ERA5_trend_masked = xr.open_dataset('data/era5_anom_trendMasked.nc')
ERA5_aa = xr.open_dataset('data/era5_anom_aa.nc')

HadCRUT5_trend_masked = xr.open_dataset('data/HadCRUT5_anom_trendMasked.nc')
HadCRUT5_aa = xr.open_dataset('data/HadCRUT5_anom_aa.nc')

Gistemp_trend_masked = xr.open_dataset('data/gistemp_anom_trendMasked.nc')
Gistemp_aa = xr.open_dataset('data/gistemp_anom_aa.nc')

BerkeleyEarth_trend_masked = xr.open_dataset('data/BerkeleyEarth_anom_trendMasked.nc')
BerkeleyEarth_aa = xr.open_dataset('data/BerkeleyEarth_anom_aa.nc')

datasets_trend = [Gistemp_trend_masked, BerkeleyEarth_trend_masked, HadCRUT5_trend_masked, ERA5_trend_masked]
datasets_aa = [Gistemp_aa, BerkeleyEarth_aa, HadCRUT5_aa, ERA5_aa]

subplot_labels = ['a) Gistemp', 'b) Berkeley Earth', 'c) HadCRUT5', 'd) ERA5']
arctic_lower = 66.5

# plot multi-panel Arctic Temperature Trend 
plot_multipanel_ArcticTrend(datasets_trend, arctic_lower, subplot_labels)
plt.savefig('figures/obs_subplots_temp_trend.png', dpi=300, bbox_inches='tight')
plt.close()

# plot multi-panel Arctic Amplification
plot_multipanel_AA(datasets_aa, arctic_lower, subplot_labels)
plt.savefig('figures/obs_subplots_aa.png', dpi=300, bbox_inches='tight')
plt.close()