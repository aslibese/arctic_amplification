#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon May 20 2024

This script calculates the warming contributions (K) of temperature, water vapour, albedo and cloud feedbacks, ERF,
AHT and OHU in each grid cell following the method from Hanh et al. (2021) and saves them as a single NetCDF file.

@author: aslibese
"""

import xarray as xr
from FeedbackUtilities import weighted_avg
import os

# load the file containing the energetic contributions of feedbacks and change in AHT and OHU
ds_path = 'data/CMIP6/MIROC6_feedbacks.nc'
ds = xr.open_dataset(ds_path)
# extract model name from the feedbacks path
model = os.path.basename(ds_path).split('_')[0]

# load Planck_lambda dataset
ds_pl = xr.open_dataset('data/CMIP6/MIROC6_Planck_lambda.nc')
# calculate global-mean Planck feedback (W m-2 K-1)
global_avg_Planck = weighted_avg(ds_pl['Planck_lambda'])

# calculate warming contributions (K) in each grid cell
LR_warming_contr = ds['LR'] / abs(global_avg_Planck)
Planck_warming_contr = ds['Planck'] / abs(global_avg_Planck)
q_lw_warming_contr = ds['q_lw'] / abs(global_avg_Planck)
q_sw_warming_contr = ds['q_sw'] / abs(global_avg_Planck)
alb_warming_contr = ds['alb'] / abs(global_avg_Planck)
ERF_lw_warming_contr = ds['ERF_lw'] / abs(global_avg_Planck)
ERF_sw_warming_contr = ds['ERF_sw'] / abs(global_avg_Planck)
OHU_warming_contr = ds['OHU'] / abs(global_avg_Planck)
AHT_warming_contr = ds['AHT'] / abs(global_avg_Planck)
cld_lw_warming_contr = ds['cld_lw'] / abs(global_avg_Planck)
cld_sw_warming_contr = ds['cld_sw'] / abs(global_avg_Planck)

# create a dataset to hold the warming contributions of all feedbacks (lat, lon)
warming_contr_ds = xr.Dataset({
    'LR': LR_warming_contr,
    'Planck': Planck_warming_contr,
    'q_lw': q_lw_warming_contr,
    'q_sw': q_sw_warming_contr,
    'alb': alb_warming_contr,
    'cld_lw': cld_lw_warming_contr,
    'cld_sw': cld_sw_warming_contr,
    'ERF_lw': ERF_lw_warming_contr,
    'ERF_sw': ERF_sw_warming_contr,
    'OHU': OHU_warming_contr,
    'AHT': AHT_warming_contr
})

variable_attrs = {
    'LR': {'long_name': "Warming contributions from lapse rate feedback", 'standard_name': "lapse_rate_feedback_warm_contr", 'units': "K"},
    'Planck': {'long_name': "Warming contributions from Planck feedback", 'standard_name': "planck_feedback_warm_contr", 'units': "K"},
    'q_lw': {'long_name': "Warming contributions from LW component of water vapour feedback", 'standard_name': "water_vapour_feedback_lw_warm_contr", 'units': "K"},
    'q_sw': {'long_name': "Warming contributions from SW component of water vapour feedback", 'standard_name': "water_vapour_feedback_sw_warm_contr", 'units': "K"},
    'alb': {'long_name': "Warming contributions from surface albedo feedback", 'standard_name': "albedo_feedback_warm_contr", 'units': "K"},
    'cld_lw': {'long_name': "Warming contributions from LW component of cloud feedback", 'standard_name': "cloud_feedback_lw_warm_contr", 'units': "K"},
    'cld_sw': {'long_name': "Warming contributions from SW component of cloud feedback", 'standard_name': "cloud_feedback_sw_warm_contr", 'units': "K"},
    'ERF_lw': {'long_name': "Warming contributions from LW effective radiative forcing", 'standard_name': "effective_radiative_forcing_lw_warm_contr", 'units': "K"},
    'ERF_sw': {'long_name': "Warming contributions from SW effective radiative forcing", 'standard_name': "effective_radiative_forcing_sw_warm_contr", 'units': "K"},
    'OHU': {'long_name': "Warming contributions from oceanic heat uptake", 'standard_name': "oceanic_heat_uptake_warm_contr", 'units': "K"},
    'AHT': {'long_name': "Warming contributions from atmospheric heat transport", 'standard_name': "atmospheric_heat_transport_warm_contr", 'units': "K"}
}

for var_name, attrs in variable_attrs.items():
    warming_contr_ds[var_name].attrs.update(attrs)

# save the dataset to a NetCDF file
warming_contr_ds.to_netcdf(f'data/CMIP6/{model}_warming_contr.nc')

print(f"data/CMIP6/{model}_warming_contr.nc saved successfully.") 