#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon May 20 2024

This script loads the CMIP6 perturbed dataset and calculates the energetic contributions (W m-2) of temperature, 
water vapour, albedo and cloud feedbacks and the change in AHT, OHT and ERF to the surface temperature warming during 1979-2023. 
It saves the energetic contributions data (lat, lon) and Planck_lambda data (lat, lon) as NetCDF files.

@author: aslibese
"""

import xarray as xr
import numpy as np
import climkern as ck 
import os
import argparse # to handle command line arguments
from FeedbackUtilities import make_clim, regress_against_ts, calculate_delta, weighted_avg

parser = argparse.ArgumentParser(description='A program that requires the file name of the CMIP6 perturbed dataset.')
parser.add_argument('filename', type=str, help='Name of the data file')

args = parser.parse_args()

print("Processing data file:", args.filename)

# open the CMIP6 perturbed dataset
pert_ds = xr.open_dataset(args.filename)

# select the variables needed for feedback calculations
pert_ta = pert_ds['ta'] # air temperature at various pressure levels
pert_ts = pert_ds['ts'] # skin (surface) temperature
pert_ps = pert_ds['ps'] # surface air pressure
pert_q = pert_ds['hus'] # specific humidity at various pressure levels
pert_trop = pert_ds['trop'] # tropopause pressure
pert_rsus = pert_ds['rsus'] # upwelling SW radiation at the surface
pert_rsds = pert_ds['rsds'] # downwelling SW radiation at the surface 

pert_rlnt = pert_ds['rlnt'] # net LW flux at TOA
pert_rlntcs = pert_ds['rlntcs'] # clear-sky net LW flux at TOA
pert_rsnt = pert_ds['rsnt'] # net SW flux at TOA
pert_rsntcs = pert_ds['rsntcs'] # clear-sky net SW flux at TOA

pert_rlns = pert_ds['rlns'] # net LW flux at the surface
pert_rsns = pert_ds['rsns'] # net SW flux at the surface
pert_hfls = pert_ds['hfls'] # surface upward latent heat flux
pert_hfss = pert_ds['hfss'] # surface upward sensible heat flux


# create control dataset (monthly climatologies) of the perturbed variables 
ctrl_ta = make_clim(pert_ta)
ctrl_ts = make_clim(pert_ts)
ctrl_ps = make_clim(pert_ps)
ctrl_q = make_clim(pert_q)
ctrl_rsus = make_clim(pert_rsus)
ctrl_rsds = make_clim(pert_rsds)

ctrl_rlnt = make_clim(pert_rlnt)
ctrl_rlntcs = make_clim(pert_rlntcs)
ctrl_rsnt = make_clim(pert_rsnt)
ctrl_rsntcs = make_clim(pert_rsntcs)

ctrl_rlns = make_clim(pert_rlns)
ctrl_rsns = make_clim(pert_rsns)
ctrl_hfls = make_clim(pert_hfls)
ctrl_hfss = make_clim(pert_hfss)

years = range(1979, 2024)
original_times = pert_ds['time']

LR_list = []
Planck_list = []
q_lw_list = []
q_sw_list = []
alb_list = []
LR_cs_list = []
Planck_cs_list = []
q_lw_cs_list = []
q_sw_cs_list = []
alb_cs_list = []
anom_ts_list = []
anom_rlntcs_list = []
anom_rsntcs_list = []
dCRE_LW_list = []
dCRE_SW_list = []

for year in years:
    pert_ta_year = pert_ta.where(pert_ta['time'].dt.year == year, drop=True)
    pert_ts_year = pert_ts.where(pert_ts['time'].dt.year == year, drop=True)
    pert_ps_year = pert_ps.where(pert_ps['time'].dt.year == year, drop=True)
    pert_trop_year = pert_trop.where(pert_trop['time'].dt.year == year, drop=True)
    pert_q_year = pert_q.where(pert_q['time'].dt.year == year, drop=True)
    pert_ps_year = pert_ps.where(pert_ps['time'].dt.year == year, drop=True)

    pert_rsus_year = pert_rsus.where(pert_rsus['time'].dt.year == year, drop=True)
    pert_rsds_year = pert_rsds.where(pert_rsds['time'].dt.year == year, drop=True)

    pert_rlntcs_year = pert_rlntcs.where(pert_rlntcs['time'].dt.year == year, drop=True)
    pert_rsntcs_year = pert_rsntcs.where(pert_rsntcs['time'].dt.year == year, drop=True)
    pert_rlntcs_year['time'] = pert_rlntcs_year['time'].dt.month # convert time to months
    pert_rsntcs_year['time'] = pert_rsntcs_year['time'].dt.month

    pert_rlnt_year = pert_rlnt.where(pert_rlnt['time'].dt.year == year, drop=True)
    pert_rsnt_year = pert_rsnt.where(pert_rsnt['time'].dt.year == year, drop=True)
    pert_rlnt_year['time'] = pert_rlnt_year['time'].dt.month
    pert_rsnt_year['time'] = pert_rsnt_year['time'].dt.month

    pert_rlns_year = pert_rlns.where(pert_rlns['time'].dt.year == year, drop=True)
    pert_rsns_year = pert_rsns.where(pert_rsns['time'].dt.year == year, drop=True)
    pert_hfls_year = pert_hfls.where(pert_hfls['time'].dt.year == year, drop=True)
    pert_hfss_year = pert_hfss.where(pert_hfss['time'].dt.year == year, drop=True)

    # calculate temperature feedbacks (LR and Planck)
    LR_year, Planck_year = ck.calc_T_feedbacks(ctrl_ta, ctrl_ts, ctrl_ps,
                                               pert_ta_year, pert_ts_year, pert_ps_year, pert_trop_year,
                                               kern="ERA5")

    # calculate water vapor feedbacks
    q_lw_year, q_sw_year = ck.calc_q_feedbacks(ctrl_q, ctrl_ta, ctrl_ps,
                                               pert_q_year, pert_ps_year, pert_trop_year,
                                               kern="ERA5", method="zelinka")

    # calculate albedo feedback
    alb_year = ck.calc_alb_feedback(ctrl_rsus, ctrl_rsds,
                                    pert_rsus_year, pert_rsds_year,
                                    kern="ERA5")

    # calculate clear-sky versions of the temperature, water vapor, and surface albedo feedbacks
    LR_cs_year, Planck_cs_year = ck.calc_T_feedbacks(ctrl_ta, ctrl_ts, ctrl_ps, 
                                                     pert_ta_year, pert_ts_year, pert_ps_year, pert_trop_year,
                                                     kern="ERA5", sky="clear-sky")

    q_lw_cs_year, q_sw_cs_year = ck.calc_q_feedbacks(ctrl_q, ctrl_ta, ctrl_ps,
                                                     pert_q_year, pert_ps_year, pert_trop_year,
                                                     kern="ERA5", method="zelinka", sky="clear-sky")

    alb_cs_year = ck.calc_alb_feedback(ctrl_rsus, ctrl_rsds,
                                       pert_rsus_year, pert_rsds_year,
                                       kern="ERA5", sky="clear-sky")
    
    # calculate clear-sky LW and SW radiation anomalies at TOA (for ERF calculations)
    anom_rlntcs_year = pert_rlntcs_year.groupby('time') - ctrl_rlntcs
    anom_rsntcs_year = pert_rsntcs_year.groupby('time') - ctrl_rsntcs
    anom_rlntcs_year['time'] = original_times.where(original_times.dt.year == year, drop=True)
    anom_rsntcs_year['time'] = original_times.where(original_times.dt.year == year, drop=True)

    # align the DataArrays
    ctrl_rlnt_aligned, pert_rlnt_year_aligned = xr.align(ctrl_rlnt, pert_rlnt_year, join='exact')
    ctrl_rsnt_aligned, pert_rsnt_year_aligned = xr.align(ctrl_rsnt, pert_rsnt_year, join='exact')
    ctrl_rlntcs_aligned, pert_rlntcs_year_aligned = xr.align(ctrl_rlntcs, pert_rlntcs_year, join='exact')
    ctrl_rsntcs_aligned, pert_rsntcs_year_aligned = xr.align(ctrl_rsntcs, pert_rsntcs_year, join='exact')

    # calculate the change in LW and SW cloud radiative effects (for cloud feedback calculations)
    dCRE_LW_year = ck.calc_dCRE_LW(ctrl_rlnt_aligned, pert_rlnt_year_aligned, ctrl_rlntcs_aligned, pert_rlntcs_year_aligned)
    dCRE_SW_year = ck.calc_dCRE_SW(ctrl_rsnt_aligned, pert_rsnt_year_aligned, ctrl_rsntcs_aligned, pert_rsntcs_year_aligned)
    dCRE_LW_year['time'] = original_times.where(original_times.dt.year == year, drop=True)
    dCRE_SW_year['time'] = original_times.where(original_times.dt.year == year, drop=True)

    # store the radiative contributions of feedbacks and anomalies (W m-2)
    LR_list.append(LR_year)
    Planck_list.append(Planck_year)
    q_lw_list.append(q_lw_year)
    q_sw_list.append(q_sw_year)
    alb_list.append(alb_year)
    LR_cs_list.append(LR_cs_year)
    Planck_cs_list.append(Planck_cs_year)
    q_lw_cs_list.append(q_lw_cs_year)
    q_sw_cs_list.append(q_sw_cs_year)
    alb_cs_list.append(alb_cs_year)
    anom_rlntcs_list.append(anom_rlntcs_year)
    anom_rsntcs_list.append(anom_rsntcs_year)
    dCRE_LW_list.append(dCRE_LW_year)
    dCRE_SW_list.append(dCRE_SW_year)


# calculate ts anomalies 
for year in years:
    pert_ts_year = pert_ts.where(pert_ts['time'].dt.year == year, drop=True)
    pert_ts_year['time'] = pert_ts_year['time'].dt.month

    # calculate 'ts' anomalies and restore the original time coordinates
    anom_ts_year = pert_ts_year.groupby('time') - ctrl_ts
    anom_ts_year['time'] = original_times.where(original_times.dt.year == year, drop=True)

    # store the anomalies 
    anom_ts_list.append(anom_ts_year) # (540, lat, lon)


# combine the feedbacks into a single xarray DataArray along a new 'year' dimension
LR_rad = xr.concat(LR_list, dim='time') # radiative contributions (W m-2)
Planck_rad = xr.concat(Planck_list, dim='time')
q_lw_rad = xr.concat(q_lw_list, dim='time')
q_sw_rad = xr.concat(q_sw_list, dim='time')
alb_rad = xr.concat(alb_list, dim='time')
LR_cs_rad = xr.concat(LR_cs_list, dim='time')
Planck_cs_rad = xr.concat(Planck_cs_list, dim='time')
q_lw_cs_rad = xr.concat(q_lw_cs_list, dim='time')
q_sw_cs_rad = xr.concat(q_sw_cs_list, dim='time')
alb_cs_rad = xr.concat(alb_cs_list, dim='time')
# combine anomalies 
anom_rlntcs = xr.concat(anom_rlntcs_list, dim='time')
anom_rsntcs = xr.concat(anom_rsntcs_list, dim='time')
dCRE_LW = xr.concat(dCRE_LW_list, dim='time')
dCRE_SW = xr.concat(dCRE_SW_list, dim='time')
anom_ts = xr.concat(anom_ts_list, dim='time')


# take the annual mean of the feedbacks and anomalies (45, 64, 128) 
# to reduce the noise associated with seasonal variations
LR_rad_ann = LR_rad.groupby('time.year').mean('time') # (year, lat, lon)
Planck_rad_ann = Planck_rad.groupby('time.year').mean('time')
q_lw_rad_ann = q_lw_rad.groupby('time.year').mean('time')
q_sw_rad_ann = q_sw_rad.groupby('time.year').mean('time')
alb_rad_ann = alb_rad.groupby('time.year').mean('time')
LR_cs_rad_ann = LR_cs_rad.groupby('time.year').mean('time')
Planck_cs_rad_ann = Planck_cs_rad.groupby('time.year').mean('time')
q_lw_cs_rad_ann = q_lw_cs_rad.groupby('time.year').mean('time')
q_sw_cs_rad_ann = q_sw_cs_rad.groupby('time.year').mean('time')
alb_cs_rad_ann = alb_cs_rad.groupby('time.year').mean('time')
anom_ts_ann = anom_ts.groupby('time.year').mean('time') 
anom_rlntcs_ann = anom_rlntcs.groupby('time.year').mean('time')
anom_rsntcs_ann = anom_rsntcs.groupby('time.year').mean('time')

# calculate the feedback parameter (λ) in each grid cell (W m-2 K-1) (lat, lon)
LR_lambda = regress_against_ts(anom_ts_ann, LR_rad_ann)
Planck_lambda = regress_against_ts(anom_ts_ann, Planck_rad_ann)
q_lw_lambda = regress_against_ts(anom_ts_ann, q_lw_rad_ann)
q_sw_lambda = regress_against_ts(anom_ts_ann, q_sw_rad_ann)
alb_lambda = regress_against_ts(anom_ts_ann, alb_rad_ann)
LR_cs_lambda = regress_against_ts(anom_ts_ann, LR_cs_rad_ann)
Planck_cs_lambda = regress_against_ts(anom_ts_ann, Planck_cs_rad_ann)
q_lw_cs_lambda = regress_against_ts(anom_ts_ann, q_lw_cs_rad_ann)
q_sw_cs_lambda = regress_against_ts(anom_ts_ann, q_sw_cs_rad_ann)
alb_cs_lambda = regress_against_ts(anom_ts_ann, alb_cs_rad_ann)

# convert Planck_lambda to a DataArray 
PL_da = xr.DataArray(
    Planck_lambda, 
    coords={'lat': Planck_rad_ann.lat, 'lon': Planck_rad_ann.lon}, 
    dims=['lat', 'lon'], 
    name='Planck_lambda',
    )
PL_da.attrs['long_name'] = 'Planck feedback parameter'
PL_da.attrs['units'] = 'W m-2 K-1'

# λp′ (lambda prime) is the difference between the regional Planck feedback (lambda) and its global-mean value (lambda bar)
global_avg_Planck = weighted_avg(PL_da) # scalar value
Planck_lambda_prime = PL_da - global_avg_Planck # (lat, lon)

# Δ is the trend in each variable during 1979-2023, multiplied by the number of years (45)
# calculate calculate Δ for ts anomalies in each grid cell 
delta_anom_ts = xr.apply_ufunc(
    calculate_delta, 
    anom_ts_ann, 
    input_core_dims=[['year']], 
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float]
    )
delta_anom_ts = xr.DataArray(delta_anom_ts, coords=[anom_ts_ann.lat, anom_ts_ann.lon], dims=['lat', 'lon'])

# calculate energetic contributions of feedbacks (W m-2) in each grid cell (lat, lon)
LR = LR_lambda * delta_anom_ts 
Planck = Planck_lambda_prime * delta_anom_ts
q_lw = q_lw_lambda * delta_anom_ts
q_sw = q_sw_lambda * delta_anom_ts
alb = alb_lambda * delta_anom_ts 
LR_cs = LR_cs_lambda * delta_anom_ts 
Planck_cs = Planck_cs_lambda * delta_anom_ts
q_lw_cs = q_lw_cs_lambda * delta_anom_ts
q_sw_cs = q_sw_cs_lambda * delta_anom_ts
alb_cs = alb_cs_lambda * delta_anom_ts 

# calculate Ocean Heat Uptake (OHU) in W m-2
# OHT at each lat and lon is calculated as the sum of 'rlns', 'rsns', 'hfls', and 'hfss'
# hfls and hfss are positive when upward, rlns and rsns are positive when downward
# to align the directions, subtract hfls and hfss 
OHU = pert_rlns + pert_rsns - pert_hfls - pert_hfss # (540, lat, lon)

# calculate Atmospheric Heat Transport (AHT) in W m-2
# AHT at each lat and lon is calculated as the difference between OHT and the sum of 'rlnt' and 'rsnt' 
# (+) values: poleward heat transport; (-) values: equatorward heat transport; zero at the poles and equator
AHT = OHU - (pert_rlnt + pert_rsnt)

# calculate Δ for AHT (W m-2)
AHT_ann = AHT.groupby('time.year').mean('time') # annual values (45, lat, lon)
delta_AHT = xr.apply_ufunc(
    calculate_delta,
    AHT_ann,
    input_core_dims =[['year']],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float]
    )
delta_AHT = xr.DataArray(delta_AHT, coords=[AHT_ann.lat, AHT_ann.lon], dims=['lat', 'lon']) # (lat, lon)

# calculate Δ for OHU (W m-2)
OHU_ann = OHU.groupby('time.year').mean('time') # annual values
delta_OHU = xr.apply_ufunc(
    calculate_delta,
    OHU_ann,
    input_core_dims =[['year']],
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float]
    )
delta_OHU = xr.DataArray(delta_OHU, coords=[OHU_ann.lat, OHU_ann.lon], dims=['lat', 'lon'])

# calculate Δ for anom_rlntcs (W m-2)
delta_anom_rlntcs = xr.apply_ufunc(
    calculate_delta, 
    anom_rlntcs_ann, 
    input_core_dims=[['year']], 
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float]
    )
delta_anom_rlntcs = xr.DataArray(delta_anom_rlntcs, coords=[anom_rlntcs_ann.lat, anom_rlntcs_ann.lon], dims=['lat', 'lon'])

# calculate Δ for anom_rsntcs (W m-2)
delta_anom_rsntcs = xr.apply_ufunc(
    calculate_delta, 
    anom_rsntcs_ann, 
    input_core_dims=[['year']], 
    vectorize=True,
    dask='parallelized',
    output_dtypes=[float]
    )
delta_anom_rsntcs = xr.DataArray(delta_anom_rsntcs, coords=[anom_rsntcs_ann.lat, anom_rsntcs_ann.lon], dims=['lat', 'lon'])

# calculate clear-sky LW and SW ERF as described in Hanh et al. (2021):
# step 1: calculate the sum of clear-sky versions of the temperature, water vapor, and surface albedo feedbacks
sum_lwcs_feedbacks = LR_cs + Planck_cs + q_lw_cs # (lat, lon)
sum_swcs_feedbacks = alb_cs + q_sw_cs 

# step 2: calculate clear-sky LW and SW ERF as the difference between TOA radiation anomalies
# and the sum of kernel-derived clear-sky feedback energetic contribution
ERF_lwcs = delta_anom_rlntcs - sum_lwcs_feedbacks # (lat, lon)
ERF_swcs = delta_anom_rsntcs - sum_swcs_feedbacks

# step 3: calculate all-sky LW and SW ERF by dividing clear-sky ERF by 1.16 (Soden et al. (2008))
ERF_lwas = ERF_lwcs / 1.16
ERF_swas = ERF_swcs / 1.16  

# calculate cloud SW and LW feedbacks
cld_lw_list = []
cld_sw_list = []
for year in years:
    LR_year = LR_rad.where(LR_rad['time'].dt.year == year, drop=True)
    Planck_year = Planck_rad.where(Planck_rad['time'].dt.year == year, drop=True)
    LR_cs_year = LR_cs_rad.where(LR_cs_rad['time'].dt.year == year, drop=True)
    Planck_cs_year = Planck_cs_rad.where(Planck_cs_rad['time'].dt.year == year, drop=True)
    q_lw_year = q_lw_rad.where(q_lw_rad['time'].dt.year == year, drop=True)
    q_lw_cs_year = q_lw_cs_rad.where(q_lw_cs_rad['time'].dt.year == year, drop=True)
    dCRE_LW_year = dCRE_LW.where(dCRE_LW['time'].dt.year == year, drop=True)

    alb_year = alb_rad.where(alb_rad['time'].dt.year == year, drop=True)
    alb_cs_year = alb_cs_rad.where(alb_cs_rad['time'].dt.year == year, drop=True)
    q_sw_year = q_sw_rad.where(q_sw_rad['time'].dt.year == year, drop=True)
    q_sw_cs_year = q_sw_cs_rad.where(q_sw_cs_rad['time'].dt.year == year, drop=True)
    dCRE_SW_year = dCRE_SW.where(dCRE_SW['time'].dt.year == year, drop=True)

    anom_rlntcs_year = anom_rlntcs.where(anom_rlntcs['time'].dt.year == year, drop=True)
    anom_rsntcs_year = anom_rsntcs.where(anom_rlntcs['time'].dt.year == year, drop=True)

    # step 1: calculate the sum of clear-sky versions of the temperature, water vapor, and surface albedo feedbacks
    sum_lwcs_feedbacks_year = LR_cs_year + Planck_cs_year + q_lw_cs_year 
    sum_swcs_feedbacks_year = alb_cs_year + q_sw_cs_year 

    # step 2: calculate clear-sky LW and SW ERF as the difference between TOA radiation anomalies
    # and the sum of kernel-derived clear-sky feedback energetic contribution
    ERF_lwcs_year = anom_rlntcs_year - sum_lwcs_feedbacks_year 
    ERF_swcs_year = anom_rsntcs_year - sum_swcs_feedbacks_year

    # step 3: calculate all-sky LW and SW ERF by dividing clear-sky ERF by 1.16 (Soden et al. (2008))
    ERF_lwas_year = ERF_lwcs_year / 1.16
    ERF_swas_year = ERF_swcs_year / 1.16  

    # step 4: calculate cloud feedbacks
    cld_lw_year = ck.calc_cloud_LW(LR_year + Planck_year, LR_cs_year + Planck_cs_year, q_lw_year, q_lw_cs_year, 
                                   dCRE_LW_year, ERF_lwas_year, ERF_lwcs_year)

    cld_sw_year = ck.calc_cloud_SW(alb_year, alb_cs_year, q_sw_year, q_sw_cs_year, 
                                   dCRE_SW_year, ERF_swas_year, ERF_swcs_year)

    cld_lw_list.append(cld_lw_year)
    cld_sw_list.append(cld_sw_year)

cld_lw_rad = xr.concat(cld_lw_list, dim='time') # radiative contributions (W m-2), (540, lat, lon)
cld_sw_rad = xr.concat(cld_sw_list, dim='time') 

cld_lw_rad_ann = cld_lw_rad.groupby('time.year').mean('time') # (45, lat, lon)
cld_sw_rad_ann = cld_sw_rad.groupby('time.year').mean('time')

cld_lw_lambda = regress_against_ts(anom_ts_ann, cld_lw_rad_ann) # (lat, lon)
cld_sw_lambda = regress_against_ts(anom_ts_ann, cld_sw_rad_ann)

cld_lw = cld_lw_lambda * delta_anom_ts
cld_sw = cld_sw_lambda * delta_anom_ts


# create a dataset to hold the energetic contributions of feedbacks (W m-2) and change in AHT and OHU
feedbacks_ds = xr.Dataset({
    'LR': LR,
    'Planck': Planck,
    'q_lw': q_lw,
    'q_sw': q_sw,
    'alb': alb,
    'OHU' : delta_OHU,
    'AHT' : delta_AHT,
    'cld_lw': cld_lw,
    'cld_sw': cld_sw,
    'ERF_lw': ERF_lwas,
    'ERF_sw': ERF_swas,
})

variable_attrs = {
    'LR': {'long_name': "Energetic contributions of lapse rate feedback", 'standard_name': "lapse_rate_feedback", 'units': "W m-2"},
    'Planck': {'long_name': "Energetic contributions of Planck feedback", 'standard_name': "planck_feedback", 'units': "W m-2"},
    'q_lw': {'long_name': "Energetic contributions of longwave water vapour feedback", 'standard_name': "water_vapour_feedback_lw", 'units': "W m-2"},
    'q_sw': {'long_name': "Energetic contributions of shortwave water vapour feedbacky", 'standard_name': "water_vapour_feedback_sw", 'units': "W m-2"},
    'alb': {'long_name': "Energetic contributions of surface albedo feedback", 'standard_name': "albedo_feedback", 'units': "W m-2"},
    'OHU' : {'long_name': "Change in Oceanic Heat Uptake", 'standard_name': "oceanic_heat_uptake", 'units': "W m-2"},
    'AHT' : {'long_name': "Change in Atmospheric Heat Transport", 'standard_name': "atmospheric_heat_transport", 'units': "W m-2"},
    'cld_lw': {'long_name': "Energetic contributions of LW cloud feedback", 'standard_name': "cloud_feedback_lw", 'units': "W m-2"},
    'cld_sw': {'long_name': "Energetic contributions of SW cloud feedback", 'standard_name': "cloud_feedback_sw", 'units': "W m-2"},
    'ERF_lw': {'long_name': "LW effective radiative forcing", 'standard_name': "effective_radiative_forcing_lw", 'units': "W m-2"},
    'ERF_sw': {'long_name': "SW effective radiative forcing", 'standard_name': "effective_radiative_forcing_sw", 'units': "W m-2"},
}
for var_name, attrs in variable_attrs.items():
    feedbacks_ds[var_name].attrs.update(attrs)

# extract base filename
base_filename = os.path.basename(args.filename)

# extract model and ensemble names from the path
model = base_filename.split('_')[0]
ensemble = base_filename.split('_')[1]

# save the datasets as NetCDF files
feedbacks_ds.to_netcdf(f'data/CMIP6/{model}_{ensemble}_feedbacks.nc') # (W m-2)
print(f"data/CMIP6/{model}_{ensemble}_feedbacks.nc saved successfully.")

PL_da.to_netcdf(f'data/CMIP6/{model}_{ensemble}_Planck_lambda.nc') # (W m-2 K-1)
print(f"data/CMIP6/{model}_{ensemble}_Planck_lambda.nc saved successfully.")



 