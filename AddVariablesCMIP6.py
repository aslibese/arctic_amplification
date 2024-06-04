#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun 2 2024

This script serves to calculate the missing variables needed for feedback calculations.

@author: aslibese
"""

import xarray as xr
import argparse
import os
from CMIP6_Utilities import make_tropo

parser = argparse.ArgumentParser(description='A program that requires the file name of the CMIP6 regridded dataset.')
parser.add_argument('filename', type=str, help='Name of the data file')

args = parser.parse_args()

print("Processing data file:", args.filename)

ds = xr.open_dataset(args.filename)

# calculate missing variables
# calculate net LW flux at TOA ('rlnt'), ensuring positive values represent downwards flux
ds['rlnt'] = -ds['rlut']
ds['rlnt'].attrs.update({
    'long_name': "TOA Net Longwave Radiative Flux",
    'standard_name': "toa_net_longwave_flux",
    'original_name': "FLNT"
})

# calculate clear-sky net LW flux at TOA ('rlntcs'), ensuring positive = downwards
ds['rlntcs'] = -ds['rlutcs']
ds['rlntcs'].attrs.update({
    'long_name': "TOA Clear-Sky Net Longwave Radiative Flux",
    'standard_name': "toa_net_longwave_flux_assuming_clear_sky",
    'original_name': "FLNTC",
    'comment': "Net clear-sky longwave radiation at top of atmosphere"
})

# calculate net SW flux at TOA ('rsnt'), ensuring positive = downwards 
ds['rsnt'] = ds['rsdt'] - ds['rsut']
ds['rsnt'].attrs.update({
    'long_name': "TOA Net Shortwave Radiative Flux",
    'standard_name': "toa_net_shortwave_flux",
    'original_name': "FSNT",
    'units': "W m-2"
})

# calculate clear-sky net SW flux at TOA ('rsntcs'), ensuring positive = downwards 
ds['rsntcs'] = ds['rsdt'] - ds['rsutcs'] 
ds['rsntcs'].attrs.update({
    'long_name': "TOA Clear-Sky Net Shortwave Radiative Flux",
    'standard_name': "toa_net_shortwave_flux_assuming_clear_sky",
    'original_name': "FSNTC",
    'units': "W m-2"
})

# calculate net surface LW flux, ensuring positive = downwards
ds['rlns'] = ds['rlds'] - ds['rlus']
ds['rlns'].attrs.update({
    'long_name': "Surface Net Longwave Radiative Flux",
    'standard_name': "surface_net_longwave_flux",
    'original_name': "FLNS",
    'units': "W m-2"
})

# calculate net surface SW flux, ensuring positive = downwards
ds['rsns'] = ds['rsds'] - ds['rsus']
ds['rsns'].attrs.update({
    'long_name': "Surface Net Shortwave Radiative Flux",
    'standard_name': "surface_net_shortwave_flux",
    'original_name': "FSNS",
    'units': "W m-2"
})

# calculate tropopause pressure
p_trop = make_tropo(ds['ps'])

# add p_trop as a new variable in the dataset 
ds['trop'] = xr.DataArray(
    p_trop,
    dims=('time', 'lat', 'lon'),
    coords={
        'time': ds['time'],
        'lat': ds['lat'],
        'lon': ds['lon']
    },
    attrs={'units': "Pa", 'long_name': "Tropopause Pressure"}
)

# extract base filename
base_filename = os.path.basename(args.filename)
# extract model and ensemble names from the path
model = base_filename.split('_')[0]
ensemble = base_filename.split('_')[1]

ds.to_netcdf(f'data/CMIP6/{model}_{ensemble}_regridded_wvars.nc')
print(f"data/CMIP6/{model}_{ensemble}_regridded_wvars.nc saved successfully.")