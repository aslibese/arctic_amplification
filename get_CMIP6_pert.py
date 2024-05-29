#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon May 20 2024

This script serves to access CMIP6 data stored on Google Cloud based on the specified models and variables 
for the historical and SSP5-8.5 experiments, concatenate them along the time axis, calculate the missing variables 
needed for feedback calculations, slice the dataset for the 1979-2023 period and download it as a NetCDF file. 

@author: aslibese
"""

import xarray as xr
import pandas as pd
import gcsfs
import numpy as np

# function to open dataset
def open_dataset(fs, url):
    store = fs.get_mapper(url) # create a mutable-mapping-style interface to the store
    ds = xr.open_zarr(store, consolidated=True) # open it using xarray and zarr
    return ds

# initialize Google Cloud Storage file system
fs = gcsfs.GCSFileSystem(token='anon')

# load the CSV containing the metadata
csv_url = 'https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv'
df = pd.read_csv(csv_url)

# function to get URL for a specific dataset
def get_url(model, experiment, variables):
    subset = df[(df['table_id'] == 'Amon') & (df['source_id'] == model) & 
                (df['experiment_id'] == experiment) & (df['variable_id'] == variables)] 
    if not subset.empty:
        return subset.iloc[0]['zstore']
    else:
        return None

# specify the model and variables
model = 'MIROC6'
variables = ['tas', 'ts', 'ta', 'ps', 'hus', 
             'rsdt', 'rsut', 'rlut', 'rsutcs', 'rlutcs',
             'rsds', 'rsus', 'rlus', 'rlds',
             'hfls', 'hfss']
# tas: near-surface (2m) air temperature
# ts: skin (surface) temperature
# ta: air temperature at various pressure levels
# ps: surface air pressure
# hus: specific humidity at various pressure levels
# rsdt: downwelling SW radiation at TOA 
# rsut: upwelling SW radiation at TOA
# rlut: upwelling LW radiation at TOA 
# rsutcs: clear-sky upwelling SW radiation at TOA
# rlutcs: clear-sky upwelling LW radiation at TOA
# rsds: downwelling SW radiation at the surface 
# rsus: upwelling SW radiation at the surface
# rlus: upwelling LW radiation at the surface
# rlds: downwelling LW radiation at the surface
# rsuscs: clear-sky upwelling SW radiation at the surface - NOT USED
# rldscs: clear-sky downwelling LW radiation at the surface - NOT USED
# hfls: surface upward latent heat flux
# hfss: surface upward sensible heat flux


# function to open and concatenate datasets for all variables 
# ADD FOR ALL MODELS 
def open_concat_datasets(model, experiment, variables):
    datasets = []
    for var in variables:
        url = get_url(model, experiment, var)
        if url:
            ds = open_dataset(fs, url)
            datasets.append(ds)
    combined_ds = xr.merge(datasets)
    return combined_ds

# open and concatenate datasets for historical and ssp585 (strongest forcing) experiments 
historical_ds = open_concat_datasets(model, 'historical', variables)
ssp585_ds = open_concat_datasets(model, 'ssp585', variables)

# concatenate historical and ssp585 along the time dimension 
concat_ds = xr.concat([historical_ds, ssp585_ds], dim='time')

# calculate missing variables
# calculate net LW flux at TOA ('rlnt'), ensuring positive values represent downwards flux
concat_ds['rlnt'] = -concat_ds['rlut']
concat_ds['rlnt'].attrs.update({
    'long_name': "TOA Net Longwave Radiative Flux",
    'standard_name': "toa_net_longwave_flux",
    'original_name': "FLNT"
})

# calculate clear-sky net LW flux at TOA ('rlntcs'), ensuring positive = downwards
concat_ds['rlntcs'] = -concat_ds['rlutcs']
concat_ds['rlntcs'].attrs.update({
    'long_name': "TOA Clear-Sky Net Longwave Radiative Flux",
    'standard_name': "toa_net_longwave_flux_assuming_clear_sky",
    'original_name': "FLNTC",
    'comment': "Net clear-sky longwave radiation at top of atmosphere"
})

# calculate net SW flux at TOA ('rsnt'), ensuring positive = downwards 
concat_ds['rsnt'] = concat_ds['rsdt'] - concat_ds['rsut']
concat_ds['rsnt'].attrs.update({
    'long_name': "TOA Net Shortwave Radiative Flux",
    'standard_name': "toa_net_shortwave_flux",
    'original_name': "FSNT",
    'units': "W m-2"
})

# calculate clear-sky net SW flux at TOA ('rsntcs'), ensuring positive = downwards 
concat_ds['rsntcs'] = concat_ds['rsdt'] - concat_ds['rsutcs'] 
concat_ds['rsntcs'].attrs.update({
    'long_name': "TOA Clear-Sky Net Shortwave Radiative Flux",
    'standard_name': "toa_net_shortwave_flux_assuming_clear_sky",
    'original_name': "FSNTC",
    'units': "W m-2"
})

# calculate net surface LW flux, ensuring positive = downwards
concat_ds['rlns'] = concat_ds['rlds'] - concat_ds['rlus']
concat_ds['rlns'].attrs.update({
    'long_name': "Surface Net Longwave Radiative Flux",
    'standard_name': "surface_net_longwave_flux",
    'original_name': "FLNS",
    'units': "W m-2"
})

# calculate net surface SW flux, ensuring positive = downwards
concat_ds['rsns'] = concat_ds['rsds'] - concat_ds['rsus']
concat_ds['rsns'].attrs.update({
    'long_name': "Surface Net Shortwave Radiative Flux",
    'standard_name': "surface_net_shortwave_flux",
    'original_name': "FSNS",
    'units': "W m-2"
})

# function to calculate the tropopause pressure at different latitudes using the method described in the ClimKern package:
# tropopause height of 10,000 Pa at the equator, descending with the cosine of latitude to 30,000 Pa at the poles
def make_tropo(da):
    p_trop = (3e4 - 2e4 * np.cos(np.deg2rad(da.lat))).broadcast_like(da)
    return p_trop

# calculate tropopause pressure
p_trop = make_tropo(concat_ds['ps'])

# add p_trop as a new variable in the dataset 
concat_ds['trop'] = xr.DataArray(
    p_trop,
    dims=('time', 'lat', 'lon'),
    coords={
        'time': concat_ds['time'],
        'lat': concat_ds['lat'],
        'lon': concat_ds['lon']
    },
    attrs={'units': "Pa", 'long_name': "Tropopause Pressure"}
)

# select the 1979-2023 time period for our analysis
sliced_ds = concat_ds.sel(time=slice('1979-01-01', '2023-12-31'))

# the sliced dataset will represent the perturbed dataset, save it as a NetCDF file
sliced_ds.to_netcdf(f'data/CMIP6/pert_{model}.nc')

print(f"{model} perturbed dataset has been successfully created and saved to data/CMIP6/.")
