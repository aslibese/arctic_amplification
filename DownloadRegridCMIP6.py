#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun 2 2024

This script serves to access CMIP6 data stored on Google Cloud based on the specified models and variables 
for the historical and SSP5-8.5 experiments, concatenate them along the time axis, slice the dataset for the 1979-2023 period,
perform horizontal regridding to match ERA5 kernels, and download the output as a NetCDF file. 

@author: aslibese
"""

import xarray as xr
import pandas as pd
import os
import sys
import numpy as np
import cftime
from dask.diagnostics import ProgressBar
from CMIP6_Utilities import open_concat_datasets, cdo_bilinear_regridding, cdo_conservative_regridding, apply_nearest_neighbour_parallel


model = 'UKESM1-0-LL'
ensemble = 'r1i1p1f2'
variables = ['ta','hus','tas', 'ts', 'ps','rsdt', 'rsut', 'rlut', 'rsutcs', 'rlutcs',
            'rsds', 'rsus', 'rlus', 'rlds', 'hfls', 'hfss']
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

with ProgressBar():
    # open and concatenate datasets for historical and ssp585 (strongest forcing) experiments 
    historical_ds = open_concat_datasets(model, 'historical', ensemble, variables)
    ssp585_ds = open_concat_datasets(model, 'ssp585', ensemble, variables)
    # concatenate historical and ssp585 along the time dimension 
    concat_ds = xr.concat([historical_ds, ssp585_ds], dim='time')

time_type = type(concat_ds['time'].values[0]) 
start_date = '1979-01-01' 
end_date = '2023-12-31'

# handle different cftime data types e.g., CanESM5 is cftime.DatetimeNoLeap
if time_type in [cftime.DatetimeNoLeap, cftime.Datetime360Day, cftime.DatetimeAllLeap, cftime.DatetimeGregorian, cftime.DatetimeProlepticGregorian, cftime.DatetimeJulian]:
    # adjust end_date for cftime.Datetime360Day calendar
    if time_type == cftime.Datetime360Day:
        end_date = '2023-12-30'
    sliced_ds = concat_ds.sel(time=slice(start_date, end_date))
    # extract datetime values from CFTimeIndex
    time_values = sliced_ds.indexes['time'].to_datetimeindex().to_pydatetime() 
    # create a pandas DateTimeIndex from the extracted datetime values
    sliced_ds['time'] = pd.DatetimeIndex(time_values)

# handle numpy datetime64 e.g., MIROC6
elif time_type == np.datetime64:
    start_date = pd.to_datetime(start_date) # convert string dates to pandas Timestamp
    end_date = pd.to_datetime(end_date)
    sliced_ds = concat_ds.sel(time=slice(start_date, end_date))
    # extract datetime values from CFTimeIndex
    time_values = sliced_ds.indexes['time'].to_pydatetime()
    # create a pandas DateTimeIndex from the extracted datetime values
    sliced_ds['time'] = pd.DatetimeIndex(time_values)

# handle pandas Timestamp
elif time_type == pd.Timestamp:
    start_date = pd.to_datetime(start_date)
    end_date = pd.to_datetime(end_date)
    sliced_ds = concat_ds.sel(time=slice(start_date, end_date))
    # ensure time is in DateTimeIndex
    if not isinstance(sliced_ds.indexes['time'], pd.DatetimeIndex):
        sliced_ds['time'] = pd.DatetimeIndex(sliced_ds.indexes['time'])

else:
    print(f'Error: datetime values in {model} is {time_type}.')
    sys.exit()  # exit the code

# save the data before interpolation and regridding
sliced_ds.to_netcdf(f'data/CMIP6/{model}_{ensemble}_raw.nc')
print(f"data/CMIP6/{model}_{ensemble}_raw.nc saved successfully.\n")

# fill the missing ta and hus values using the nearest neighbour interpolation method
with ProgressBar():
    for var in ['ta', 'hus']:
        if np.any(np.isnan(sliced_ds[var].values)): 
            var_filled = apply_nearest_neighbour_parallel(sliced_ds, var)
            sliced_ds[var].values = var_filled
            print('Interpolation step for missing {var} values completed.')
            # check if there are any NaN values left after interpolation
            if np.any(np.isnan(sliced_ds[var].values)): 
                nan_indices = np.where(np.isnan(sliced_ds[var].values))
                nan_indices_combined = list(zip(*nan_indices))
                print(f'NaN {var} values found after interpolation: {nan_indices_combined}')


# bilinear interpolation for scalar variables, conservative interpolation for fluxes 
bilinear_var = ['tas', 'ts', 'ta', 'ps', 'hus']
conservative_var = ['rsdt', 'rsut', 'rlut', 'rsutcs', 'rlutcs', 'rsds', 'rsus', 'rlus', 'rlds', 'hfls', 'hfss']

# create a new dataset with time coordinate
new_ds = xr.Dataset(coords={'time': sliced_ds['time']})

# perform horizontal regridding
for var in variables:
    # for ta and hus, perform horizontal regridding in each plev separately and then combine
    if var in ['ta', 'hus']:  
        pressure_levels = sliced_ds[var]['plev'].values
        regridded_levels = []
        for plev in pressure_levels:
            input_file = f'{var}_input_{plev}.nc'
            output_file = f'{var}_output_{plev}.nc'

            # save the variable at the specific pressure level to a NetCDF file
            sliced_ds[var].sel(plev=plev).to_netcdf(input_file)

            cdo_bilinear_regridding(input_file, output_file)

            # load the interpolated data back into the dataset
            regridded_var = xr.open_dataset(output_file)[var]
            regridded_levels.append(regridded_var)

            # remove the intermediate files
            os.remove(input_file)
            os.remove(output_file)

        # combine all pressure levels into a single DataArray
        combined_regridded = xr.concat(regridded_levels, dim='plev')

        # add plev as a coordinate if not already present
        if 'plev' not in new_ds.coords:
            new_ds = new_ds.assign_coords(plev=combined_regridded['plev'])
            new_ds['plev'].attrs['standard_name'] = "air_pressure"
            new_ds['plev'].attrs['long_name'] = "Pressure Level"
            new_ds['plev'].attrs['units'] = "Pa"
        
        new_ds[var] = combined_regridded

    else:  # variables without the 'plev' dimension
        input_file = f'{var}_input.nc'
        output_file = f'{var}_output.nc'

        # save the variable to a NetCDF file
        sliced_ds[var].to_netcdf(input_file)

        if var in bilinear_var:  
            cdo_bilinear_regridding(input_file, output_file)
        elif var in conservative_var:
            cdo_conservative_regridding(input_file, output_file)

        # load the interpolated data back into the dataset
        regridded_var = xr.open_dataset(output_file)[var]

        if 'lat' not in new_ds.coords:
            regridded_coords = xr.open_dataset(output_file)
            new_ds = new_ds.assign_coords(lat=regridded_coords['lat'], lon=regridded_coords['lon'])

        # Add the regridded variable to the new dataset
        new_ds[var] = regridded_var

        # Remove the intermediate files
        os.remove(input_file)
        os.remove(output_file)


time = new_ds['time']
lat = new_ds['lat']
lon = new_ds['lon']
plev = new_ds['plev']

for var in new_ds.variables:
    if 'time' in new_ds[var].dims:
        new_ds[var] = new_ds[var].assign_coords(time=time)
    if 'lat' in new_ds[var].dims:
        new_ds[var] = new_ds[var].assign_coords(lat=lat)
    if 'lon' in new_ds[var].dims:
        new_ds[var] = new_ds[var].assign_coords(lon=lon)
    if 'plev' in new_ds[var].dims:
        new_ds[var] = new_ds[var].assign_coords(plev=plev)

new_ds.to_netcdf(f'data/CMIP6/{model}_{ensemble}_regridded.nc')
print(f"data/CMIP6/{model}_{ensemble}_regridded.nc saved successfully.")

print(new_ds.coords)