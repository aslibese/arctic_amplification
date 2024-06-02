#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun 2 2024

This script serves to access CMIP6 data stored on Google Cloud based on the specified models and variables 
for the historical and SSP5-8.5 experiments, concatenate them along the time axis, slice the dataset for the 1979-2023 period,
perform vertical and horizontal regridding, and download the output as a NetCDF file. 

@author: aslibese
"""

import xarray as xr
import pandas as pd
import os
from CMIP6_Utilities import open_concat_datasets, cdo_vertical_interpolation, cdo_bilinear_regridding, cdo_conservative_regridding


model = 'CanESM5'
plev_vars = ['ta', 'hus']
other_vars = ['tas', 'ts', 'ps','rsdt', 'rsut', 'rlut', 'rsutcs', 'rlutcs',
              'rsds', 'rsus', 'rlus', 'rlds', 'hfls', 'hfss']
all_vars = plev_vars + other_vars
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

# get plev variables
# open and concatenate datasets for historical and ssp585 (strongest forcing) experiments 
historical_ds_plev = open_concat_datasets(model, 'historical', 'r1i1p1f1', plev_vars)
ssp585_ds_plev = open_concat_datasets(model, 'ssp585', 'r1i1p1f1', plev_vars)
# concatenate historical and ssp585 along the time dimension 
concat_ds_plev = xr.concat([historical_ds_plev, ssp585_ds_plev], dim='time')
# select the 1979-2023 time period for our analysis
sliced_ds_plev = concat_ds_plev.sel(time=slice('1979-01-01', '2023-12-31'))

# get other variables
# open and concatenate datasets for historical and ssp585 (strongest forcing) experiments 
historical_ds = open_concat_datasets(model, 'historical', 'r1i1p1f1', other_vars)
ssp585_ds = open_concat_datasets(model, 'ssp585', 'r1i1p1f1', other_vars)
# concatenate historical and ssp585 along the time dimension 
concat_ds = xr.concat([historical_ds, ssp585_ds], dim='time')
# select the 1979-2023 time period for our analysis
sliced_ds = concat_ds.sel(time=slice('1979-01-01', '2023-12-31'))

# extract datetime values from CFTimeIndex
time_values = sliced_ds_plev.indexes['time'].to_datetimeindex().to_pydatetime()
# create a pandas DateTimeIndex from the extracted datetime values
sliced_ds_plev['time'] = pd.DatetimeIndex(time_values)
sliced_ds['time'] = pd.DatetimeIndex(time_values)

# target plev matched ERA5 kernels
target_plev = [100000, 97500, 95000, 92500, 90000, 87500, 85000, 82500, 80000, 77500, 
                75000, 70000, 65000, 60000, 55000, 50000, 45000, 40000, 35000, 30000, 
                25000, 22500, 20000, 17500, 15000, 12500, 10000, 7000, 5000, 3000, 2000, 
                1000, 700, 500, 300, 200, 100]

# create a temporary dataset for plev with time coordinate
temp_ds = xr.Dataset(coords={'time': sliced_ds_plev['time']})

for var in plev_vars:
    # perform vertical interpolation 
    input_file = f'{var}_input.nc'
    output_file = f'{var}_output.nc'

    sliced_ds_plev[var].to_netcdf(input_file)
    cdo_vertical_interpolation(input_file, output_file, target_plev)

    interpolated_data = xr.open_dataset(output_file)[var]
    # add the interpolated variable to the new dataset
    temp_ds[var] = interpolated_data
    # assign the correct plev values
    temp_ds[var] = temp_ds[var].assign_coords(plev=target_plev)

    # remove the intermediate files
    os.remove(input_file)
    os.remove(output_file)


bilinear_var = ['tas', 'ts', 'ta', 'ps', 'hus']
conservative_var = ['rsdt', 'rsut', 'rlut', 'rsutcs', 'rlutcs', 'rsds', 'rsus', 'rlus', 'rlds', 'hfls', 'hfss']

# create a new dataset with time coordinate
new_ds = xr.Dataset(coords={'time': sliced_ds['time']})

for var in all_vars:
    # for ta and hus, perform horizontal regridding in each plev separately and then combine
    if var in ['ta', 'hus']:  
        regridded_levels = []
        for plev in target_plev:
            input_file = f'{var}_input_{plev}.nc'
            output_file = f'{var}_output_{plev}.nc'

            # save the variable at the specific pressure level to a NetCDF file
            temp_ds[var].sel(plev=plev).to_netcdf(input_file)

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
            new_ds = new_ds.assign_coords(plev=target_plev)
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

new_ds.to_netcdf(f'data/CMIP6/{model}_regridded.nc')

print(new_ds.coords)