#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sun Jun 2 2024

This utilities script contains the functions used in the DownloadRegridCMIP6.py and AddVariablesCMIP6.py scripts.

@author: aslibese
"""

import pandas as pd
import xarray as xr
import numpy as np
from scipy.interpolate import griddata
import gcsfs
import os
import warnings
import concurrent.futures
import dask.array as da
warnings.filterwarnings("ignore", message="Converting a CFTimeIndex with dates from a non-standard calendar")

# initialize Google Cloud Storage file system
fs = gcsfs.GCSFileSystem(token='anon')

# load the CSV containing the metadata
csv_url = 'https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv'
df = pd.read_csv(csv_url)

# function to get URL for a specific dataset
def get_url(model, experiment, ensemble, variables):
    subset = df[(df['table_id'] == 'Amon') & (df['source_id'] == model) & 
                (df['experiment_id'] == experiment) & (df['member_id'] == ensemble) & (df['variable_id'] == variables)] 
    if not subset.empty:
        return subset.iloc[0]['zstore']
    else:
        return None

# function to open dataset
def open_dataset(fs, url):
    store = fs.get_mapper(url) # create a mutable-mapping-style interface to the store
    ds = xr.open_zarr(store, consolidated=True) # open it using xarray and zarr
    return ds

# function to open and concatenate datasets for all variables 
def open_concat_datasets(model, experiment, ensemble, variables):
    datasets = []
    for var in variables:
        url = get_url(model, experiment, ensemble, var)
        if url:
            ds = open_dataset(fs, url)
            datasets.append(ds)
    try:
        combined_ds = xr.merge(datasets)  # try standard merge first

    # some models have conflicts in 'bnds' across different experiments (historical and ssp585)
    except xr.MergeError as e:
        if 'conflicting values for variable' in str(e):
            print("Conflicting values detected, using compat='override'")
            combined_ds = xr.merge(datasets, compat='override')
        else:
            raise e  # re-raise the error if it's not the specific MergeError we expect

    return combined_ds

# function to assign the nearest value to ta and hus values at missing plev
def nearest_neighbour(var_data, plev):
    var_filled = griddata(plev[np.isfinite(var_data)], var_data[np.isfinite(var_data)], plev, method="nearest", fill_value="extrapolate")
    return var_filled

# function to apply nearest neighbour interpolation to each time, lat, lon combination
def apply_nearest_neighbour(ds, var):
    var_filled = np.empty_like(ds[var].values)
    for i in range(ds[var].shape[0]):  # iterate over time dimension
        for j in range(ds[var].shape[2]):  # iterate over lat dimension
            for k in range(ds[var].shape[3]):  # iterate over lon dimension
                var_slice = ds[var].values[i, :, j, k]
                var_filled[i, :, j, k] = nearest_neighbour(var_slice, ds['plev'].values)
    return var_filled

# horiztontal regridding; bilinear for continous variables and conservative for fluxes to ensure energy conservation
def cdo_bilinear_regridding(input_file, output_file):
    command = f'cdo remapbil,target_grid.txt {input_file} {output_file}'
    return os.system(command)

def cdo_conservative_regridding(input_file, output_file):
    command = f'cdo remapcon,target_grid.txt {input_file} {output_file}'
    return os.system(command)

# function to calculate the tropopause pressure at different latitudes using the method described in the ClimKern package:
# tropopause height of 10,000 Pa at the equator, descending with the cosine of latitude to 30,000 Pa at the poles
def make_tropo(da):
    p_trop = (3e4 - 2e4 * np.cos(np.deg2rad(da.lat))).broadcast_like(da)
    return p_trop

# function to perform nearest neighbour interpolation on a chunk of the dataset
def interpolate_chunk(chunk, plev_values):
    def interpolate(slice_values):
        finite_indices = np.isfinite(slice_values)
        if finite_indices.sum() > 0:  # ensure there are finite values to interpolate
            return griddata(plev_values[finite_indices], slice_values[finite_indices], plev_values, method="nearest", fill_value="extrapolate")
        else:
            return slice_values  # return original values if no finite values found

    # apply the interpolation function along the 'plev' axis (axis 1)
    for t in range(chunk.shape[0]):
        for y in range(chunk.shape[2]):
            for x in range(chunk.shape[3]):
                chunk[t, :, y, x] = interpolate(chunk[t, :, y, x])
    return chunk

# fucntion to apply the interpolation function to each chunk in parallel
def apply_nearest_neighbour_parallel(ds, var):
    var_values = ds[var].values
    plev_values = ds['plev'].values

    # specify chunk size and chunck the data
    chunk_size = 10  
    chunks = [var_values[i:i + chunk_size] for i in range(0, var_values.shape[0], chunk_size)]

    # process each chunk in parallel to fill missing values using nearest neighbour interpolation
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(interpolate_chunk, chunk, plev_values) for chunk in chunks]
        results = [future.result() for future in concurrent.futures.as_completed(futures)] # processed data

    # combine the results back into a single array
    var_filled = np.concatenate(results, axis=0)
    return var_filled

