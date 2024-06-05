import pandas as pd
import xarray as xr
import numpy as np
import gcsfs
import os
import warnings
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
    # combined_ds = xr.merge(datasets)
    combined_ds = xr.merge(datasets, compat='override')
    return combined_ds

# function to use CDO for vertical interpolation
def cdo_vertical_interpolation(input_file, output_file, target_plev):
    # join the list of pressure levels into a comma-separated string
    target_plevs_str = ','.join(map(str, target_plev))
    command = f'cdo intlevel,{target_plevs_str} {input_file} {output_file}'
    return os.system(command)

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