# -*- coding: utf-8 -*-

# @author: rantanem
# script to fix Berkeley Earth data time axis

import sys
import pandas as pd
import xarray as xr
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='A program that requires two data file names: input path and output path.')
parser.add_argument('filename', type=str, nargs=2, help='Path of the input and output files')

args = parser.parse_args()

inputpath = args.filename[0]
outputpath = args.filename[1]

# Open the file as a dataset
ds = xr.open_dataset(inputpath)


# Create a new time axis which better follows the netCDF standards
years = np.arange(1850,2024)

t_axis = pd.date_range(str(years[0])+'-01-01',pd.Timestamp(str(years[-1])+'-12-31'), freq='1M')


# Change the time variable inside the dataset and rewrite to a temp file

ds['time'] = t_axis
ds.to_netcdf(outputpath, format='NETCDF4', encoding={'time': {'dtype': 'i4'}})