#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat May 11 2024

This utilities script contains the functions serving to process the data:
select the period of interest, compute annual mean of temperature anomalies, select the Arctic region, 
perform linear regression on annual mean temp anomalies for each grid cell, 
compute the weighted average of the annual mean temperature anomalies. 

@author: aslibese
"""

import numpy as np   
import scipy.stats as st
import xarray as xr 

# function to slice the dataset based on the period of interest  
def selectPeriod(ds, start_year, end_year):
	# slice it based on start and end years
	ds_period = ds.sel(time = slice(str(start_year),str(end_year)))
	return ds_period

# function to find variables by dimensions
def find_variable_by_dims(ds, acceptable_dims):
	for var_name, variable in ds.data_vars.items():
		var_dims = set(variable.dims)
		if any(var_dims == dims for dims in acceptable_dims):
			return var_name
	return None

# function to compute the annual-mean of temperature anomalies 
def annualMean(ds, var_name):
	# make a copy
	ds_annualMean = ds.copy()
	# compute the annual-mean of tempanomaly
	ds_annualMean['annual_mean_tempanomaly'] = ds_annualMean[var_name].groupby('time.year').mean('time')
	return ds_annualMean

# function to slice the Arctic dataset based on the lower boundary 
def selectRegion(ds, lower_boundary):
	# determine the correct latitude dimension name
	lat_dim = 'lat' if 'lat' in ds.dims else 'latitude'
	ds_arctic = ds.sel(**{lat_dim: slice(float(lower_boundary), 90)})
	return ds_arctic

# custom function to perform linear regression using st.lingress
def performLinearRegression(y):
	x = np.arange(len(y))
	# perform linear regression
	slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
	return slope, intercept, p_value 

# function to perform linear regression for each grid cell and save slope and p_values to the dataset
def gridCellRegression(ds):
	# use xr.apply_ufunc() to apply our custom function to xarray objects
	result = xr.apply_ufunc(
		performLinearRegression, # custom function to apply 
		ds['annual_mean_tempanomaly'], # input: annual mean temperature anomalies
		vectorize = True, # apply function to each point in the dataset independently
		dask = 'parallelized', # use Dask for parallel computing if the data is chunked, to speed up computation
		input_core_dims = [['year']], # 'year' as the core dimension for the input
		output_core_dims = [[], [], []], # the outputs ('slope' and 'p_value') will have no core dimensions (it's scalar)
		output_dtypes = ['float64','float64','float64'] # data types of the output
		# maybe add .type()
	)

	# assign the results back to the dataset
	ds['slope'], ds['intercept'], ds['p_value'] = result

	return ds 


# function to compute the weighted average of annual mean temperature anomaly 
def weightedAverage(ds):
	# Determine the correct latitude dimension name
	lat_dim = 'lat' if 'lat' in ds.dims else 'latitude'
	lon_dim = 'lon' if 'lon' in ds.dims else 'longitude'

	annual_mean_tempanomaly = ds['annual_mean_tempanomaly']

	# weight each grid cell by its area to have an accurate representation of how much each grid cell contributes to the area average
	weights = np.cos(np.deg2rad(ds[lat_dim]))
	weights /= weights.mean() # normalizing weights so they sum to the total number of latitude points

	weighted_tempanomaly = annual_mean_tempanomaly * weights

	# apply the weights to the annual mean temperature anomaly
	# sum over the lat and lon dimensions to get a weighted global average
	# normalize it by the total weight to get a meaningful average
	weighted_sum = weighted_tempanomaly.sum(dim=[lat_dim, lon_dim])
	total_weights = weights.sum(dim=lat_dim) * len(ds[lon_dim])  # assuming lon_dim is evenly spaced and fully covered
	weighted_average = weighted_sum / total_weights

	return weighted_average








