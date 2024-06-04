#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon May 20 2024

Utilities script containing the functions used in CalculateFeedbacks.py

@author: aslibese
"""

import numpy as np
import scipy.stats as st

# function to create monthly climatology
def make_clim(da):
    "Produce monthly climatology of model field."
    clim = (
        da.groupby(da.time.dt.month)
        .mean(dim="time", skipna=True)
        .rename({"month": "time"})
    )
    return clim

# function to regress radiative contributions of feedbacks (y) against 'ts' anomalies (x) in each grid cell 
# to quantify the change in radiative forcing per unit change in surface temperature: feedback parameter (λ)
def regress_against_ts(ts_anom, feedback):
    n_time, n_lat, n_lon = feedback.shape
    slopes = np.zeros((n_lat, n_lon))
    for i in range(n_lat):
        for j in range(n_lon):
            ts_anom_gridcell = ts_anom[:, i, j] 
            feedback_gridcell = feedback[:, i, j]
            slopes[i, j] = np.polyfit(ts_anom_gridcell, feedback_gridcell, 1)[0]
    return slopes
    

# Δ is the trend in each variable during 1979-2023, multiplied by the number of years (45)
def calculate_delta(y):
    years = np.arange(1979, 2024)
    slope, intercept, r_value, p_value, std_err = st.linregress(years, y)
    return slope * len(years)


# function to compute the spatial average while weighting for cos(latitude), adapted from the ClimKern package
def weighted_avg(data, lat_bound_s=-90, lat_bound_n=90):
    # constrain area of interest and take zonal mean
    data = data.sel(lat=slice(lat_bound_s,lat_bound_n)).mean(dim='lon')

    # compute weights
    weights = np.cos(np.deg2rad(data.lat))/np.cos(np.deg2rad(data.lat)).sum()
  
    # compute average
    avg = (data * weights).sum(dim='lat')

    # return avg: new spatially averaged DataArray with lat and long coordinates now removed
    return avg
