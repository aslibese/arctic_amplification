#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Jun 7 2024

This script serves to find the list of CMIP6 models and their ensemble members containing the specified 
variables both in historical and ssp585 runs and saves this list to a CSV file. 

@author: aslibese
"""

import intake
import pandas as pd

# load the CMIP6 data catalog
col_url = "https://storage.googleapis.com/cmip6/pangeo-cmip6.json"
col = intake.open_esm_datastore(col_url)

# define the variables and experiments 
variables = ['ta', 'hus', 'tas', 'ts', 'ps', 'rsdt', 'rsut', 'rlut', 'rsutcs', 'rlutcs', 'rsds', 'rsus', 'rlus', 'rlds', 'hfls', 'hfss']
experiments = ['historical', 'ssp585']

# query the catalog for datasets that match the criteria
query = col.search(experiment_id=experiments, variable_id=variables)

# convert to a DataFrame
df = query.df

# group by model and member_id, and filter those that have all required variables
required_vars_set = set(variables)
grouped = df.groupby(['source_id', 'member_id']).agg({'variable_id': lambda x: set(x)})
filtered = grouped[grouped['variable_id'].apply(lambda x: required_vars_set.issubset(x))]

# save the results to a CSV file
filtered.reset_index().to_csv('ensemble_members.csv', index=False)
