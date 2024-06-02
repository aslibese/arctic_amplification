#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 2024

This driver script reads in the CMIP6 data from the Earth System Grid Federation (ESGF) archive,
finds the CMIP6 models that contain the specified variables in the specified runs and their ensemble members,
and saves them in .csv file.

@author: aslibese
"""
from pyesgf.search import SearchConnection
import concurrent.futures
import pandas as pd

# performs a search for specific variables and records all model and run combinations that contain the variable
def fetch_model_runs(variable, experiment_ids, project='CMIP6', frequency='mon'):
    conn = SearchConnection('https://esgf-node.llnl.gov/esg-search', distrib=True)
    model_runs = {}

    for experiment in experiment_ids:
        ctx = conn.new_context(
            project=project,
            experiment_id=experiment,
            frequency=frequency,
            variable=variable,
        )
        
        for result in ctx.search():
            model_id = result.json['source_id'][0] if isinstance(result.json['source_id'], list) else result.json['source_id']
            member_id = result.json['variant_label'][0] if isinstance(result.json['variant_label'], list) else result.json['variant_label']

            if model_id not in model_runs:
                model_runs[model_id] = set()
            model_runs[model_id].add(member_id)

    return model_runs

# intersects these results to find common model IDs that contain all variables
def search_cmip6_models(variables, experiment_ids):
    variable_model_runs = {var: {} for var in variables}
    common_model_runs = {}

    # perform searches in parallel to speed up the process
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {}
        for variable in variables:
            future = executor.submit(fetch_model_runs, variable, experiment_ids)
            futures[future] = variable

        for future in concurrent.futures.as_completed(futures):
            variable_model_runs[futures[future]] = future.result()

    # Initialize common model_runs with the first variable's runs
    common_model_runs = variable_model_runs[variables[0]]

    # Intersect with the runs of other variables
    for variable in variables[1:]:
        new_common = {}
        for model in common_model_runs:
            if model in variable_model_runs[variable]:
                new_common[model] = common_model_runs[model] & variable_model_runs[variable][model]
        common_model_runs = new_common

    return common_model_runs


variables = ['ta', 'hus','tas', 'ts', 'ps','rsdt', 'rsut', 'rlut', 'rsutcs', 'rlutcs',
             'rsds', 'rsus', 'rlus', 'rlds', 'hfls', 'hfss']
experiments = ['historical', 'ssp585']
models_and_runs = search_cmip6_models(variables, experiments)

df = pd.DataFrame(models_and_runs.items(), columns=['Model', 'Runs'])

file_path = 'data/CMIP6/CMIP6_models_runs_list.csv'
df.to_csv(file_path, index=False)

print(f"Data saved to {file_path}")
