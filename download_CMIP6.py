#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 2024

This driver script reads in the CMIP6 data from the Earth System Grid Federation (ESGF) archive based on specified variables. 

@author: aslibese
"""
from pyesgf.search import SearchConnection
import concurrent.futures

# function that performs a search for a specific variable and records all model and run combinations that contain the variable
def fetch_model_runs(variable, experiment_ids, project='CMIP6', frequency='mon', start_date='1979-01-01', end_date='2023-12-31'):
    conn = SearchConnection('https://esgf-node.llnl.gov/esg-search', distrib=True)
    model_runs = {}

    for experiment in experiment_ids:
        ctx = conn.new_context(
            project=project,
            experiment_id=experiment,
            frequency=frequency,
            variable=variable,
            start=start_date,
            end=end_date
        )
        
        for result in ctx.search():
            model_id = result.json['source_id'][0] if isinstance(result.json['source_id'], list) else result.json['source_id']
            member_id = result.json['variant_label'][0] if isinstance(result.json['variant_label'], list) else result.json['variant_label']

            # Ensure model_id is a string, especially if it comes as a list
            if isinstance(model_id, list):
                model_id = model_id[0]

            # Ensure member_id is a string, especially if it comes as a list
                if isinstance(member_id, list):
                    member_id = member_id[0]  # Taking the first element as the member_id

            # Ensure model_id is a string, especially if it comes as a list
                if model_id not in model_runs:
                    model_runs[model_id] = set()
                model_runs[model_id].add(member_id)

    return model_runs


# after collecting the results for each variable, the script finds the intersection of these results
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

variables = ['tas', 'tos', 'siconc', 'sftof']
experiments = ['historical', 'ssp585']
models_and_runs = search_cmip6_models(variables, experiments)

print("\nModels and their runs containing all specified variables:")
for model, runs in models_and_runs.items():
    print(f"Model: {model}, Runs: {list(runs)}")

    