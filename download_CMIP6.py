#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 2024

This driver script reads in the CMIP6 data from the Earth System Grid Federation (ESGF) archive based on specified variables. 

@author: aslibese
"""
from pyesgf.search import SearchConnection

def search_cmip6_models(variables, experiment_ids, project='CMIP6', frequency='mon'):
    conn = SearchConnection('https://esgf-node.llnl.gov/esg-search', distrib=True)
    common_models = None

    for variable in variables:
        variable_model_runs = {}
        for experiment in experiment_ids:
            ctx = conn.new_context(
                project=project,
                experiment_id=experiment,
                frequency=frequency,
                variable=variable
            )
            print(f"Searching {experiment} for {variable}, found {ctx.hit_count} datasets.")

            for result in ctx.search():
                model_id = result.json['source_id'][0] if isinstance(result.json['source_id'], list) else result.json['source_id']
                member_id = result.json['variant_label'][0] if isinstance(result.json['variant_label'], list) else result.json['variant_label']

                if model_id not in variable_model_runs:
                    variable_model_runs[model_id] = set()
                variable_model_runs[model_id].add(member_id)

        if common_models is None:
            common_models = variable_model_runs
        else:
            # Intersect with the previous results
            common_models = {model: runs & common_models.get(model, set()) for model, runs in variable_model_runs.items() if model in common_models}

    return common_models

variables = ['tas', 'tos', 'sic', 'sftotf']
experiments = ['historical', 'ssp245']

models_and_runs = search_cmip6_models(variables, experiments)
print("\nModels and their runs containing all specified variables:")
for model, runs in models_and_runs.items():
    print(f"Model: {model}, Runs: {list(runs)}")


    