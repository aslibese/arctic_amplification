#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon May 20 2024

This script calculates the Arctic and Tropical warming contributions (K) of individual feedback mechanisms, ERF, OHT and AHT
for 1979-2023 and creates a scatter plot. 

@author: aslibese
"""

import xarray as xr
import numpy as np   
from FeedbackUtilities import weighted_avg
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import os

# load the warming contributions dataset
ds_path = 'data/CMIP6/MIROC6_warming_contr.nc'
warming_contr_ds = xr.open_dataset(ds_path)

# extract model name from the path
model = os.path.basename(ds_path).split('_')[0]

# define the boundaries for the Arctic and Tropics
# tropics (30°N–30°S) definition follow Pithan & Mauritsen (2014) and Soden et al. (2008)
arctic_lower = 66.5 
tropics_lower = -30.0
tropics_upper = 30.0

# calculate Arctic warming contributions; ck.spat_avg computes the spatial average while weighting for cos(latitude)
arctic_LR = weighted_avg(warming_contr_ds['LR'], arctic_lower)
arctic_Planck = weighted_avg(warming_contr_ds['Planck'], arctic_lower)
arctic_q_lw = weighted_avg(warming_contr_ds['q_lw'], arctic_lower)
arctic_q_sw = weighted_avg(warming_contr_ds['q_sw'], arctic_lower)
arctic_alb = weighted_avg(warming_contr_ds['alb'], arctic_lower)
arctic_cld_lw = weighted_avg(warming_contr_ds['cld_lw'], arctic_lower)
arctic_cld_sw = weighted_avg(warming_contr_ds['cld_sw'], arctic_lower)
arctic_ERF_lw = weighted_avg(warming_contr_ds['ERF_lw'], arctic_lower)
arctic_ERF_sw = weighted_avg(warming_contr_ds['ERF_sw'], arctic_lower)
arctic_OHU = weighted_avg(warming_contr_ds['OHU'], arctic_lower)
arctic_AHT = weighted_avg(warming_contr_ds['AHT'], arctic_lower)

# calculate Tropical warming contributions
tropics_LR = weighted_avg(warming_contr_ds['LR'], tropics_lower, tropics_upper)
tropics_Planck = weighted_avg(warming_contr_ds['Planck'], tropics_lower, tropics_upper)
tropics_q_lw = weighted_avg(warming_contr_ds['q_lw'], tropics_lower, tropics_upper)
tropics_q_sw = weighted_avg(warming_contr_ds['q_sw'], tropics_lower, tropics_upper)
tropics_alb = weighted_avg(warming_contr_ds['alb'], tropics_lower, tropics_upper)
tropics_cld_lw = weighted_avg(warming_contr_ds['cld_lw'], tropics_lower, tropics_upper)
tropics_cld_sw = weighted_avg(warming_contr_ds['cld_sw'], tropics_lower, tropics_upper)
tropics_ERF_lw = weighted_avg(warming_contr_ds['ERF_lw'], tropics_lower, tropics_upper)
tropics_ERF_sw = weighted_avg(warming_contr_ds['ERF_sw'], tropics_lower, tropics_upper)
tropics_OHU = weighted_avg(warming_contr_ds['OHU'], tropics_lower, tropics_upper)
tropics_AHT = weighted_avg(warming_contr_ds['AHT'], tropics_lower, tropics_upper)


# create a data dictionary
data = {
    "LR": {"x": tropics_LR, "y": arctic_LR, "color": "green"},
    "P'": {"x": tropics_Planck, "y": arctic_Planck, "color": "goldenrod"},
    "WV$_{LW}$": {"x": tropics_q_lw, "y": arctic_q_lw, "color": "blue"},
    "WV$_{SW}$": {"x": tropics_q_sw, "y": arctic_q_sw, "color": "cornflowerblue"},
    "A": {"x": tropics_alb, "y": arctic_alb, "color": "red"},
    "C$_{LW}$": {"x": tropics_cld_lw, "y": arctic_cld_lw, "color": "teal"},
    "C$_{SW}$": {"x": tropics_cld_sw, "y": arctic_cld_sw, "color": "darkturquoise"},
    "ERF$_{LW}$": {"x": tropics_ERF_lw, "y": arctic_ERF_lw, "color": "darkorange"},
    "ERF$_{SW}$": {"x": tropics_ERF_sw, "y": arctic_ERF_sw, "color": "sienna"},
    "OHU" : {"x": tropics_OHU, "y": arctic_OHU, "color" : "midnightblue"},
    "AHT" : {"x": tropics_AHT, "y": arctic_AHT, "color" : "magenta"},
}

# create figure and axis
fig, ax = plt.subplots()

# add lines
ax.plot([-2.0, 2.0], [-2.0, 2.0], color='lightgrey', linestyle='-', zorder=1)  # 1:1 line
ax.plot([-2.0, 1.0], [-1.0, 2.0], color='lightgrey', linestyle='--', zorder=1)  # 1:2 line

# plot data points 
texts = []
for label, coord in data.items():
    ax.scatter(coord["x"], coord["y"], label=label, color=coord["color"], s=100)
    
    if label == "LR":
        texts.append(ax.text(coord["x"] - 0.1, coord["y"], label, fontsize=12, ha='right', color=coord["color"])) 
    elif label == "P'":
        texts.append(ax.text(coord["x"] - 0.08, coord["y"], label, fontsize=12, ha='right', color=coord["color"])) 
    elif label == "WV$_{LW}$":
        texts.append(ax.text(coord["x"] - 0.1, coord["y"], label, fontsize=12, ha='right', color=coord["color"])) 
    elif label == "WV$_{SW}$":
        texts.append(ax.text(coord["x"] - 0.1, coord["y"], label, fontsize=12, ha='right', color=coord["color"])) 
    elif label == "A":
        texts.append(ax.text(coord["x"], coord["y"] + 0.13, label, fontsize=12, ha='center', color=coord["color"])) 
    elif label == "C$_{LW}$":
        texts.append(ax.text(coord["x"] + 0.1, coord["y"], label, fontsize=12, ha='left', color=coord["color"])) 
    elif label == "C$_{SW}$":
        texts.append(ax.text(coord["x"] - 0.1, coord["y"], label, fontsize=12, ha='right', color=coord["color"])) 
    elif label == "ERF$_{LW}$":
        texts.append(ax.text(coord["x"] + 0.1, coord["y"], label, fontsize=12, ha='left', color=coord["color"])) 
    elif label == "ERF$_{SW}$":
        texts.append(ax.text(coord["x"] + 0.1, coord["y"], label, fontsize=12, ha='left', color=coord["color"])) 
    elif label == "OHU":
        texts.append(ax.text(coord["x"], coord["y"] + 0.13, label, fontsize=12, ha='center', color=coord["color"])) 
    elif label == "AHT":
        texts.append(ax.text(coord["x"], coord["y"] - 0.3, label, fontsize=12, ha='center', color=coord["color"])) 

ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.xaxis.set_tick_params(top=True, labeltop=False)
ax.yaxis.set_tick_params(right=True, labelright=False)
ax.minorticks_on()

# determine the min and max values for x and y data
x_min = min(coord["x"] for coord in data.values())
x_max = max(coord["x"] for coord in data.values())
y_min = min(coord["y"] for coord in data.values())
y_max = max(coord["y"] for coord in data.values())

# set major and minor ticks 
ax.set_xticks(np.arange(-2.0, 2.5, 0.5))
ax.set_yticks(np.arange(-2.0, 2.5, 0.5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.minorticks_on()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')

# add x and y axis limit
ax.set_xlim(-2.0, 2.0)
ax.set_ylim(-2.0, 2.0)

# add gridlines
ax.set_axisbelow(True)
ax.grid(True, linestyle='--', linewidth=0.5)

# set labels and title
ax.set_xlabel('Tropical warming (K)')
ax.set_ylabel('Arctic warming (K)')
ax.set_title(f'Contributions to Arctic vs. Tropical Warming for 1979-2023 in {model}', pad=20)

# save plot
plt.savefig(f'figures/{model}_warming_contr.png', dpi=300)