#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 14:57:41 2023

@author: silvia

@description:
    diversity: color
    community type: shape
    ejes: dilfactor, success
    size: tamaño de objeto
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import decimal
import os
import matplotlib.markers as mmarkers
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
from matplotlib import colormaps as cmaps
from matplotlib.patches import Patch # legend labels


output_folder = "/home/silvia/repos/predicting_fixation/figures/0__target_variable_space/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

perc = 0.5
mydilfactors = [0.00025, 0.0005, 0.001, 0.0025, 0.004, 0.005, 0.008, 0.025]

# read data ===============================================================
csv = pd.read_csv(f"/home/silvia/repos/predicting_fixation/1_datasets/simulation_results/processed_data_simcomms_{perc}_full_jun")
mydilfactors = csv.dilfactor.unique() # DEBUG

prefix = str(decimal.Decimal(perc*100)) + "%_"

csv = csv.dropna() # remove samples that have not even reached fixation in 1000 transfers

# substitute 1.00E+06 with 1000000
csv = csv.replace(to_replace="1.00E+06", value=100000)

# factorize distrib
csv.distrib = [1 if i=="uniform" else 0 for i in csv.distrib]

# filter by dilution factor ===============================================
csv = csv[csv.dilfactor.isin(mydilfactors)]
csv.index = range(len(csv))
# important: fix indeces!!!
        
        
## plot ===================================================================
csv["community_type"] = csv.apply(lambda row: '_'.join([str(row['distrib']), str(row['size']), str(row['richness'])]), axis=1)
df_mapping = {}
for n, df in enumerate(mydilfactors):
    df_mapping[df] = n
    
csv['dilfactor2'] = csv['dilfactor'].replace(df_mapping)

# ## markers
# unique_values = csv['combined'].unique()
# marker_shapes = list(mmarkers.MarkerStyle.markers.keys())[:len(unique_values)]  # Adjust the marker shapes as needed
# mapping = dict(zip(unique_values, marker_shapes))

# plt.figure(figsize=(10, 8), dpi=80)

# for dun, du in enumerate(csv.combined.unique()):
#     selection = csv.loc[csv.combined==du]
#     plt.scatter(selection.dilfactor2,
#                 selection.success,
#                 c=selection.shannon,
#                 facecolors="none",
#                 s = selection["size"]**(2/5),
#                 marker=mapping[du],
#                 alpha=0.5)
# plt.colorbar()


# var1 = "distrib"
# var2 = "community_type"

var1 = "community_type"
var2 = "shannon"


color_palette = cmaps.get_cmap('tab10')

if var2 == "shannon":
    ## colors
    normalized_values = (csv[var2] - np.min(csv[var2])) / (np.max(csv[var2]) - np.min(csv[var2]))
    colors = [color_palette(value) for value in normalized_values]

    ## markers
    unique_values = csv[var1].unique()
    marker_shapes = list(mmarkers.MarkerStyle.markers.keys())[:len(unique_values)]  # Adjust the marker shapes as needed
    mapping = dict(zip(unique_values, marker_shapes))
    
    plt.figure(figsize=(14, 10), dpi=80)
    
    for dun, du in enumerate(csv[var1].unique()):
    
        selection = csv.loc[csv[var1]==du]
        
        plt.scatter(selection.dilfactor2,
                    selection.success,
                    c=[colors[i] for i in selection.index],
                    facecolors="none",
                    s = (1/2)*selection["size"]**(2/5),
                    marker=mapping[du],
                    alpha=0.5)


    
    plt.colorbar()

else:
    ## markers
    unique_values = csv[var1].unique()
    marker_shapes = list(mmarkers.MarkerStyle.markers.keys())[:len(unique_values)]  # Adjust the marker shapes as needed
    mapping = dict(zip(unique_values, marker_shapes))
    
    ## colors
    color_mapping = {}
    for us, string in enumerate(csv[var2].unique()):
        color_mapping[string] = color_palette(us % color_palette.N)
    colors = [mcolors.rgb2hex(color_mapping[i]) for i in csv[var2]]
        
    plt.figure(figsize=(14, 10), dpi=80)
    
    for dun, du in enumerate(csv[var1].unique()):
    
        selection = csv.loc[csv[var1]==du]
        
    
        plt.scatter(selection.dilfactor2,
                    selection.success,
                    c=[colors[i] for i in selection.index],
                    facecolors="none",
                    s = (1/3)*selection["size"]**(2/5),
                    marker=mapping[du],
                    alpha=0.5)
    
    # unique_labels = csv[var2].unique()
    # legend_handles = [Patch(color=color_mapping[s]) for s in unique_labels]
    # legend = plt.legend(handles=legend_handles, labels=unique_labels, loc='best', frameon=False)
    
    
plt.title(f"Markers shapes: {var1}\n\nColors: {var2}\n\nDot size: community size")
plt.ylabel("Dilution-transfer cycle")
plt.xticks(ticks=csv.dilfactor2.unique(), labels = csv.dilfactor.unique())
plt.xlabel("← Dilfactor →")    

for c in csv["sample"].unique():
    plt.plot(csv.loc[csv["sample"] == c]["dilfactor2"], csv.success[csv["sample"] == c], linestyle='-', marker='', alpha=0.5, linewidth=0.5)
    
plt.savefig("/home/silvia/test.png")
plt.show()
plt.close()


for c in csv["sample"].unique()[15:17]:
    plt.plot(csv.loc[csv["sample"] == c]["dilfactor2"], csv.success[csv["sample"] == c], linestyle='-', marker='', alpha=0.5, linewidth=0.5)