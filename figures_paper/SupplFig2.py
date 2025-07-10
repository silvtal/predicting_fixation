#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 30 14:39:05 2023

@author: silvia
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import decimal
import os
from matplotlib.patches import Patch # legend labels
import colorsys # version B
import random # randomizar la paleta de colores

output_folder = "/home/silvia/repos/predicting_fixation/figures_paper"

# color
color_palette = plt.cm.get_cmap('tab20')
color_palette = [color_palette(i) for i in range(color_palette.N)]
random.shuffle(color_palette)
color_palette = colors.ListedColormap(color_palette)
    
for perc in [0.50, 0.90]:
    
    # read data ===============================================================
    csv = pd.read_csv(f"/home/silvia/repos/predicting_fixation/1_datasets/simulation_results/processed_data_simcomms_{perc}_full_jul")
    
    prefix = str(decimal.Decimal(perc*100)) + "%_"
    
    csv = csv.dropna() # remove samples that have not even reached fixation in 1000 transfers
    csv = csv.reset_index(drop=True) #> para que no pete después

    # substitute 1.00E+06 with 1000000
    csv = csv.replace(to_replace="1.00E+06", value=100000)
    
    # suitability of dilution factors =========================================
    for fd in np.sort(np.unique(csv.dilfactor)):
        successes = list(csv.success[csv.dilfactor==fd])
        print(fd, "-- total communities:", str(len(successes)), " -- mean:", str(round(np.mean(successes), 3)))
            # 0.00025 -- total communities: 125  -- mean: 49.432
            # 0.0004 -- total communities: 334  -- mean: 56.973
            # 0.0005 -- total communities: 330  -- mean: 90.194
            # 0.001 -- total communities: 329  -- mean: 175.152
            # 0.0025 -- total communities: 229  -- mean: 77.131
            # 0.004 -- total communities: 201  -- mean: 90.02
            # 0.005 -- total communities: 201  -- mean: 110.433
            # 0.008 -- total communities: 195  -- mean: 176.021
            # 0.01 -- total communities: 195  -- mean: 222.015
            # 0.025 -- total communities: 96  -- mean: 145.417
            # 0.04 -- total communities: 88  -- mean: 196.705
            # 0.05 -- total communities: 89  -- mean: 231.472
            # 0.1 -- total communities: 66  -- mean: 420.364
            # 0.25 -- total communities: 8  -- mean: 719.625
            
            
    # make plots ==============================================================

    ## Plot together every 30 replicates of a single community type
    ## ============================================================
    csv['combined'] = csv.apply(lambda row: '_'.join([str(row['distrib']), str(row['size']), str(row['richness'])]), axis=1)

    unique_strings = np.unique(csv["combined"])

    # fill AND border for each point
    color_mapping_bright = {}
    color_mapping_dark = {}
    
    for us, string in enumerate(unique_strings):
        base_color = color_palette(us % color_palette.N)
        r, g, b = base_color[:3]
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        darker_rgb = colorsys.hls_to_rgb(h, max(0, l * 0.6), s)
    
        color_mapping_bright[string] = base_color[:3]  # relleno con alpha fijo
        color_mapping_dark[string] = (*darker_rgb, 1.0)        # borde más oscuro
    
    colors_edge = [color_mapping_dark[c] for c in csv['combined']]
    # colors_bright = [color_mapping_bright[c] for c in csv['combined']] # para el alpha fijo
    grouped = csv.groupby('combined')['dilfactor'] #> como el grouping de R
    min_dil = grouped.transform('min'); max_dil = grouped.transform('max')
    normalized_alpha = (csv['dilfactor'] - min_dil) / (max_dil - min_dil + 1e-9) #> Evitar división por cero (si todos los valores de dilfactor son iguales)
    scaled_alpha = 0.3 + 0.7 * normalized_alpha #> importante: que los alpha estén entre 0.3 y 1.0 para que se vean todos
    colors_bright = [
        (*color_mapping_bright[csv.loc[i, 'combined']][:3], scaled_alpha[i]) for i in range(len(csv))
        
    ]

    legend_elements = []
    for comb in unique_strings:
        color = color_mapping_bright[comb]
        patch = Patch(color=color, label=comb)
        legend_elements.append(patch)

    plt.figure(figsize=(8.25, 8), dpi=100)
    plt.scatter(
        range(csv.shape[0]),
        csv.success,
        c=colors_bright,
        edgecolors=colors_edge,
        linewidths=0.5,
        s=40
    )

    # plt.title("Replicates of simulated communities (fixation threshold of " + str(decimal.Decimal(perc*100)) + "%)", fontsize=14)
    plt.ylabel("Dilution-growth cycles until fixation", fontsize=13)
    plt.xlabel("\u2190 Simulated communities \u2192", fontsize=13)
    plt.xticks([])

    plt.legend(handles=legend_elements,
    title="Distribution_Size_Richness",
    # loc='upper left',
    bbox_to_anchor=(0.48, -0.3),
    ncol=3,
    loc="lower center",
    title_fontsize=13,
    fontsize=11,
    frameon=True)
    
    plt.tight_layout()
    plt.savefig(output_folder + "/" + prefix + "SupplFig2.png")
    plt.show()
    plt.close()


