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

output_folder = "/home/silvia/repos/predicting_fixation/figures/0__target_variable_space/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for perc in [0.50, 0.90]:
    
    # read data ===============================================================
    csv = pd.read_csv(f"/home/silvia/repos/predicting_fixation/1_datasets/simulation_results/processed_data_simcomms_{perc}_full_jul")
    
    prefix = str(decimal.Decimal(perc*100)) + "%_"
    
    csv = csv.dropna() # remove samples that have not even reached fixation in 1000 transfers
    
    ## Plot ocurrences
    ## ============================================================
    csv['param_combination'] = csv.apply(lambda row: '_'.join(map(str, [row['distrib'], row['size'], row['richness'], row['dilfactor']])), axis=1)
    
    param_counts = csv['param_combination'].value_counts()
    csv['count'] = csv['param_combination'].map(param_counts)
    
    for n, i in enumerate(['count']):
        cmap = None
        mynorm = None
        
        plt.figure(figsize=(8, 8), dpi=100)
        
        scatter = plt.scatter(
            range(csv.shape[0]),
            csv.success,
            c=csv[i],
            cmap='cividis',
            norm=mynorm,
            facecolors="none",
            alpha=0.5
        )
        plt.colorbar(scatter, label='Número de comunidades con mismos parámetros', orientation='horizontal', pad=.05)
        # plt.title("Número de comunidades exitosas")
        plt.ylabel("Ciclo dilución-crecimiento", fontsize=13)
        plt.xticks([])
        plt.xlabel("← Comunidades →", fontsize=13)
        
        # Guardar la figura
        output_folder = "/home/silvia/repos/predicting_fixation/figures/0__target_variable_space/"
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, prefix + "success_colored_by_count.png"))
        plt.show()
