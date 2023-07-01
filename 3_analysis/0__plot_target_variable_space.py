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

for perc in [0.5, 0.9]:
    
    # read data ===============================================================
    csv = pd.read_csv(f"/home/silvia/repos/predicting_fixation/1_datasets/simcomms/processed_data_simcomms_{perc}_full_jun", index_col="sample")
    
    prefix = str(decimal.Decimal(perc*100)) + "%_"
    
    csv = csv.dropna() # remove samples that have not even reached fixation in 1000 transfers
    
    # substitute 1.00E+06 with 1000000
    csv = csv.replace(to_replace="1.00E+06", value=100000)
    
    # factorize distrib
    csv.distrib = [1 if i=="uniform" else 0 for i in csv.distrib]

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
    for n, i in enumerate(['distrib',
                           'size',
                           'richness',
                           'dilfactor',
                           'shannon',
                           'evenness',
                           'gini']):
        name = ["Abundance distribution", "Total community size", "Richness", "Dilution factor", "Shannon diversity", "Evenness", "Gini index"][n]
        
        cmap = None
        mynorm = None
        ## Logarithmic color scale for these
        if i in ["richness", "dilfactor"]:
            mynorm = colors.LogNorm()
        ## Not a gradient needed for these
        if i in ["distrib", "size"]:
            cmap = colors.ListedColormap(["#541352FF", "yellow"])  # Binary colormap (viridis)
            bounds = [np.unique(csv[i]).min(), np.unique(csv[i]).max()]  # Boundaries for colormap
            
            plt.figure(figsize=(10, 8), dpi=80)
            plt.scatter(range(csv.shape[0]),
                        csv.success,
                        c=csv[i],
                        cmap=cmap,
                        norm=None,
                        facecolors="none",
                        alpha=0.5)
            cbar = plt.colorbar(ticks=bounds)
            if i == "distrib":
                cbar.ax.set_yticklabels(["Log-normal", "Uniform"])
            else:
                cbar.ax.set_yticklabels(["10000", "10⁶"])
        else:
            plt.figure(figsize=(10, 8), dpi=80)
            plt.scatter(range(csv.shape[0]),
                        csv.success,
                        c=csv[i],
                        norm=mynorm,
                        facecolors="none",
                        alpha=0.5)
            plt.colorbar()
        
    
        plt.title(name + " (" + str(decimal.Decimal(perc*100)) + "% threshold)")
        plt.ylabel("Dilution-transfer cycle")
        plt.xticks([])
        plt.xlabel("← Communities →")
        plt.savefig("../figures/0__target_variable_space_" + prefix + i + ".png")
        plt.show()
        plt.close()
