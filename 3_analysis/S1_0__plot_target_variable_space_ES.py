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

output_folder = "/home/silvia/repos/predicting_fixation/figures/0__target_variable_space/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for perc in [0.50, 0.90]:
    
    # read data ===============================================================
    csv = pd.read_csv(f"/home/silvia/repos/predicting_fixation/1_datasets/simulation_results/processed_data_simcomms_{perc}_full_jul")
    
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
        name = ["Distribución de abundancias", "Tamaño total de la comunidad", "Riqueza", "Factor de dilución", "Diversidad de Shannon", "Homogeneidad de Pielou", "Índice de Gini"][n]
        
        cmap = None
        mynorm = None
        ## Logarithmic color scale for these
        if i in ["richness", "dilfactor"]:
            mynorm = colors.LogNorm()
        
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
                cbar.ax.set_yticklabels(["Log-normal", "Uniforme"])
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
        
    
        plt.title(name + " (umbral: " + str(decimal.Decimal(perc*100)) + "%)")
        plt.ylabel("Ciclo dilución-crecimiento")
        plt.xticks([])
        plt.xlabel("← Comunidades →")
        plt.savefig(output_folder + prefix + i + ".png")
        plt.show()
        plt.close()

    ## Plot together every 30 replicates of a single community type
    ## ============================================================
    csv['combined'] = csv.apply(lambda row: '_'.join([str(row['distrib']), str(row['size']), str(row['richness']), str(row['dilfactor'])]), axis=1)

    unique_strings = np.unique(csv["combined"])
    color_mapping = {}
    color_palette = plt.cm.get_cmap('tab20')
    for us, string in enumerate(unique_strings):
        color_mapping[string] = color_palette(us % color_palette.N)
    my_c = [color_mapping[s] for s in list(csv["combined"])]

    plt.figure(figsize=(10, 8), dpi=80)
    plt.scatter(range(csv.shape[0]),
                csv.success,
                c=my_c, # use mapping
                facecolors="none",
                alpha=0.5)
    plt.colorbar()
    
    plt.title("all 30 replicates" + " (" + str(decimal.Decimal(perc*100)) + "% threshold)")
    plt.ylabel("Dilution-transfer cycle")
    plt.xticks([])
    plt.xlabel("← Communities →")
    plt.savefig(output_folder + "/Combined/" + prefix + "Combined" + ".png")
    plt.show()
    plt.close()
    
    ## Same but each 30-replicate set generates a separate plot
    for us, string in enumerate(unique_strings):
        selection = csv.loc[csv["combined"]==string]
        my_c = [color_mapping[s] for s in list(selection["combined"])]
        plt.figure(figsize=(10, 8), dpi=80)
        plt.scatter(selection.index,
                    selection.success,
                    c=my_c, # use mapping
                    facecolors="none",
                    alpha=0.5)
        plt.colorbar()
        
        plt.title("Sample" + " (" + str(decimal.Decimal(perc*100)) + "% threshold)")
        plt.ylabel("Dilution-transfer cycle")
        plt.xticks([])
        plt.ylim(0, 1000)
        plt.xlim(0, len(csv.index))
        plt.xlabel("← Communities →")
        plt.savefig(output_folder + "/Combined/" + prefix + string + ".png")
        plt.show()
        plt.close()
        

    ## Now we make a subplot-plot for each , each subplots being a unique
    ## possible value 
    ## ============================================================
    variables = ['distrib', 'richness', 'size', 'dilfactor']
    unique_values = {var: csv[var].unique() for var in variables}
    
    # Loop through each variable
    for var in variables:
        # Get unique values for the current variable
        values = unique_values[var]
    
        # Calculate the number of rows and columns for the grid based on the number of unique values
        num_rows = int(len(values) / 3) + 1  # Adjust the number of columns as per your preference
        num_cols = min(len(values), 3)
    
        # Create a new figure and subplots grid
        if var == "dilfactor":
            fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, 28), dpi=260)
            plot_legend = True
        else:
            fig, axs = plt.subplots(num_rows, num_cols, figsize=(12, 7), dpi=180)
            plot_legend = False
        # Flatten the subplots grid
        axs = axs.flatten()
    
        # Loop through each unique value and create a scatter plot in a subplot
        for i, value in enumerate(values):
            # Filter the DataFrame based on the current variable and value
            selection = csv.loc[csv[var] == value]
    
            # Create a scatter plot in the current subplot
            for dun, du in enumerate(selection.distrib.unique()):
                selection2 = selection.loc[selection.distrib==du]                
                my_c = [color_mapping[s] for s in selection2["combined"]] # Create a color list for the scatter plot using the 'combined' column
                marker = [".", ","][dun] # And, choose marker based on the 'distrib' column
                axs[i].scatter(selection2.index, selection2['success'], c=my_c, facecolors='none', alpha=0.5, marker=marker)
            axs[i].set_title(f'{var} = {value}')
            axs[i].set_ylim(0, 1000)
            
            # Legend
            unique_labels = list(selection["combined"].unique())
            legend_handles = [Patch(color=color_mapping[s]) for s in unique_labels]
            if plot_legend:
                legend = axs[i].legend(handles=legend_handles, labels=unique_labels, loc='best', frameon=False)
                legend.get_frame().set_alpha(0.0)  # Make the background of the legend box invisible
                legend.get_texts()[0].set_fontsize(8)  # Set the font size for the legend labels


        # Remove any unused subplots
        if len(values) < num_rows * num_cols:
            for j in range(len(values), num_rows * num_cols):
                fig.delaxes(axs[j])
        
        # # Legend
        # legend_handles = [Patch(color=color_mapping[s]) for s in color_mapping]
        # legend_labels = list(color_mapping.keys())
        # fig.legend(handles=legend_handles, labels=legend_labels, loc = "lower center")
        
        fig.tight_layout() # Adjust the spacing between subplots
        plt.savefig(output_folder + "/Combined/" + prefix + str(var) + ".png")  # Change the file extension to .pdf for PDF files
    
    



    


    ## Save the legend as its own file
    ## ===============================
    # from matplotlib.patches import Patch # plotting the legend
    # fig_legend = plt.figure(figsize=(4, 18))
    # ax_legend = fig_legend.add_subplot(111)
    # # Create a list of proxy artists for the legend
    # legend_elements = [Patch(facecolor=color_mapping[key], label=key) for key in color_mapping]
    # ax_legend.legend(handles=legend_elements, title="Legend")
    # ax_legend.axis('off') # Remove ticks and spines
    # fig_legend.savefig(output_folder + "/Combined/LEGEND.png")
        
    ## Sample plots (how each unique replicate changes over time)
    ## ============
    # unique_strings = np.unique(csv["sample"])
    # color_mapping = {}
    # color_palette = plt.cm.get_cmap('tab20')
    # for us, string in enumerate(unique_strings):
    #     color_mapping[string] = color_palette(us % color_palette.N)
    # for us, string in enumerate(unique_strings):
    #     selection = csv.loc[csv["sample"]==string]
    #     my_c = [color_mapping[s] for s in list(selection["sample"])]
    #     plt.figure(figsize=(10, 8), dpi=80)
    #     plt.scatter(range(selection.shape[0]),
    #                 selection.success,
    #                 c=my_c, # use mapping
    #                 facecolors="none",
    #                 alpha=0.5)
    #     plt.colorbar()
        
    #     plt.title("Sample" + " (" + str(decimal.Decimal(perc*100)) + "% threshold)")
    #     plt.ylabel("Dilution-transfer cycle")
    #     plt.xticks([])
    #     plt.xlabel("← Communities →")
    #     plt.savefig(output_folder + prefix + string + ".png")
    #     plt.show()
    #     plt.close()
