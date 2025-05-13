#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import decimal
import os

# Leer datos
csv = pd.read_csv("/home/silvia/repos/predicting_fixation/1_datasets/simulation_results/processed_data_simcomms_0.5_full_jul")
csv = csv.dropna()
csv = csv.replace(to_replace="1.00E+06", value=100000)
csv['distrib'] = [1 if i == "uniform" else 0 for i in csv['distrib']]

# Crear una columna combinada con los parámetros relevantes
csv['param_combination'] = csv.apply(lambda row: '_'.join(map(str, [row['distrib'], row['size'], row['richness'], row['dilfactor']])), axis=1)

# Contar ocurrencias de cada combinación
param_counts = csv['param_combination'].value_counts()
csv['count'] = csv['param_combination'].map(param_counts)

# Crear el scatter plot
for n, i in enumerate(['count']):
    
    cmap = None
    mynorm = None
    
    plt.figure(figsize=(10, 8), dpi=80)
    plt.scatter(range(csv.shape[0]),
                csv.success,
                c=csv[i],
                cmap='viridis_r',
                norm=mynorm,
                facecolors="none",
                alpha=0.5)
    plt.colorbar()
    
plt.colorbar(scatter, label='Número de comunidades con mismos parámetros')
plt.title("Número de comunidades exitosas")
plt.ylabel("Ciclo dilución-crecimiento")
plt.xticks([])
plt.xlabel("← Comunidades →")

# Guardar la figura
output_folder = "/home/silvia/repos/predicting_fixation/figures/0__target_variable_space/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
plt.savefig(os.path.join(output_folder, "success_colored_by_count.png"))
plt.show()

