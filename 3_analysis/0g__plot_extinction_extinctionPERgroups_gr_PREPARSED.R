# HAY OTRO SCRIPT PARA LOS 'PARSED', LOS GENERADOS A PARTIR DE 'generate_data(...)'

# 'PREPARSED' : datos mejor separados y tal
# generated with ~/repos/predicting_fixation/3_analysis/g_plot_results_facet_by_dil.r
# No hace falta tener incluido en el df el número de ciclo en el que se ha dado
# el éxito, así que se pueden usar tablas donde no se sepa.

library(tidyverse)
library(flexplot)
#> run from script's directory !!
setwd("~/repos/predicting_fixation/3_analysis/")

# OPTIONS ----------------------------------------------------------------------
#> only 50% fixation threshold (% abundance for an OTU to be fixated)
threshold <- '0.50'
# threshold <- '0.30'
#> 95% of iterations needed (% iterations that need to be successful in order to consider the community successful)
percN <- '0.95'
# percN <- '0.90'

# no of dil-transf cycles (for the y axis limits)
n_cycles <- 200

#> different output folder ("_groups")
out_folder = "../figures_groups/"
out_folder <- paste0(out_folder, "/p", percN, "_f", threshold, "/")
if (!file.exists(out_folder)) {system(paste("mkdir -p", out_folder))}

threshold <- as.numeric(threshold); percN <- as.character(percN) #> numeric
prefix <- paste0(threshold*100, "_")

plot_bars <- TRUE #> whether to plot niche size bars at the background

# LOAD DATA --------------------------------------------------------------------
# Read the final processed table (~/repos/predicting_fixation/3_analysis/0g__plot_fixation_fixationPERgroups_gr_PREPARSED.R)
# Like "all_information.csv" but with two extra cols
csv_groups <- read_csv(paste0("group_PREPARSED_analysis_table_p", percN, "_f", threshold, ".csv"))

# PLOT - fixation --------------------------------------------------------------
#> EvenGroups; GROUP fixation (functional group no. vs. dilfactor effect)---
# we're going to plot the frequency of every nicheN-Percentage of Groups with Success combination.
# (!) "group_extinctions" es en csv_groups (comprobar aún así) un string con ";" separando ! 
# porcentaje de los grupos que sí tienen éxito

if (plot_bars) {
  df_long <- csv_groups %>%
    mutate(group_extinctions = strsplit(group_extinctions, ";")) %>%  # Separar los valores
    unnest(group_extinctions) %>%  # Expandir filas
    mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convertir a numérico
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
    mutate(group_number = row_number()) %>%  # Asignar el número de grupo
    mutate(abundances = 100 * str_split(abundances, ";")[[1]][group_number] %>% as.numeric) %>% 
    ungroup()
  
  df_summary <- df_long %>% # Calcular la media de los valores acumulativos para cada número de grupos
    group_by(nicheN, group_number, dilfactor, nichedist, richness,
             abundances) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
    summarise(mean_success = mean(group_extinctions), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
  
} else {
  df_long <- csv_groups %>%
    mutate(group_extinctions = strsplit(group_extinctions, ";")) %>%  # Separar los valores
    unnest(group_extinctions) %>%  # Expandir filas
    mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convertir a numérico
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
    mutate(cumulative_count = row_number()) %>%  # Asignar el número acumulado dentro de cada sample
    ungroup()
  
  df_summary <- df_long %>% # Calcular la media de los valores acumulativos para cada número de grupos
    group_by(nicheN, cumulative_count, dilfactor, nichedist, richness) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
    summarise(mean_success = mean(group_extinctions), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
}

# plot(s)
for (nichedistvalue in c("EvenGroups", "SkewedGroups")){
  for (richnessvalue in c(100, 1000)){
    df_summary_filtered <- df_summary[((df_summary$nichedist == nichedistvalue) &
                                        (df_summary$richness == richnessvalue)), ] %>%
      mutate(nicheN = as.factor(nicheN))
    
    # plots[[nichedistvalue]] <- ggplot(df_summary_filtered, aes(x = group_number, y = mean_success, color = as.factor(nicheN))) +
    p <- ggplot(df_summary_filtered, aes(x = group_number, y = mean_success, color = nicheN, fill = nicheN)) +
      geom_point(size = 1) +
      labs(x = "Número de grupo",
           y = "Tasa de extinción (%) / Tamaño de nicho (%)",
           color = "Número de grupos (nicheN)") +
      theme_minimal() +
      # facet_wrap(~nicheN, scales = "free_x") +
      # facet_grid(dilfactor ~ nichedist + nicheN, 
      # facet_wrap(. ~ nichedist + nicheN + dilfactor,
      facet_wrap(
                 . ~ dilfactor + nicheN,
                 axis.labels = "margins",
                 strip.position = "top",
                 scales = "free_x",
                 # ncol = 2,
                 labeller = labeller(nicheN = function(x) "", dilfactor = label_value) # otros: label:both
                 ) +
      theme(strip.text.x = element_text(angle = 0, vjust = -10, hjust = 1)) +
      coord_cartesian(ylim = c(0, 100)) +
      # theme(legend.position = "bottom") +
      theme(
        legend.position = "none",
        panel.grid.major = element_line(color = "lightgray", size = 0.25),  # Major grid lines
        panel.grid.minor = element_line(color = "gray90", size = 0.1)    # Minor grid lines
      ) +
      # scale_color_brewer(palette = "Paired", direction = -1) # https://www.datanovia.com/en/wp-content/uploads/dn-tutorials/ggplot2/figures/029-r-color-palettes-rcolorbrewer-palettes-1.png
      # scale_color_manual(values = c("#A6CEE3", "#B2DF8A", "#1F78B4", "#33A02C", "#FB9A99", "#FDBF6F", "#E31A1C", "#FF7F00"))
      scale_color_manual(values = c("#E31A1C", "#33A02C")) +
      scale_fill_manual(values = c("#FB9A99", "#B2DF8A")) +
      scale_x_continuous(breaks = function(x) seq(0, max(x), by = 1))  # Solo números naturales hasta nicheN
    
    if (plot_bars) {
      p <- p + 
        geom_bar(aes(x = group_number, y = abundances),
                 stat = "identity", position = "stack", alpha = 0.2,
                 color = NA)
    }

    ggsave(
      filename = paste0(out_folder,"/0g_extinction_per_group_", nichedistvalue, "_", richnessvalue, ".png"),
      plot = p,
      height = 1.2 * 2000, width = 2 * 1100,
      units = "px"
    )
  }
}
