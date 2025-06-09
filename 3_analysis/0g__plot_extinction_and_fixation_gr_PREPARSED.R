# CURSOR MADE
# Script that combines both extinction and fixation data
# Based on 0g__plot_extinction_extinctionPERgroups_gr_PREPARSED.R and 0g__plot_fixation_fixationPERgroups_gr_PREPARSED.R

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
# Read both extinction and fixation data
csv_groups_extinction <- read_csv(paste0("group_PREPARSED_analysis_table_p", percN, "_f", threshold, ".csv"))
csv_groups_fixation <- read_csv(paste0("group_PREPARSED_analysis_table_p", percN, "_f", threshold, ".csv"))

# PLOT - extinction and fixation --------------------------------------------------------------
#> Process extinction data
if (plot_bars) {
  df_long_extinction <- csv_groups_extinction %>%
    mutate(group_extinctions = strsplit(group_extinctions, ";")) %>%  # Separar los valores
    unnest(group_extinctions) %>%  # Expandir filas
    mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convertir a numérico
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
    mutate(group_number = row_number()) %>%  # Asignar el número de grupo
    mutate(abundances = 100 * str_split(abundances, ";")[[1]][group_number] %>% as.numeric) %>% 
    ungroup()
  
  df_summary_extinction <- df_long_extinction %>% # Calcular la media de los valores acumulativos para cada número de grupos
    group_by(nicheN, group_number, dilfactor, nichedist, richness,
             abundances) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
    summarise(mean_success = mean(group_extinctions), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
  
} else {
  df_long_extinction <- csv_groups_extinction %>%
    mutate(group_extinctions = strsplit(group_extinctions, ";")) %>%  # Separar los valores
    unnest(group_extinctions) %>%  # Expandir filas
    mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convertir a numérico
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
    mutate(cumulative_count = row_number()) %>%  # Asignar el número acumulado dentro de cada sample
    ungroup()
  
  df_summary_extinction <- df_long_extinction %>% # Calcular la media de los valores acumulativos para cada número de grupos
    group_by(nicheN, cumulative_count, dilfactor, nichedist, richness) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
    summarise(mean_success = mean(group_extinctions), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
}

# Process fixation data
if (plot_bars) {
  df_long_fixation <- csv_groups_fixation %>%
    mutate(group_success = strsplit(group_success, ";")) %>%  # Separar los valores
    unnest(group_success) %>%  # Expandir filas
    mutate(group_success = 100 * as.numeric(group_success)) %>%  # Convertir a numérico
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
    mutate(group_number = row_number()) %>%  # Asignar el número de grupo
    mutate(abundances = 100 * str_split(abundances, ";")[[1]][group_number] %>% as.numeric) %>% 
    ungroup()
  
  df_summary_fixation <- df_long_fixation %>% # Calcular la media de los valores acumulativos para cada número de grupos
    group_by(nicheN, group_number, dilfactor, nichedist, richness,
             abundances) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
    summarise(mean_success = mean(group_success), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
  
} else {
  df_long_fixation <- csv_groups_fixation %>%
    mutate(group_success = strsplit(group_success, ";")) %>%  # Separar los valores
    unnest(group_success) %>%  # Expandir filas
    mutate(group_success = 100 * as.numeric(group_success)) %>%  # Convertir a numérico
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
    mutate(group_number = row_number()) %>%  # Asignar el número de grupo
    ungroup()
  
  df_summary_fixation <- df_long_fixation %>% # Calcular la media de los valores acumulativos para cada número de grupos
    group_by(nicheN, group_number, dilfactor, nichedist, richness) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
    summarise(mean_success = mean(group_success), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
}

# plot(s)
for (nichedistvalue in c("EvenGroups", "SkewedGroups")){
  for (richnessvalue in c(100, 1000)){
    df_summary_filtered_extinction <- df_summary_extinction[((df_summary_extinction$nichedist == nichedistvalue) &
                                        (df_summary_extinction$richness == richnessvalue)), ] %>%
      mutate(nicheN = as.factor(nicheN))
    
    df_summary_filtered_fixation <- df_summary_fixation[((df_summary_fixation$nichedist == nichedistvalue) &
                                        (df_summary_fixation$richness == richnessvalue)), ] %>%
      mutate(nicheN = as.factor(nicheN))
    
    p <- ggplot() +
      # Plot fixation data in light gray
      geom_point(data = df_summary_filtered_fixation,
                aes(x = group_number, y = mean_success),
                color = "gray80", size = 1) +
       # Plot extinction data in color
      geom_point(data = df_summary_filtered_extinction, 
                aes(x = group_number, y = mean_success, color = nicheN, fill = nicheN),
                size = 1, shape = 1) +

      labs(x = "Número de grupo",
           y = "Tasa de extinción o fijación (%) / Tamaño de nicho (%)",
           color = "Número de grupos (nicheN)") +
      theme_minimal() +
      facet_wrap(
        . ~ dilfactor + nicheN,
        axis.labels = "margins",
        strip.position = "top",
        scales = "free_x",
        labeller = labeller(nicheN = function(x) "", dilfactor = label_value)
      ) +
      theme(strip.text.x = element_text(angle = 0, vjust = -10, hjust = 1)) +
      coord_cartesian(ylim = c(0, 100)) +
      theme(
        legend.position = "none",
        panel.grid.major = element_line(color = "lightgray", size = 0.25),
        panel.grid.minor = element_line(color = "gray90", size = 0.1)
      ) +
      scale_color_manual(values = c("#E31A1C", "#33A02C")) +
      scale_fill_manual(values = c("#FB9A99", "#B2DF8A")) +
      scale_x_continuous(breaks = function(x) seq(0, max(x), by = 1))  # Solo números naturales hasta nicheN
    
    if (plot_bars) {
      p <- p + 
        geom_bar(data = df_summary_filtered_extinction,
                 aes(x = group_number, y = abundances, fill = nicheN),
                 stat = "identity", position = "stack", alpha = 0.2,
                 color = NA)
    }

    ggsave(
      filename = paste0(out_folder,"/0g_extinction_and_fixation_per_group_", nichedistvalue, "_", richnessvalue, ".png"),
      plot = p,
      height = 1.2 * 2000, width = 2 * 1100,
      units = "px"
    )
  }
} 