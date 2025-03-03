# Comparamos la tasa de extinción por grupo con y sin interacciones de varios tipos.
# 'PREPARSED' : datos mejor separados y tal
# generated with ~/repos/predicting_fixation/3_analysis/g_plot_results_facet_by_dil.r
# No hace falta tener incluido en el df el número de ciclo en el que se ha dado
# el éxito, así que se pueden usar tablas donde no se sepa.

library(tidyverse)
library(flexplot)
#> run from script's directory !!
setwd("~/repos/predicting_fixation/3_analysis/")

# FIXED OPTIONS ----------------------------------------------------------------
#> only 50% fixation threshold (% abundance for an OTU to be fixated)
threshold <- '0.50'
# threshold <- '0.30'
#> 95% of iterations needed (% iterations that need to be successful in order to consider the community successful)
percN <- '0.95'
# percN <- '0.90'
# no of dil-transf cycles (for the y axis limits)
n_cycles <- 200

# OTHER OPTIONS ----------------------------------------------------------------
threshold <- as.numeric(threshold); percN <- as.character(percN) #> numeric
prefix <- paste0(threshold*100, "_")

# FILE OPTIONS ----------------------------------------------------------------
#> different output folder ("_groups")
out_folder = "../figures_groups_INTERACTIONS/"
out_folder <- paste0(out_folder, "/p", percN, "_f", threshold, "/")
if (!file.exists(out_folder)) {system(paste("mkdir -p", out_folder))}
# ---
# no interactions
NO_INTER_FILE <- paste0("group_PREPARSED_analysis_table_p", percN, "_f", threshold, ".csv")
INTER_FILE <- "SkewedGroups_N10_100sp"

# READ FILES -------------------------------------------------------------------
csv_groups <- read_csv(NO_INTER_FILE)

# with interactions
folder <- paste0("../../2024-10-15_deriva_3/results_01_", INTER_FILE, "/stats/")
csv_interactions <- tibble()
for (f in list.files(folder, full.names = F)) {
  read_f <- read_csv(paste0(folder, f)) %>% 
    mutate(nicheN = str_split_i(INTER_FILE, "_", 2) %>% str_split_i(., "N", 2),
           nichedist = str_split_i(INTER_FILE, "_", 1),
           richness = str_split_i(INTER_FILE, "_", 3) %>% str_split_i(., "sp", 1))
  read_f$group_numbers <- paste(1:unique(read_f$nicheN), collapse = ";") # IMPORTANTE para la posterior identificación # TODO es diferente para csv_groups pero no da problemas
  csv_interactions <- csv_interactions %>% bind_rows(read_f)
}

# PARSE : no inter -------------------------------------------------------------
df_long <- csv_groups %>%
  mutate(group_success = strsplit(group_success, ";")) %>%  # Separar los valores
  unnest(group_success) %>%  # Expandir filas
  mutate(group_success = 100 * as.numeric(group_success)) %>%  # Convertir a numérico
  group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
  mutate(group_number = row_number()) %>%  # Asignar el número de grupo
  ungroup()

df_summary <- df_long %>% # Calcular la media de los valores acumulativos para cada número de grupos
  group_by(nicheN, group_number, dilfactor, nichedist, richness) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
  summarise(mean_success = mean(group_success), .groups = "drop") %>%
  mutate(above_percN = mean_success > (100 * as.numeric(percN)),  # TRUE si está por encima, FALSE si no
         inter_type = NA) #> to separate this dataset from the interactions data

# PARSE : inter ----------------------------------------------------------------
df_long_i <- csv_interactions %>%
  mutate(group_number = strsplit(group_numbers, ";")) %>%  # Separar los valores
  unnest(group_number) %>%  # Expandir filas
  mutate(group_number = as.numeric(group_number)) %>% 
  mutate(group_success_rate = map2_chr(group_success_rate, group_number, ~ str_split_i(.x, ";", .y))) %>% # Seleccionar el índice correspondiente
  mutate(group_success_rate = 100 * as.numeric(group_success_rate)) %>%  # Convertir a numérico
  # group_by(nicheN, dilfactor, nichedist, richness, Interaction_presence, Interaction_value, Interaction_sign, Mutual_interactions,
           # shannons, group_extinctions, Fixation_threshold, `Number of iterations`, `Number of dilution cycles`, success_rate, group_success_rate) %>%  # Mantener todas estas variables para luego
  mutate(inter_type = paste(Interaction_presence, Interaction_value, Interaction_sign, Mutual_interactions, sep = "_")) %>% 
  ungroup()

df_summary_i <- df_long_i %>% # Calcular la media de los valores acumulativos para cada número de grupos
  group_by(nicheN, group_number, dilfactor, nichedist, richness,
           inter_type, # importante
           Interaction_presence, Interaction_value, Interaction_sign, Mutual_interactions, # interaction specifics...
           ) %>%  # se conservan para el plot luego
  summarise(mean_success = mean(group_success_rate), .groups = "drop") %>%
  mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no

# PLOT - fixation --------------------------------------------------------------
# we're going to plot the frequency of every nicheN-Percentage of Groups with Success combination.
# (!) "group_extinctions" es en csv_groups (comprobar aún así) un string con ";" separando ! 
# porcentaje de los grupos que sí tienen éxito
for (case in c(1, 2)) {
  if (case == 1) 
  {
    nichedistvalue <- "SkewedGroups"
    dilfactorvalue <- 0.05
    richnessvalue <- 100
    nicheNvalue <- 10
  }
  else if (case == 2)
  {
    nichedistvalue <- "SkewedGroups"
    dilfactorvalue <- 0.1
    richnessvalue <- 100
    nicheNvalue <- 10
  }
  # no inter
  df_summary_filtered <- df_summary[((df_summary$nichedist == nichedistvalue) &
                                       (df_summary$dilfactor == dilfactorvalue) &
                                       (df_summary$nicheN == nicheNvalue) &
                                       (df_summary$richness == richnessvalue)), ] %>%
    mutate(nicheN = as.factor(nicheN))
  
  # inter
  df_summary_filtered_i <- df_summary_i[((df_summary_i$nichedist == nichedistvalue) &
                                           (df_summary_i$dilfactor == dilfactorvalue) &
                                           (df_summary_i$nicheN == nicheNvalue) &
                                           (df_summary_i$richness == richnessvalue)), ] %>%
    mutate(nicheN = as.factor(nicheN))
  
  # plots[[nichedistvalue]] <- ggplot(df_summary_filtered, aes(x = group_number, y = mean_success, color = as.factor(nicheN))) +
  p <- ggplot(df_summary_filtered_i, aes(x = group_number, y = mean_success, color = nicheN, fill = nicheN)) +
    geom_point(size = 1) +
    labs(x = "Número de grupo",
         y = "Tasa de extinción (%) / abundancia relativa",
         color = "Número de grupos (nicheN)") +
    theme_minimal() +
    # facet_grid(. ~ Interaction_presence + Interaction_value + Interaction_sign + Mutual_interactions) +
    # facet_wrap(
    #            . ~ inter_type,
    #            axis.labels = "margins",
    #            strip.position = "top",
    #            scales = "free_x",
    #            ) +
    facet_grid(Interaction_presence + Interaction_value ~ Interaction_sign + Mutual_interactions) +
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

  p
  
  ggsave(
    filename = paste0(out_folder,"/0g_extinction_per_group_", nichedistvalue, "_", richnessvalue, "_", nicheNvalue, ".png"),
    plot = p,
    height = 1.2 * 2000, width = 2 * 1100,
    units = "px"
  )
}
