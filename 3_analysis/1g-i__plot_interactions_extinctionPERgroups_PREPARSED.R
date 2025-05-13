# Comparamos la tasa de extinción por grupo con y sin interacciones de varios tipos.
# 'PREPARSED' : datos mejor separados y tal

library(tidyverse)
library(flexplot)
#> run from script's directory !!
setwd("~/repos/predicting_fixation/3_analysis/")

# FIXED OPTIONS ----------------------------------------------------------------
# (solo sirven para el output_folder, quitar)
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
IS_THE_ORDER_CORRECT_inter <- FALSE # set to FALSE if : 1, 10, 2, 3... in shannons, group_success and group_extinctions

# FILE OPTIONS ----------------------------------------------------------------
#> different output folder ("_groups")
out_folder = "../figures_groups_INTERACTIONS/"
out_folder <- paste0(out_folder, "/p", percN, "_f", threshold, "/")
if (!file.exists(out_folder)) {system(paste("mkdir -p", out_folder))}
# ---
# no interactions
NO_INTER_FILE <- paste0("group_PREPARSED_analysis_table_p", percN, "_f", threshold, ".csv")
INTER_FOLDERS <- c("SkewedGroups_N10_100sp", "EvenGroups_N3_100sp")
                 

# READ FILES -------------------------------------------------------------------
csv_groups <- read_csv(NO_INTER_FILE)

# with interactions
csv_interactions <- tibble()
for (INTER_FOLDER in INTER_FOLDERS) {
  folder <- paste0("../../2024-10-15_deriva_3/results_01_", INTER_FOLDER, "/stats/")
  for (f in list.files(folder, full.names = F)) {
    read_f <- read_csv(paste0(folder, f)) %>% 
      mutate(nicheN = str_split_i(INTER_FOLDER, "_", 2) %>% str_split_i(., "N", 2),
             nichedist = str_split_i(INTER_FOLDER, "_", 1),
             richness = str_split_i(INTER_FOLDER, "_", 3) %>% str_split_i(., "sp", 1))
    read_f$group_numbers <- paste(1:unique(read_f$nicheN), collapse = ";") # IMPORTANTE para la posterior identificación # TODO es diferente para csv_groups pero no da problemas
    csv_interactions <- csv_interactions %>% bind_rows(read_f)
  }
}

# PARSE : no inter -------------------------------------------------------------
df_long <- csv_groups %>%
  mutate(group_extinctions = strsplit(group_extinctions, ";")) %>%  # Separar los valores
  unnest(group_extinctions) %>%  # Expandir filas
  mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convertir a numérico
  group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
  mutate(group_number = row_number()) %>%  # Asignar el número de grupo
  ungroup()

df_summary <- df_long %>% # Calcular la media de los valores acumulativos para cada número de grupos
  group_by(nicheN, group_number, dilfactor, nichedist, richness) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
  summarise(mean_extinction = mean(group_extinctions), .groups = "drop") %>%
  mutate(above_percN = mean_extinction > (100 * as.numeric(percN)),  # TRUE si está por encima, FALSE si no
         inter_type = NA) #> to separate this dataset from the interactions data

# PARSE : inter ----------------------------------------------------------------
df_long_i <- csv_interactions %>%
  mutate(group_number = strsplit(group_numbers, ";")) %>%  # Separar los valores
  unnest(group_number) %>%  # Expandir filas
  mutate(group_number = as.numeric(group_number)) %>% 
  mutate(group_extinctions = map2_chr(group_extinctions, group_number, ~ str_split_i(.x, ";", .y))) %>% # Seleccionar el índice correspondiente
  mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convertir a numérico
  mutate(inter_type = paste(Interaction_presence, Interaction_value, Interaction_sign, Mutual_interactions, sep = "_")) %>% 
  ungroup()

df_summary_i <- df_long_i %>% # Calcular la media de los valores acumulativos para cada número de grupos
  group_by(nicheN, group_number, dilfactor, nichedist, richness,
           inter_type, # importante
           Interaction_presence, Interaction_value, Interaction_sign, Mutual_interactions, # interaction specifics...
           ) %>%  # se conservan para el plot luego
  summarise(mean_extinction = mean(group_extinctions), .groups = "drop") %>%
  mutate(above_percN = mean_extinction > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no

if (!IS_THE_ORDER_CORRECT_inter) {
  df_summary_i <- df_summary_i %>%
    mutate(group_number = ifelse(nicheN == 10 & group_number == 9, 2,
                                 ifelse(nicheN == 10 & group_number == 2, 9, group_number)))
  IS_THE_ORDER_CORRECT_inter <- TRUE
}

# PLOT - fixation --------------------------------------------------------------
# Cada plot recoge una combinación de valores de 4 variables
for (INTER_FOLDER in INTER_FOLDERS) {
  {
    nicheNvalue = str_split_i(INTER_FOLDER, "_", 2) %>% str_split_i(., "N", 2)
    nichedistvalue = str_split_i(INTER_FOLDER, "_", 1)
    richnessvalue = str_split_i(INTER_FOLDER, "_", 3) %>% str_split_i(., "sp", 1)
    
    if (nicheNvalue == 3) {
      case_scale <- RColorBrewer::brewer.pal(9, "Reds")[(9 - length(unique(df_summary_i$Interaction_value))):9]
    } else {
      case_scale <- RColorBrewer::brewer.pal(9, "Greens")[(9 - length(unique(df_summary_i$Interaction_value))):9]
    }
  }
  
  for (dilfactorvalue in c(0.05, 0.1)) {
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
      mutate(nicheN = as.factor(nicheN)) %>%
      mutate(group_number = as.numeric(group_number)) 
    
    # plots[[nichedistvalue]] <- ggplot(df_summary_filtered, aes(x = group_number, y = mean_extinction, color = as.factor(nicheN))) +
    p <- ggplot(df_summary_filtered_i, aes(x = group_number, y = mean_extinction, color = as.factor(Interaction_value), group = Interaction_value)) +  # next
      geom_line(linewidth = .3) +
      geom_point(size = 1) +
      labs(title = paste(INTER_FOLDER, dilfactorvalue),
           x = "Número de grupo",
           y = "Tasa de extinción (%) / abundancia relativa",
           color = "Valor de la interacción") +
      theme_minimal() +
      facet_grid(
        rows = vars(Interaction_presence),
        cols = vars(Interaction_sign),
        axis.labels = "all",
        labeller = label_both # otros: label:both
        # facet_wrap(
        # Interaction_presence + Interaction_sign ~ ., # TODO si hubiera Mutual, se pondría aquí.
      ) +
      coord_cartesian(ylim = c(0, 100)) +
      # theme(legend.position = "bottom") +
      theme(
        legend.position = "bottom",
        panel.grid.major = element_line(color = "lightgray", size = 0.25),  # Major grid lines
        panel.grid.minor = element_line(color = "gray90", size = 0.1),      # Minor grid lines
        plot.background = element_rect(fill = "gray90")
      ) +
      scale_color_manual(values = case_scale) +
      scale_x_continuous(breaks = function(x) seq(0, max(x), by = 1)) + # Solo números naturales hasta nicheN
      # añado sin inter:  
      geom_point(data = df_summary_filtered, aes(x = group_number, y = mean_extinction), size = 1, inherit.aes = F) +
      geom_line(data = df_summary_filtered, aes(x = group_number, y = mean_extinction), size = .5, inherit.aes = F)
    
    ggsave(
      filename = paste0(out_folder,"/0g_extinction_per_group_", INTER_FOLDER, "_", dilfactorvalue, ".png"),
      plot = p,
      height = 1.2 * 2000, width = 2 * 1100,
      units = "px"
    )
  }
}

