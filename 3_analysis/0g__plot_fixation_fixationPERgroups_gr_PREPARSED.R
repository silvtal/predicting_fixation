# 'PREPARSED' : datos mejor separados // better separated data
# generated with ~/repos/predicting_fixation/3_analysis/g_plot_results_facet_by_dil.r

# THERE'S ANOTHER SCRIPT FOR 'PARSED', DATA GENERATED FROM 'generate_data(...)'

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

#> different input folder(s!)
nogroups_file <- paste0("../1_datasets/simulation_results/processed_data_simcomms_", threshold, "_full_jul") # this one's 'PARSED', always.
# groups_file   <- paste0("../1_datasets/simulation_results/2024_parsed_results_GROUPS_p", percN, "_f", threshold, "_all.csv")
groups_file   <- paste0("../1_datasets/simulation_results/preparsed/all_information.csv") # calculado para threshold = 0.50

threshold <- as.numeric(threshold); percN <- as.character(percN) #> numeric
prefix <- paste0(threshold*100, "_")

plot_bars <- TRUE #> whether to plot niche size bars at the background

# LOAD DATA --------------------------------------------------------------------
#> Load the N > 1 data
csv_groups <- read.csv(groups_file, fileEncoding = "UTF-8", sep = ",")
print(paste("Muestras con grupos procesadas:", nrow(csv_groups)))

#> FIRST, adapt the groups csv to another format, compatible with older scripts
colnames(csv_groups) <- c("sample", # FROM: community    # TODO debería ser 1_lognorm_1e+06_1000
                          "group_success", # FROM: "group_success_rate", # nos ahorra calcularlo luego
                          "success_rate", # same
                          "group_extinctions", # same
                          "num_groups", # nos ahorra sacarlo desde el nombre
                          "richness",
                          "nichedist", # FROM: "group_distribution",
                          "dilfactor", # FROM: "dilution_factor",
                          "distrib", # FROM: distribution",
                          "community_type", # ej. "10_Even_100"
                          "comm_number", # ej. "17"
                          "niter", # FROM: "length",
                          "gini", # FROM: "mean_gini",
                          "evenness", # FROM: "mean_evenness",
                          "shannon", # FROM: "mean_shannon",
                          "ginis",
                          "evennesses",
                          "shannons")
# Add some fields
csv_groups$nichedist <- paste0(csv_groups$nichedist, "Groups")
csv_groups$nicheN <- lapply(csv_groups$community_type, FUN = function(ct) stringr::str_split_1(ct, "_")[[1]][1]) %>% as.numeric
csv_groups["logdf"] <- log(as.numeric(csv_groups$dilfactor))

# FIX -- TODO -- esto es para mi set de datos concreto !!
# For additional 
csv_groups[csv_groups$community_type == "3_Even_100", "abundances"] <- paste(c(0.33, 0.33, 0.34), collapse = ";")
csv_groups[csv_groups$community_type == "3_Even_1000", "abundances"] <- paste(c(0.33, 0.33, 0.34), collapse = ";")
csv_groups[csv_groups$community_type == "3_Skewed_100", "abundances"] <- paste(c(0.6, 0.3, 0.1), collapse = ";")
csv_groups[csv_groups$community_type == "3_Skewed_1000", "abundances"] <- paste(c(0.6, 0.3, 0.1), collapse = ";")
csv_groups[csv_groups$community_type == "10_Even_100", "abundances"] <- paste(rep(0.1, 10), collapse = ";")
csv_groups[csv_groups$community_type == "10_Even_1000", "abundances"] <- paste(rep(0.1, 10), collapse = ";")
csv_groups[csv_groups$community_type == "10_Skewed_100", "abundances"] <- paste(c(0.3, 0.2, 0.15, 0.12, 0.09, 0.07, 0.04, 0.02, 0.009, 0.001), collapse = ";")
csv_groups[csv_groups$community_type == "10_Skewed_1000", "abundances"] <- paste(c(0.3, 0.2, 0.15, 0.12, 0.09, 0.07, 0.04, 0.02, 0.009, 0.001), collapse = ";")

# Save the final processed table for any subsequent analyses
write_csv(csv_groups, paste0("group_PREPARSED_analysis_table_p", percN, "_f", threshold, ".csv"))

# If our table includes the cycle in which success happens
if ("success" %in% colnames(csv_groups)) {
  #> When no success, encode that as 0 instead of NA
  csv_groups$success       <- csv_groups$success %>% replace(is.na(.), 0)
  csv_groups$group_success <- csv_groups$group_success %>% str_replace_all("NA", "0")
  csv_groups$group_success <- str_split(csv_groups$group_success, ";")
}

# PLOT - success ---------------------------------------------------------------
#> Now load the single-group data, if it exists, and create the nicheN plots only in this case
#> "SUCCESS" NEEDS TO BE IN THE GROUPS CSV
if (file.exists(nogroups_file) && ("success" %in% colnames(csv_groups))) {
  csv_nogroups <- read.csv(nogroups_file, fileEncoding = "UTF-8", sep = ",")
  csv_nogroups["logdf"] <- log(csv_nogroups$dilfactor)
  
  #> When no success, encode that as 0 instead of NA
  csv_nogroups$success <- csv_nogroups$success %>% replace(is.na(.), 0)
  
  #> IMPORTANT to remove data that we have no equivalent for in the multiple-group
  #> results
  csv_nogroups <- csv_nogroups[csv_nogroups$richness %in% unique(csv_groups$richness) & 
                                 csv_nogroups$size == 10000 & 
                                 csv_nogroups$distrib == "lognorm", ]  # lognorm only !!!
  
  #> Add fields absent in the single-group dataset
  csv_nogroups$group_success <- csv_nogroups$success %>% as.character() %>% as.list
  csv_nogroups$nicheN <- 1
  csv_nogroups$nichedist <- "NoGroups"
  
  csv_g_and_ng <- bind_rows(csv_groups, csv_nogroups)
  
  #> EvenGroups; success (functional group no. vs. dilfactor effect)---
  png(paste0(out_folder, "/0.2_success_v_nicheN_Even.png"),
      width = 800, height = 600)
  # f <- flexplot(success~dilfactor + nicheN, data=csv_groups)
  f <- flexplot(success~dilfactor + nicheN,
                data=csv_g_and_ng[csv_g_and_ng$nichedist %in% c("EvenGroups", "NoGroups"), ])
  plot(f)
  dev.off()
  
  #> SkewedGroups; success (functional group no. vs. dilfactor effect)---
  png(paste0(out_folder, "/0.2_success_v_nicheN_Skewed.png"),
      width = 800, height = 600)
  # f <- flexplot(success~dilfactor + nicheN, data=csv_groups)
  f <- flexplot(success~dilfactor + nicheN,
                data=csv_g_and_ng[csv_g_and_ng$nichedist %in% c("SkewedGroups", "NoGroups"), ])
  plot(f)
  dev.off()
}

# PLOT - fixation --------------------------------------------------------------
#> EvenGroups; GROUP fixation (functional group no. vs. dilfactor effect)---
# we're going to plot the frequency of every nicheN-Percentage of Groups with Success combination.
# (!) "group_extinctions" es un string con ";" como separador
# porcentaje de los grupos que sí tienen éxito

if (plot_bars) {
  df_long <- csv_groups %>%
    mutate(group_success = strsplit(group_success, ";")) %>%  # Separar los valores
    unnest(group_success) %>%  # Expandir filas
    mutate(group_success = 100 * as.numeric(group_success)) %>%  # Convertir a numérico
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Mantener todas estas variables para luego
    mutate(group_number = row_number()) %>%  # Asignar el número de grupo
    mutate(abundances = 100 * str_split(abundances, ";")[[1]][group_number] %>% as.numeric) %>% 
    ungroup()
  
  df_summary <- df_long %>% # Calcular la media de los valores acumulativos para cada número de grupos
    group_by(nicheN, group_number, dilfactor, nichedist, richness,
             abundances) %>%  # dilfactor, richness y nichedist se conservan para el plot luego
    summarise(mean_success = mean(group_success), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
  
} else {
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
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE si está por encima, FALSE si no
}

# plot(s)
for (nichedistvalue in c("EvenGroups", "SkewedGroups")){
  for (richnessvalue in c(100, 1000)){
    df_summary_filtered <- df_summary[((df_summary$nichedist == nichedistvalue) &
                                        (df_summary$richness == richnessvalue)), ]
    
    p <- ggplot(df_summary_filtered,
                aes(x = group_number, y = mean_success,
                    color = interaction(as.factor(nicheN), above_percN), fill = interaction(as.factor(nicheN), above_percN))) +
      geom_point(size = 1) +
      geom_hline(yintercept = 100 * as.numeric(percN), color = "coral", linewidth = .41, linetype = "dashed") +  # Línea roja horizontal
      labs(x = "Número de grupo", 
           y = "Tasa de fijación (%) / Tamaño de nicho (%)",
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
      # scale_y_continuous(
      #   sec.axis = sec_axis(~., name = "Tamaño de nicho (%)")
      # ) +
      theme(
        legend.position = "none",
        panel.grid.major = element_line(color = "lightgray", size = 0.25),
        panel.grid.minor = element_line(color = "gray90", size = 0.1)
      ) +
      scale_color_manual(values = c("#FB9A99", "#B2DF8A", "#E31A1C", "#33A02C")) +
      scale_fill_manual(values = c("#FB9A99", "#B2DF8A", "#E31A1C", "#33A02C")) +
      scale_x_continuous(breaks = function(x) seq(0, max(x), by = 1))  # Solo números naturales hasta nicheN
  
    if (plot_bars) {
      p <- p + 
        geom_bar(aes(x = group_number, y = abundances),
                 stat = "identity", position = "stack", alpha = 0.2)
    }
    
    ggsave(
      filename = paste0(out_folder,"/0g_success_per_group_", nichedistvalue, "_", richnessvalue, ".png"),
      plot = p,
      height = 1.2 * 2000, width = 2 * 1100,
      units = "px"
    )
  }
}
