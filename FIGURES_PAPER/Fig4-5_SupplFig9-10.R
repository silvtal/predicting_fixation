# CURSOR MADE
# Script that combines both extinction and fixation data
# Based on 0g__plot_extinction_extinctionPERgroups_gr_PREPARSED.R and 0g__plot_fixation_fixationPERgroups_gr_PREPARSED.R

library(tidyverse)
library(flexplot)
#> run from script's directory !!

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
out_folder <- "."
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
    mutate(group_extinctions = strsplit(group_extinctions, ";")) %>%  # Separate values
    unnest(group_extinctions) %>%  # Expand rows
    mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convert to numeric
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Keep all these variables for later
    mutate(group_number = row_number()) %>%  # Assign group number
    mutate(abundances = 100 * str_split(abundances, ";")[[1]][group_number] %>% as.numeric) %>% 
    ungroup()
  
  df_summary_extinction <- df_long_extinction %>% # Calculate mean for each group number
    group_by(nicheN, group_number, dilfactor, nichedist, richness,
             abundances) %>%  # dilfactor, richness, and nichedist kept for plot
    summarise(mean_success = mean(group_extinctions), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE if above, FALSE if not
  
} else {
  df_long_extinction <- csv_groups_extinction %>%
    mutate(group_extinctions = strsplit(group_extinctions, ";")) %>%  # Separate values
    unnest(group_extinctions) %>%  # Expand rows
    mutate(group_extinctions = 100 * as.numeric(group_extinctions)) %>%  # Convert to numeric
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Keep all these variables for later
    mutate(cumulative_count = row_number()) %>%  # Assign cumulative number within each sample
    ungroup()
  
  df_summary_extinction <- df_long_extinction %>% # Calculate mean for each group number
    group_by(nicheN, cumulative_count, dilfactor, nichedist, richness) %>%  # dilfactor, richness, and nichedist kept for plot
    summarise(mean_success = mean(group_extinctions), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE if above, FALSE if not
}

# Process fixation data
if (plot_bars) {
  df_long_fixation <- csv_groups_fixation %>%
    mutate(group_success = strsplit(group_success, ";")) %>%  # Separate values
    unnest(group_success) %>%  # Expand rows
    mutate(group_success = 100 * as.numeric(group_success)) %>%  # Convert to numeric
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Keep all these variables for later
    mutate(group_number = row_number()) %>%  # Assign group number
    mutate(abundances = 100 * str_split(abundances, ";")[[1]][group_number] %>% as.numeric) %>% 
    ungroup()
  
  df_summary_fixation <- df_long_fixation %>% # Calculate mean for each group number
    group_by(nicheN, group_number, dilfactor, nichedist, richness,
             abundances) %>%  # dilfactor, richness, and nichedist kept for plot
    summarise(mean_success = mean(group_success), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE if above, FALSE if not
  
} else {
  df_long_fixation <- csv_groups_fixation %>%
    mutate(group_success = strsplit(group_success, ";")) %>%  # Separate values
    unnest(group_success) %>%  # Expand rows
    mutate(group_success = 100 * as.numeric(group_success)) %>%  # Convert to numeric
    group_by(sample, nicheN, dilfactor, nichedist, richness) %>%  # Keep all these variables for later
    mutate(group_number = row_number()) %>%  # Assign group number
    ungroup()
  
  df_summary_fixation <- df_long_fixation %>% # Calculate mean for each group number
    group_by(nicheN, group_number, dilfactor, nichedist, richness) %>%  # dilfactor, richness, and nichedist kept for plot
    summarise(mean_success = mean(group_success), .groups = "drop") %>%
    mutate(above_percN = mean_success > (100 * as.numeric(percN)))  # TRUE if above, FALSE if not
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
    
    p <- ggplot()

      if (plot_bars) { 
      p <- p + geom_bar(data = df_summary_filtered_fixation,
                      aes(x = group_number, y = abundances, fill = interaction(nicheN, above_percN), color = after_scale(fill)),
                      stat = "identity", position = "stack", alpha = 0.2)
      }
      # Plot extinction data in darker gray
      p <- p +
      geom_point(data = df_summary_filtered_extinction,
                aes(x = group_number, y = mean_success),
                color = "#333333", size = 1) +
      # Plot fixation data in color, fully opaque, color depends on above_percN
      geom_point(data = df_summary_filtered_fixation, 
                aes(x = group_number, y = mean_success, color = interaction(nicheN, above_percN), fill = interaction(nicheN, above_percN)),
                size = 1, shape = 21, alpha = 1) +
      geom_hline(yintercept = 100 * as.numeric(percN), color = "coral", linewidth = .41, linetype = "dashed") +

      labs(x = "Group number",
           y = "Extinction or fixation rate (%) / Niche size (%)",
           color = "Number of groups (nicheN)") +
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
      scale_color_manual(values = c("#FB9A99", "#B2DF8A", "#E31A1C", "#33A02C")) +
      scale_fill_manual(values = c("#FB9A99", "#B2DF8A", "#E31A1C", "#33A02C")) +
      scale_x_continuous(breaks = function(x) seq(0, max(x), by = 1))  # Only natural numbers up to nicheN

    # Update output filenames to new numbering
    fig_prefix <- ifelse(nichedistvalue == "EvenGroups" && richnessvalue == 100, "Fig6",
                   ifelse(nichedistvalue == "EvenGroups" && richnessvalue == 1000, "SupplFig7",
                   ifelse(nichedistvalue == "SkewedGroups" && richnessvalue == 100, "Fig7",
                   ifelse(nichedistvalue == "SkewedGroups" && richnessvalue == 1000, "SupplFig8", ""))))
    ggsave(
      filename = paste0(out_folder, "/", fig_prefix, "_0g_extinction_and_fixation_per_group_", nichedistvalue, "_", richnessvalue, ".png"),
      plot = p,
      height = 1.2 * 2000, width = 2 * 1100,
      units = "px"
    )
  }
} 