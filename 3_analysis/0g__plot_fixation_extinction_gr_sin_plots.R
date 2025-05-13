# Notas
# - El success más bajo posible para single-group y para multiple-group data es
#   2. Porque suma un ciclo por la cara. O sea que resto 1 a todos los valores.
#   Esto se hace en la parte de análisis, no al cargar los datos.
# - The multiple-group simulations all fail with a dilution factor of 0.25, so
#   we only have data for 8 different values instead of 9.

library(tidyverse)
library(flexplot)
#> run from script's directory

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
groups_file   <- paste0("../1_datasets/simulation_results/2024_parsed_results_GROUPS_p", percN, "_f", threshold, "_all.csv")
nogroups_file <- paste0("../1_datasets/simulation_results/parsed_results_1group_p", percN, "_f", threshold, "_all.csv")
threshold <- as.numeric(threshold); percN <- as.character(percN) #> numeric"
prefix <- paste0(threshold*100, "_")

# LOAD DATA --------------------------------------------------------------------
#> Load the 3 and 10 groups data
csv_groups <- read.csv(groups_file, fileEncoding = "UTF-8", sep = ",")
csv_groups["logdf"] <- log(csv_groups$dilfactor)

#> for filtering
my_dilfactors <- csv_groups$dilfactor %>% unique #> interested in all the common dilfactors (no filter)
# only 1-group --> 0.00025 0.00040 0.00400 0.00800 0.04000 0.25000 0.50000 0.80000
# 3- and 10-group samples have (global) success with the same 8 values.

#> When no success, encode that as 0 instead of NA
csv_groups$success       <- csv_groups$success %>% replace(is.na(.), 0)
csv_groups$group_success <- csv_groups$group_success %>% str_replace_all("NA", "0")
csv_groups$group_success <- str_split(csv_groups$group_success, ";")

#> Now load the single-group data, if it exists
if (file.exists(nogroups_file)) {
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
  
  csv <- bind_rows(csv_groups, csv_nogroups)
} else {
  #> Proceed with only the groups data if no groups file is absent
  csv <- csv_groups
}

print(paste("Muestras procesadas:", nrow(csv)))

# for debugging:
# csv <- csv %>%
#   na.omit()
# print(paste("Eliminando los NA:", nrow(csv)))

csv <- csv[csv$dilfactor %in% my_dilfactors, ]
print(paste("Filtrando por dilution factor ", nrow(csv)))

# !! The number of available replicates could be different from sample to sample
# We add a column with the number of iterations for each sample.
csv$niter <- csv$group_sizes %>% stringr::str_count(";") + 1

# Save the final processed table for any subsequent analyses
write_csv(csv, paste0("group_analysis_table_p", percN, "_f", threshold, ".csv"))

# PLOT - success ---------------------------------------------------------------
#> EvenGroups; success (functional group no. vs. dilfactor effect)---

# png(paste0(out_folder, "/0.2_success_v_nicheN.png"),
#     width = 800, height = 600)
# f <- flexplot(success~dilfactor + nicheN, data=csv)
f <- flexplot(success~dilfactor + nicheN,
              data=csv[csv$nichedist %in% c("EvenGroups", "NoGroups"), ])
plot(f)
# dev.off()

# PLOT - fixation --------------------------------------------------------------
#> EvenGroups; GROUP fixation (functional group no. vs. dilfactor effect)---
# we're going to plot the frequency of every nicheN-Percentage of Groups with Success combination.
csv$perc_groups_success <- 100 * sapply(csv$group_success, function(x) sum(x!="0")) / csv$nicheN
csv <- csv %>%
  group_by(nicheN, perc_groups_success) %>%
  mutate(total_samples_with_partial_success = n())

# Plot 1/2
csv_selection <- csv[
  # select Even and no groups
  (csv$nichedist %in% c("EvenGroups", "NoGroups")) &
    # don't show 0s (only count samples with SOME partial success)
    (csv$perc_groups_success > 0), ]

if (nrow(csv_selection) > 0) {
  # png(paste0(out_folder, "/0.2_group_success_v_nicheN_v_dilfactor_EVEN.png"),
  # width = 800, height = 600)
  
  p <- ggplot(csv_selection,
              aes(x = nicheN, y = total_samples_with_partial_success, color = perc_groups_success, shape = as.factor(richness))) + 
    geom_point(size = 4) +
    labs(title = "Percentage of fixated groups (Even Groups)",
         x = "Niche N",
         y = "Occurrences / Samples",
         color = "Percentage of Groups with Success (%)",
         shape = "Richness") + scale_shape_manual(values = c(1, 4)) +
    facet_wrap(~ dilfactor) +
    theme(legend.position = "bottom") +
    scale_color_gradient(low = "lightgreen", high = "darkred", limits = c(0, 100)) +
    scale_x_continuous(breaks = unique(csv$nicheN)) +
    ylim(0, n_cycles)
  
  plot(p)
  
  # dev.off()
} else {
  message(paste0("no samples with Even groups with partial success (threshold: ", threshold, ", percN: ", percN, ")"))
}


# Plot 2/2
#> SkewedGroups; GROUP fixation (functional group no. vs. dilfactor effect)---
csv_selection <- csv[
  # select Uneven
  (csv$nichedist %in% c("SkewedGroups", "NoGroups")) &
    # don't show 0s (only count samples with SOME partial success)
    (csv$perc_groups_success > 0), ]

if (nrow(csv_selection) > 0) {
  png(paste0(out_folder, "/0.2_group_success_v_nicheN_v_dilfactor_SKEWED.png"),
      width = 800, height = 600)
  
  p <- ggplot(csv_selection,
              aes(x = nicheN, y = total_samples_with_partial_success, color = perc_groups_success, shape = as.factor(richness))) +
    geom_point(size = 4) +
    labs(title = "Percentage of fixated groups (Uneven Groups)",
         x = "Niche N",
         y = "Occurrences / Samples",
         color = "Percentage of Groups with Success (%)",
         shape = "Richness") + scale_shape_manual(values = c(1, 4)) +
    facet_wrap(~ dilfactor) +
    theme(legend.position = "bottom") +
    scale_color_gradient(low = "lightgreen", high = "darkred", limits = c(0, 100)) +
    scale_x_continuous(breaks = unique(csv$nicheN)) +
    ylim(0, n_cycles)
  
  plot(p)
  
  dev.off()
} else {
  message(paste0("no samples with uneven groups with partial success (threshold: ", threshold, ", percN: ", percN, ")"))
}

# START ------------------------------------------------------------------------
#> Check for extinction (Even; N 3 & 10)
#> $total_iters_with_extinction --> number of iterations where at least one group is lost (no matter how many)
#> $mean_extinct_groups_per_iter_with_extinction --> mean number of groups lost in these iterations
for (i in 1:nrow(csv)) {
  x <- csv[i, ]$group_sizes
  
  extinct_groups_per_iter <-
    x %>% str_split(";") %>% .[[1]] %>% str_remove_all(" ") %>% str_split(",") %>% 
    sapply(., function(vec) sum(vec == "0"))
  
  csv[i, "total_iters_with_extinction"] <-
    sum(extinct_groups_per_iter!=0) #> !=0 means there's at least 1 extinct group
  
  csv[i, "mean_extinct_groups_per_iter_with_extinction"] <-
    sum(extinct_groups_per_iter) / csv[i, "total_iters_with_extinction"]
}

#> Plot 1/2
# png(paste0(out_folder, "/0.2_extinctions_v_nicheN_v_dilfactor_EVEN.png"),
# width = 800, height = 600)

ggplot(csv[csv$nichedist %in% c("EvenGroups") &
             # don't show 0s (only count samples with SOME extinction)
             (csv$total_iters_with_extinction > 0), ],
       aes(x = nicheN,
           y = 100 * total_iters_with_extinction / niter, # important to divide by niter!
           color = mean_extinct_groups_per_iter_with_extinction,
           shape = as.factor(richness))) +
  geom_point(size = 4) +
  labs(title = "Extinction rate vs niche number vs dilution factor (Even Groups)",
       x = "Dilution factor",
       y = "Iterations with extinct groups (%)",
       color = "Mean # of extinct groups",
       shape = "Richness") + scale_shape_manual(values = c(1, 4)) +
  facet_wrap(~ dilfactor) +
  theme(legend.position = "bottom") +
  scale_color_gradient(low = "lightgreen", high = "darkred") +
  scale_x_continuous(breaks = unique(csv$nicheN))

# dev.off()

# Plot 2/2
# png(paste0(out_folder, "/0.2_extinctions_v_nicheN_v_dilfactor_SKEWED.png"),
# width = 800, height = 600)

ggplot(csv[csv$nichedist %in% c("SkewedGroups") &
             # don't show 0s (only count samples with SOME extinction)
             (csv$total_iters_with_extinction > 0), ],
       aes(x = nicheN,
           y = 100 * total_iters_with_extinction / niter, # important to divide by niter!
           color = mean_extinct_groups_per_iter_with_extinction,
           shape = as.factor(richness))) +
  geom_point(size = 4) +
  labs(title = "Extinction rate vs niche number vs dilution factor (Uneven Groups)",
       x = "Dilution factor",
       y = "Iterations with extinct groups (%)",
       color = "Mean # of extinct groups",
       shape = "Richness") + scale_shape_manual(values = c(1, 4)) +
  facet_wrap(~ dilfactor) +
  theme(legend.position = "bottom") +
  scale_color_gradient(low = "lightgreen", high = "darkred") +
  scale_x_continuous(breaks = unique(csv$nicheN))

# dev.off()