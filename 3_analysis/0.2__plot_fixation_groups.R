# Notas
# - El success más bajo posible para single-group y para multiple-group data es
#   2. Porque suma un ciclo por la cara. O sea que resto 1 a todos los valores.
#   Esto se hace en la parte de análisis, no al cargar los datos.
# - The multiple-group simulations all fail with a dilution factor of 0.25, so
#   we only have data for 8 different values instead of 9.

library(tidyverse)
library(flexplot)

#> run from script's directory
#> different results folder ("_groups")
out_folder = "../figures_groups/"; if (!file.exists(out_folder)) {system(paste("mkdir -p", out_folder))}

# OPTIONS ----------------------------------------------------------------------
#> only 50% fixation threshold
threshold <- 0.5
#> different results folder(s!)
groups_file   <- paste0("../1_datasets/simulation_results/2024_parsed_results_GROUPS_", threshold, "_all.csv")
nogroups_file <- paste0("../1_datasets/simulation_results/processed_data_simcomms_", threshold, "_full_jul")
prefix <- paste0(threshold*100, "_")

#> for filtering
my_dilfactors <- csv_groups$dilfactor %>% unique #> interested in all the common dilfactors (no filter)
# only 1-group --> 0.00025 0.00040 0.00400 0.00800 0.04000 0.25000 0.50000 0.80000
# 3- and 10-group samples have (global) success with the same 8 values.

# LOAD DATA --------------------------------------------------------------------
#> Load the 3 and 10 groups data
csv_groups <- read.csv(groups_file, fileEncoding = "UTF-8", sep = ",")
csv_groups["logdf"] <- log(csv_groups$dilfactor)

#> When no success, encode that as 0 instead of NA
csv_groups$success       <- csv_groups$success %>% replace(is.na(.), 0)
csv_groups$group_success <- csv_groups$group_success %>% str_replace_all("NA", "0")
csv_groups$group_success <- str_split(csv_groups$group_success, ";")

#> Now load the single-group data
csv_nogroups <- read.csv(nogroups_file, fileEncoding = "UTF-8", sep = ",")
csv_nogroups["logdf"] <- log(csv_nogroups$dilfactor)
#> When no success, encode that as 0 instead of NA
csv_nogroups$success       <- csv_nogroups$success %>% replace(is.na(.), 0)

#> IMPORTANT to remove data that we have no equivalent for in the multiple-group
#> results
csv_nogroups <- csv_nogroups[csv_nogroups$richness %in% unique(csv_groups$richness) & # 100 and 1000 only !!!
                               csv_nogroups$size == 10000 &
                               csv_nogroups$distrib == "lognorm", ]      # lognorm only !!!

#> Also add fields absent in the single-group dataset
csv_nogroups$group_success <- csv_nogroups$success %>% as.character() %>% as.list # same object class as the nicheN>1 ones
csv_nogroups$nicheN <- 1
csv_nogroups$nichedist <- "NoGroups"


csv <- bind_rows(csv_groups, csv_nogroups)
print(paste("Muestras procesadas:", nrow(csv)))

# for debugging:
# csv <- csv %>%
#   na.omit()
# print(paste("Eliminando los NA:", nrow(csv)))

csv <- csv[csv$dilfactor %in% my_dilfactors, ]
print(paste("Filtrando por dilution factor ", nrow(csv)))


# START ------------------------------------------------------------------------
#> EvenGroups; success (functional group no. vs. dilfactor effect)---

png(paste0(out_folder, "/0.2_success_v_nicheN.png"),
    width = 800, height = 600)
# f <- flexplot(success~dilfactor + nicheN, data=csv)
# plot(f)
# dev.off()

flexplot(success~dilfactor + nicheN,
         data=csv[csv$nichedist %in% c("EvenGroups", "NoGroups"), ])
dev.off()

# START ------------------------------------------------------------------------
#> EvenGroups; GROUP fixation (functional group no. vs. dilfactor effect)---
# we're going to plot the frequency of every nicheN-percentage of success combination.
csv$perc_success <- 100 * sapply(csv$group_success, function(x) sum(x!="0")) / csv$nicheN
csv <- csv %>%
  group_by(nicheN, perc_success) %>%
  mutate(countsuccess = n())

png(paste0(out_folder, "/0.2_group_success_v_nicheN_v_dilfactor_EVEN.png"),
    width = 800, height = 600)


ggplot(csv[csv$nichedist %in% c("EvenGroups", "NoGroups"), ],
       aes(x = nicheN, y = perc_success, color = countsuccess)) +
  geom_point() +
  labs(title = "Percentage of fixated groups (Even Groups)",
       x = "Niche N",
       y = "Percentage of Success (%)",
       color = "Occurrences") +
  facet_wrap(~ dilfactor) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_gradient(low = "lightgreen", high = "darkred") +
  scale_x_continuous(breaks = unique(csv$nicheN))

dev.off()

# START ------------------------------------------------------------------------
#> SkewedGroups; GROUP fixation (functional group no. vs. dilfactor effect)---
png(paste0(out_folder, "/0.2_group_success_v_nicheN_v_dilfactor_SKEWED.png"),
    width = 800, height = 600)

ggplot(csv[csv$nichedist %in% c("SkewedGroups", "NoGroups"), ],
       aes(x = nicheN, y = perc_success, color = countsuccess)) +
  geom_point() +
  labs(title = "Percentage of fixated groups (Uneven Groups)",
       x = "Niche N",
       y = "Percentage of Success (%)",
       color = "Occurrences") +
  facet_wrap(~ dilfactor) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_gradient(low = "lightgreen", high = "darkred") +
  scale_x_continuous(breaks = unique(csv$nicheN))

dev.off()

# START ------------------------------------------------------------------------
#> Check for extinction (Even; N 3 & 10)
csv$perc_extinction <- 100 * sapply(csv$group_sizes, function(x) str_count(x, " 0")) / csv$nicheN
csv <- csv %>%
  group_by(nicheN, perc_success) %>%
  mutate(countextinction = n())

png(paste0(out_folder, "/0.2_extinctions_v_nicheN_v_dilfactor_EVEN.png"),
    width = 800, height = 600)

ggplot(csv[csv$nichedist %in% c("EvenGroups"), ],
       aes(x = dilfactor, y = perc_extinction, color = countextinction)) +
  geom_point() +
  labs(title = "Percentage of fixated groups (Even Groups)",
       x = "Dilution factor",
       y = "Percentage of Success (%)",
       color = "Occurrences") +
  facet_wrap(~ nicheN) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_gradient(low = "lightgreen", high = "darkred") +
  scale_x_continuous(breaks = unique(csv$nicheN))

dev.off()

# START ------------------------------------------------------------------------
#> Check for extinction (Even; N 3 & 10)
csv$perc_extinction <- 100 * sapply(csv$group_sizes, function(x) str_count(x, " 0")) / csv$nicheN
csv <- csv %>%
  group_by(nicheN, perc_success) %>%
  mutate(countextinction = n())

png(paste0(out_folder, "/0.2_extinctions_v_nicheN_v_dilfactor_SKEWED.png"),
    width = 800, height = 600)

ggplot(csv[csv$nichedist %in% c("SkewedGroups"), ],
       aes(x = dilfactor, y = perc_extinction, color = countextinction)) +
  geom_point() +
  labs(title = "Percentage of fixated groups (Uneven Groups)",
       x = "Dilution factor",
       y = "Percentage of Success (%)",
       color = "Occurrences") +
  facet_wrap(~ nicheN) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_gradient(low = "lightgreen", high = "darkred") +
  scale_x_continuous(breaks = unique(csv$nicheN))

dev.off()
