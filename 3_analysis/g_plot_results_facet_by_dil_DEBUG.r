# Script similar al 03_: no vamos a filtrar por SUCCESS_PERC sino que vamos a
# plotear el success_perc de cada muestra de las 30 (200 simuls cada una). 
# > O, el success_perc medio de las 30 (6000 simuls de cada combinación de
# valores de:
# - factor de dilucion
# - tamaño de poblacion
# - numero de grupos
# - patrón de capacidades de carga
# - ...)
# Cada una de estas comunidades se someterá a 45 procesos simulados de
# dilución-crecimiento, comprendiendo 9 factores de dilución distintos y 5
# configuraciones de grupos funcionales (1, 3 o 10 grupos, con dos posibilidades
# para los casos de 3 o 10 grupos: capacidades de carga idénticas o diferentes).

library(tidyverse)
library(parallel)
library(dilgrowth)

setwd("~/repos/predicting_fixation/1_datasets/simulation_results/preparsed/")
source_folder <- "new_2025_simcomms_WITH_GROUPS/"

FIXATION_THRESHOLD <- 0.5

output_figure_folder <- paste0("../../../figures_groups/preparsed_", FIXATION_THRESHOLD)

CORES <- 16

# ==============================================================================
# (1) take datos reorganizados (¡¡muchos!!)
folders <- list.files(source_folder, full.names = FALSE)
# folders <- grep("_1$", folders, value = TRUE) # DEBUG

# PARA CADA CARPETA (combinacion 1000 características + 1:30 comms):
success_rates <- data.frame(matrix(ncol=0, nrow=0))

for (folder in folders)  { # cada combi de características x30
  
  #> después leo las tablas
  files <- list.files(paste0(source_folder, "/", folder), full.names = FALSE)
  
  for (f in files) {
    csv <- read.csv(paste0(source_folder, "/", folder, "/", f))[-c(1, 2)]
    # check NA rows...
    total_rows <- nrow(csv)
    not_NA <- which(rowSums(is.na(csv)) == 0)
    
    if (sum(dim(csv)) > 1) { # DEBUG: 1 or 0
      if (length(not_NA) < tail(not_NA, 1)) {
        warning(paste("There are missing time points for", f))
        write(paste0(folder, "/", f), "~/temp_MISSINGS", append = TRUE) # DEBUG
      } else {
        warning(paste0("OK: ", folder, "/", f))
      }
      # } # debug
    } else {
      warning(paste0(folder, "/", f, " is empty"))
      write(paste0(substr(folder, 5, nchar(folder)), "/", f), "~/temp_EMPTIES", append = TRUE) # DEBUG
    }
  }
}

# ==============================================================================

# # (1) take datos reorganizados (¡¡muchos!!)
# folders <- list.files(source_folder, full.names = FALSE)
# folders <- grep("_1$", folders, value = TRUE) # DEBUG
# 
# # PARA CADA CARPETA (combinacion 1000 características + 1:30 comms):
# success_rates <- data.frame(matrix(ncol=0, nrow=0))
# 
# for (folder in folders)  { # cada combi de características x30
#   write(folder, "~/temp_output", append = TRUE) # DEBUG
#   
#   #> después leo las tablas
#   simulations <- list()  # --> cada uno de estos vectores finales será un elemento de una lista de 200 elementos
#   files <- list.files(paste0(source_folder, "/", folder), full.names = FALSE)
#   write(files, "~/temp_output", append = TRUE) # DEBUG
#   
#   for (f in files) {
#     write(f, "~/temp_output", append = TRUE) # DEBUG
#     
#     csv <- read.csv(paste0(source_folder, "/", folder, "/", f))[-c(1, 2)]
#     # check NA rows...
#     total_rows <- nrow(csv)
#     not_NA <- which(rowSums(is.na(csv)) == 0)
#     write(length(not_NA), "~/temp_output", append = TRUE) # DEBUG
#     
#     if (sum(dim(csv)) > 1) { # DEBUG: 1 or 0
#       if (length(not_NA) < tail(not_NA, 1)) {
#         warning(paste("There are missing time points for", f))
#       }
#       # else { # debug
#       # (2) guardar solo el último timestep/row de cada iteración que no esté vacío (pre-extinción o post-fijación)
#       simulations[[f]] <- csv[tail(not_NA, 1), ]
#       non_empty_f <- f # for parsing the PCG table later without errors!
#       # } # debug
#     } else {
#       warning(paste0(folder, "/", f, " is empty"))
#     }
#   }
#   
#   if (length(simulations) == 0) {
#     warning(paste0("No results at all for ", folder))
#     # if no simulations at all for this combination of features, go to the next one
#   } else {
#     #> obtengo el metadata
#     groups  <- as.numeric(sub(".*_SIMS_(\\d+)_.*", "\\1", folder))
#     species <- as.numeric(sub(".*_sp_(\\d+)_.*", "\\1", folder))
#     group_distr <- sub(".*_(Even|Skewed)Groups_.*", "\\1", folder)
#     dilution_factor <- as.numeric(sub(".*Groups_(\\d+\\.\\d+)_.*", "\\1", folder))
#     distr <- sub(".*_(\\w+)_\\d+_sp_.*", "\\1", folder)
#     richness <- as.numeric(sub(".*_(\\d+)_sp_.*", "\\1", folder))
#     community_type <- sub("_(\\d+\\.\\d+)_.*_(\\d+)$", "", folder) # FIX DEBUG TODO
#     comm_number <- as.numeric(sub(".*_sp_\\d+_(\\d+)$", "\\1", folder))
#     
#     #> cargo la tabla de PCGs necesaria
#     pcg_table_file <- paste0("~/repos/predicting_fixation/1_datasets/PCGtables/",
#                              group_distr,
#                              "Groups_N",
#                              groups,
#                              "_",
#                              species,
#                              "sp.csv")
#     
#     write(pcg_table_file, "~/temp_output", append = TRUE) # DEBUG
#     
#     pcg_table <- read_tsv(pcg_table_file)
#     pcg_table <- pcg_table[1:(nrow(pcg_table) - 1 ), ] # remove last row (general info, not core info)
#     
#     #> configurar las capacidades de carga en base a la PCG_table
#     # (solo para la última comunidad de este tipo, porque todas [a menos que
#     # estén vacías, pero esas se saltan] tienen las mismas OTUs en el mismo orden)
#     carrying_capacities <- rep(0, ncol(simulations[[non_empty_f]]))
#     for (group in 1:nrow(pcg_table)) {
#       leaves <- strsplit(pcg_table$Leaves[group], ";")[[1]]
#       names(carrying_capacities)[colnames(simulations[[non_empty_f]]) %in% leaves] <- pcg_table$Core[group]
#       carrying_capacities[colnames(simulations[[non_empty_f]]) %in% leaves]        <- (pcg_table$Average[group] * richness) %>% round
#     }
#     # (3) mirar si hay fijación o no (separando por grupos)
#     # (4) se obtiene un porcentaje de éxito para cada CARPETA
#     success_rates_folder <- mclapply(simulations,
#                                      FUN = function(iter) {
#                                        check_for_fixation(iter,
#                                                           carrying_capacities = carrying_capacities,
#                                                           fixation_at = FIXATION_THRESHOLD)
#                                      },
#                                      mc.cores = CORES) %>% bind_rows()
#     success_rates_folder <- bind_cols(community=folder,
#                                       success_rate=(apply(success_rates_folder, 1, FUN = function(row) all(row == TRUE)) %>% sum) / length(simulations),
#                                       num_groups=groups,
#                                       richness=species,
#                                       group_distribution=group_distr,
#                                       dilution_factor=dilution_factor,
#                                       distribution=distr,
#                                       community_type=community_type,
#                                       comm_number=comm_number
#     )
#     
#     success_rates <- bind_rows(success_rates, bind_cols(success_rates_folder,
#                                                         length=length(simulations))) # record how many simulations we have for each folder.
#   }
# }