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