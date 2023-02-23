## For --input == results_null_model/simcomms_2023-02-17/
## NOTE: no cores. For redefining "success" as "success for all cores", copy code from create_data.R (tomato)
## NOTE: this script should be launched for different inputs separately (and then put all tables together)
## NOTE: the complete filenames are saved into the main data matrix, unlike create_data.R did

## will run for inputs:
## - 0.8 0.5 0.1 0.05 0.01 0.008 0.005 0.001 0.25 0.05 0.04 0.025 0.005 0.004 0.0025 0.0005 0.0004 0.00025
## - size 10000 10e-06
## - uniform lognormal
## - sp in 10 100 1000
## - sa 1..30 ?

## Simuls folder format::
## simul_folder => ../results_null_model/simcomms_2023-02-17/
## filenames    => "SIMCOMM_SIMS_0.8_uniform_10000_sp_10_9/simul_9_t_86.csv"
##                                 3       4     5     7 8           11

library(stringr)
library(tidyverse)
library(gridExtra)
library(parallel)
library(optparse)
library(vegan) # shannon + pielou custom function
library(DescTools) # gini

parser <- OptionParser(option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL),
  make_option(c("-n", "--numsamples"), type="integer", default=NULL)
)
); opt <- parse_args(parser)

## Data ########################################################################
cores <- 16
## Some info about the simulation parameters
fixation_threshold <- 0.5
percN <- 0.95 # portion of simulations we're going to check for fixation

## output
output_folder <- paste0("processed_data_simcomms_", fixation_threshold)
output_name <- paste0("RESULT_", opt$input)

## input
simuls_folder  <- opt$input
simuls_folder  <- "../results_null_model/simcomms_2023-02-17/SIMCOMM_SIMS_0.001_uniform_10000_sp_10_" # DEBUG 1
num_of_samples <- opt$num_of_samples
num_of_samples  <- 30 # DEBUG 1

if (is.null(num_of_samples)) {
  stop("You must specify the number of samples (30?)")
}
## Functions ###################################################################
## Adapted from simul_fixation_functions.R -- mais pour un seule core
source("simul_fixation_functions.R")

read_simul_data <- function(simuls_folder, num_of_samples) {
  filenames <- c()
  for (samplenum in 1:num_of_samples) {
    filenames <- c(filenames, list.files(path = paste0(simuls_folder, samplenum),
                                         pattern = "*.csv",
                                         recursive = TRUE,
                                         full.names = TRUE)
                   )
  }
  print(filenames) # DEBUG 1
  simul_data  <- mapply(filenames, FUN=function(full_name){
    name  <- paste0("SIM", str_split(full_name, "/SIM")[[1]][2]) # it's a bit ugly but: select only this part for splitting. Specific to this format.
    split <- str_split(name, "_")[[1]]
    return(c(#"core" = split[5],
      "distrib" = split[4],
      "size" = split[5],
      "richness" = split[7],
      "transfer" = str_split(split[11], ".csv")[[1]][1],
      "dilfactor" = split[3],
      "filename" = full_name,
      "sample" = paste(
                   str_split(split[8], "/")[[1]][1],
                   split[4],
                   split[5],
                   split[7], sep = "_"
                 )
      )
   )
  }) %>% t %>% as_tibble()
  simul_data$transfer <- as.numeric(simul_data$transfer)
  simul_data$dilfactor <- as.numeric(simul_data$dilfactor)
  simul_data <- arrange(simul_data, transfer)
  total_transfers <- max(simul_data$transfer)
  return(simul_data)
}

record_fixation <- function(processed_data) {
  reached_fixation_at <- which((processed_data$perc_N %>% as.numeric)>0.9)[1]
  return(reached_fixation_at)
}


## Load all data and metadata ##################################################
options(scipen=10)
simul_data <- read_simul_data(simuls_folder, num_of_samples)
message(paste0("Read data!
----------               
nrow>> ", nrow(simul_data), "
ncol>> ", ncol(simul_data),"

colnames>> ", paste(colnames(simul_data), collapse=', '), "

dilution factors>> ", paste(unique(simul_data$dilfactor), collapse=', '),"

community sizes >> ", paste(unique(simul_data$size), collapse=', '),"

samples>> ", paste(unique(simul_data$sample), collapse = ', ')))

if (!file.exists(output_folder)) {system(paste("mkdir -p", output_folder))}


# processed_data + sample_info #################################################
################################################################################
metadata <- simul_data

message("creating 'processed_data'...") # DEBUG
all_processed_data <- list() # for plots [4] onwards
for (sa in unique(metadata$sample)) {
  # "all_processed_data":  max_abunds + fixation_N...
  all_processed_data[[sa]] <- create_processed_data(metadata[metadata$sample==sa,],
                                                    fixation_threshold,
                                                    percN = percN, cores=cores)
  all_processed_data[[sa]] <- arrange(all_processed_data[[sa]]) # importante ordenar
}


## NEW COLUMNS; let's do a first filter to reduce the number of rows
metadata <- metadata[metadata$transfer==0,]

## reached_fixation_at
## success: fixation for *all cores* (there's one core only) in a sample for a given dilution factor
################################################################################
## !! -- important to assign each fixation value to its corresponding core !
metadata$reached_fixation_at <- NA
metadata$success <- NA 
metadata$final_size <- NA
metadata$initial_size <- NA

metadata$filt_shannon <- NA
metadata$raw_shannon <- NA
metadata$filt_even <- NA
metadata$raw_even <- NA

# metadata$sorensen <- NA
metadata$gini <- NA

for (sa in unique(names(all_processed_data))) {
  for (df in unique(all_processed_data[[sa]]$dilfactor)) {
    selec <- all_processed_data[[sa]][all_processed_data[[sa]]$dilfactor==df,][c("perc_N",
                                                                                 "dilfactor",
                                                                                 "transfer",
                                                                                 "size",
                                                                                 "filename")]
    fix <- record_fixation(selec)
    # final_size
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$final_size <- selec[selec$transfer==1,]$size
    # initial_size
    initab <- selec[selec$transfer==0,]$size
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$initial_size <- initab
    
    # reached_fixation_at
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$reached_fixation_at <- fix
    # success
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$success <- fix    
    
    # absolute abundances (initial)
    abs_ab <- data.table::fread(selec[selec$transfer==0,]$filename[[1]],
                                header = T,
                                drop = 1,
                                nrows = 1) %>% as_tibble()
    rel_ab <- abs_ab/initab
    filt   <- abs_ab[rel_ab>=df] # filter: discard those OTUs that are going to
                                 # get diluted out by that dilution factor
    
    # filtro: que tras diluir queden al menos 10 bichos de ese otu
    # filt <- filt[filt>10]
    if (!!length(filt)) {
      metadata[metadata$sample==sa & metadata$dilfactor==df,]$filt_even    <- pielou(filt)
      metadata[metadata$sample==sa & metadata$dilfactor==df,]$raw_even     <- pielou(abs_ab)
      metadata[metadata$sample==sa & metadata$dilfactor==df,]$filt_shannon <- vegan::diversity(abs_ab)
      metadata[metadata$sample==sa & metadata$dilfactor==df,]$raw_shannon  <- vegan::diversity(filt)
      # metadata[metadata$sample==sa & metadata$dilfactor==df,]$sorensen     <- vegan::betadiver(abs_ab[!!abs_ab], "w") # <- dissimilarity between 2 vectors. Not one (comm) or many (OTUs)
      metadata[metadata$sample==sa & metadata$dilfactor==df,]$gini         <- Gini(abs_ab)
    } else {
      message(paste0("Vacío (sin fijación): sa ", sa, ", dil ", df)) # DEBUG
    }
  }
}

message("saving...") # DEBUG
df <- metadata[metadata$transfer==0, ] %>% as.data.frame()
df <- apply(df, 2, as.character)
write.csv(x = df, file = paste0(output_folder, "/", output_name, ".csv"), row.names = F)
message("saved!") # DEBUG

