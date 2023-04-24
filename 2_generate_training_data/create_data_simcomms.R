### This script must be executed from the command line. It requires the
### following arguments to be passed at runtime:
  # -i or --input: The path to the directory containing the simulation data
      ## simul_folder format => ../results_null_model/simcomms_2023-02-17/
      ## filenames format    => "SIMCOMM_SIMS_0.8_uniform_10000_sp_10_9/simul_9_t_86.csv"
      ##                                        3       4     5     7 8           11
  # -n or --numsamples: The number of simulation samples.
  # pcgtable: Table with information about groups, including relative abundance
  #           and OTU names for each group. If it's NULL, it's assumed there are
  #           no groups (==the community is a single group)

### It also has two optional arguments: pcgtable (path, only needed if there are
### multiple functional groups in our data) and successperc (see --help)

### Example:
### "~/R-4.0.5/bin/Rscript create_data_simcomms.R --input ../results_null_model/simcomms_"$date"/SIMCOMM_SIMS_"$fdil"_"$distrib"_"$size"_sp_"$sp"_ -n "$n
### (launch_create_data_simcomms, changes -i y -n in a loop)

################################################################################

### This script was tailored to the simcomm data, but it's actually suitable for
### other community types (unline the rhizosohere script)

################################################################################

library(stringr)
library(tidyverse)
library(gridExtra)
library(parallel)
library(optparse)
library(vegan) # shannon + pielou custom function
library(DescTools) # gini
library(data.table)

parser <- OptionParser(option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL),
  make_option(c("-n", "--numsamples"), type="integer", default=NULL),
  make_option(c("-p", "--pcgtable"), type="character", default=NULL),
  make_option(c("-s", "--successperc"), type="double", default=0.95,
              help = "Used to define 'success' within a transfer for a sample. Success is considered to have happened when fixation is reached at percN*100% of simulations for that sample and transfer.")
)
); opt <- parse_args(parser)

## Data ########################################################################
# The following lines define some simulation parameters, such as the fixation
# threshold, the portion of simulations to be checked for fixation, and the
# number of cores to be used. It also defines the output folder and the output
# name based on the input directory and the fixation threshold.
cores <- 16
fixation_threshold <- 0.5

## input
simuls_folder  <- opt$input
num_of_samples <- opt$numsamples
pcgtable       <- opt$pcgtable
percN          <- opt$successperc # portion of simulations we're going to check
                                  # for fixation; if percN*100% of simulations 
                                  # from one sample have reached fixation at
                                  # transfer T, we define success at transfer T. 

## output
output_folder <- paste0("processed_data_simcomms_", fixation_threshold)
output_name <- paste0("RESULT_", basename(simuls_folder))

if (is.null(num_of_samples)) {
  stop("You must specify the number of samples")
}

## Functions ###################################################################
source("simul_fixation_functions.R")
  # pielou()
  # create_processed_data(pcgtable)

read_simul_data <- function(simuls_folder, num_of_samples) {
  ## This function reads the filenames, but doesn't really open the files!
  filenames <- c()
  for (samplenum in 1:num_of_samples) {
    filenames <- c(filenames, list.files(path = paste0(simuls_folder, samplenum),
                                         pattern = "*.csv",
                                         recursive = TRUE,
                                         full.names = TRUE)
                   )
  }
  
  simul_data <- mapply(filenames, FUN=function(full_name) {
    last_dir <- tail(strsplit(dirname(full_name), "/")[[1]], 1)
    name  <- paste(last_dir, basename(full_name), sep = "/")
    name  <- str_sub(name, end=-5) # remove ".csv"
    split <- str_split(name, "_")[[1]]
    return(list(
      "distrib" = split[4],
      "size" = split[5],
      "richness" = split[7],
      "transfer" = as.numeric(split[11]),
      "dilfactor" = as.numeric(split[3]),
      "filename" = full_name,
      "sample" = paste(str_split(split[8], "/")[[1]][1],
                       split[4], split[5], split[7], sep = "_")
    ))
  }) %>% t %>% as_tibble()
}


record_success <- function(processed_data, percN=0.9, groups=FALSE) {
  message(paste0("(!) Success is considered to happen when fixation is reached at ", percN*100, "% of simulations"))
  if (groups) {
    combined_df <- do.call(rbind, processed_data$perc)
    reached_fixation_at <- lapply(combined_df, function(g) which((g %>% as.numeric)>=percN)[1]) %>% as_tibble()
  } else {
    reached_fixation_at <- which((processed_data$perc %>% as.numeric)>=percN)[1]
  }
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

if (!is.null(pcgtable)) {
  pcg_info <- fread(pcgtable)
  pcg_info <- pcg_info[1:(nrow(pcg_info)-1), ]
} else {
  message("No PCG table")
}

# processed_data + sample_info #################################################
################################################################################
metadata      <- simul_data[order(as.numeric(simul_data$transfer)),] # ordered by transfer
metadata$size <- as.numeric(metadata$size)

message("Creating 'processed_data'...")
all_processed_data <- list() # for plots [4] onwards
for (sa in unique(metadata$sample)) {
  # "all_processed_data":  fixated + perc + size
  all_processed_data[[sa]] <- create_processed_data(metadata[metadata$sample==sa,],
                                                    fixation_threshold,
                                                    cores=cores,
                                                    pcg_info = pcg_info)
}


## NEW COLUMNS; let's do a preliminary filter to reduce the number of rows to 1
## (one row per sample)
metadata <- metadata[metadata$transfer==0,]

## reached_fixation_at
## success: fixation for *all cores* in a sample for a given dilution factor
################################################################################
## !! -- important to assign each fixation value to its corresponding core !
message("Computing diversity measures..")
metadata$reached_fixation_at <- NA
metadata$group_success <- NA 
metadata$success <- NA 
metadata$final_size <- NA
metadata$initial_size <- NA

metadata$filt_shannon <- NA
metadata$raw_shannon <- NA
metadata$filt_even <- NA
metadata$raw_even <- NA

metadata$gini <- NA

for (sa in unique(names(all_processed_data))) {
  for (df in unique(all_processed_data[[sa]]$dilfactor)) {
    selec <- all_processed_data[[sa]][all_processed_data[[sa]]$dilfactor==df,][c("perc",
                                                                                 "dilfactor",
                                                                                 "transfer",
                                                                                 "size",
                                                                                 "filename")]
    successes <- record_success(selec, percN, groups=!is.null(pcgtable))
    
    # final_size
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$final_size <- selec[selec$transfer==1,]$size
    # initial_size
    initab <- selec[selec$transfer==0,]$size
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$initial_size <- initab
    
    # group_success
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$group_success <- successes %>% paste(collapse = ",")
    # success ==> if there are multiple groups, it happens when all the groups
    #             reach fixation (in percN of simulations)
    metadata[metadata$sample==sa & 
             metadata$dilfactor==df,]$success <- max(successes)    
    
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
      metadata[metadata$sample==sa & metadata$dilfactor==df,]$gini         <- Gini(abs_ab)
    } else {
      message(paste0("Vacío (sin fijación): sa ", sa, ", dil ", df))
    }
  }
}

message("saving...") # DEBUG
write.csv(x = metadata %>% apply(.,2,as.character), file = paste0(output_folder, "/", output_name, ".csv"), row.names = F)
message("saved!") # DEBUG
