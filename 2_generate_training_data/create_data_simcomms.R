### This script must be executed from the command line. It requires the
### following arguments to be passed at runtime:
# -i or --input: The path to the directory containing the simulation data
## simul_folder format => ../results_null_model/simcomms_2023-02-17/
## filenames format    => "SIMCOMM_SIMS_0.8_uniform_10000_sp_10_9/simul_9_t_86.csv"
##                                        3       4     5     7 8           11
## or                  => SIMCOMM_SIMS_"$g"_"$nichedist"_"$fdil"_"$distrib"_"$size"_sp_"$sp"_$sa/simul_<i>_t_$transfer.csv
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
### other community types (unlike the rhizosohere script)

################################################################################
library(tidyverse)
library(gridExtra)
library(parallel)
library(optparse)
library(vegan) # shannon + pielou custom function
library(DescTools) # gini
library(data.table)

parser <- OptionParser(option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL),
  make_option(c("-o", "--output_folder"), type="character", default=NULL),
  make_option(c("-n", "--numsamples"), type="integer", default=NULL),
  make_option(c("-p", "--pcgtable"), type="character", default=NULL),
  make_option(c("-c", "--cores"), type="integer", default=NULL),
  make_option(c("-f", "--fixation_threshold"), type="double", default=0.5,
              help = "Portion of simulations to be checked for fixation; 0.5 (50%) by default"),
  make_option(c("-s", "--successperc"), type="double", default=0.95,
              help = "Used to define 'success' within a transfer for a sample. Success is considered to have happened when fixation is reached at percN*100% of simulations for that sample and transfer.")
)
); opt <- parse_args(parser)

## Data ########################################################################
# The following lines define some simulation parameters, such as the fixation
# threshold, the portion of simulations to be checked for fixation, and the
# number of cores to be used. It also defines the output folder and the output
# name based on the input directory and the fixation threshold.
fixation_threshold <- opt$fixation_threshold

## input
cores <- opt$cores
simuls_folder  <- opt$input
num_of_samples <- opt$numsamples
pcgtable       <- opt$pcgtable
percN          <- opt$successperc # portion of simulations we're going to check
# for fixation; if percN*100% of simulations
# from one sample have reached fixation at
# transfer T, we define success at transfer T.

if (is.null(num_of_samples)) {
  stop("You must specify the number of samples")
}

## output
if (is.null(opt$output_folder)) {
    output_folder <- paste0("processed_data_simcomms_", fixation_threshold)
} else {
    output_folder <- opt$output_folder
}
if (!file.exists(output_folder)) {system(paste("mkdir -p", output_folder))}
output_name <- paste0("RESULT_", basename(simuls_folder))

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

# parsing of filenames is different...
read_simul_data_groups <- function(simuls_folder, num_of_samples) {
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
      "distrib" = split[6],
      "size" = split[7],
      "richness" = split[9],
      "dilfactor" = as.numeric(split[5]),
      "filename" = full_name,
      "sample" = paste(str_split(split[10], "/")[[1]][1],
                       split[6], split[7], split[9], sep = "_"),
      "nicheN" = as.numeric(split[3]),
      "nichedist" = split[4],
      "transfer" = as.numeric(split[13])
    ))
  }) %>% t %>% as_tibble()
}

record_success <- function(processed_data, percN=0.95, groups=FALSE) {
  message(paste0("(!) Success is considered to happen when fixation is reached at ", percN*100, "% of simulations"))
  if (groups) {
    combined_df <- do.call(rbind, processed_data$perc)
    # must substract 1: position 2 contains the result of the 1st cycle
    reached_fixation_at <- lapply(combined_df, function(g) which((g %>% as.numeric)>=percN)[1] - 1 )%>% as_tibble()
  } else {
    reached_fixation_at <- which((processed_data$perc %>% as.numeric)>=percN)[1] - 1
  }
  return(reached_fixation_at)
}

## Load all data and metadata ##################################################
options(scipen=10)
if (!is.null(pcgtable)) {
  simul_data <- read_simul_data_groups(simuls_folder, num_of_samples)
} else {
  simul_data <- read_simul_data(simuls_folder, num_of_samples)
}

message(paste("INPUT:", simuls_folder))

message(paste0("Read data!
----------
nrow>> ", nrow(simul_data), "
ncol>> ", ncol(simul_data),"

colnames>> ", paste(colnames(simul_data), collapse=', '), "

dilution factors>> ", paste(unique(simul_data$dilfactor), collapse=', '),"

community sizes >> ", paste(unique(simul_data$size), collapse=', '),"

samples>> ", paste(unique(simul_data$sample), collapse = ', ')))

if (!is.null(pcgtable)) {
  pcg_info <- fread(pcgtable)
  pcg_info <- pcg_info[1:(nrow(pcg_info)-1), ]
} else {
  pcg_info <- NULL
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

## $success: fixation for *all cores* in a sample for a given dilution factor
################################################################################
## !! -- important to assign each fixation value to its corresponding core !
message("Computing diversity measures..")
if (!is.null(pcgtable)) {
  metadata$group_success <- NA
  metadata$group_sizes   <- list()
}
metadata$success <- NA

metadata$final_size <- NA # $size is initial_size

metadata$evenness <- NA
metadata$shannon <- NA

metadata$gini <- NA

for (sa in unique(names(all_processed_data))) {
  for (d in unique(all_processed_data[[sa]]$dilfactor)) {
    selec <- all_processed_data[[sa]][unlist(all_processed_data[[sa]]$dilfactor)==d,][c("perc",
                                                                                         "dilfactor",
                                                                                         "transfer",
                                                                                         "size",
                                                                                         "filename",
                                                                                         "group_sizes")]
    successes <- record_success(selec, percN, groups=!is.null(pcgtable))

    # final_size
    metadata[metadata$sample==sa &
               metadata$dilfactor==d, "final_size"] <- selec[selec$transfer==1, "size"]

    # group_sizes (to check for extinction later; the last non-NA group sizes are saved)
    df <- selec$group_sizes %>%
      data.frame() %>%
      mutate(result = apply(., 1, function(row) {
        # Get the last value in the row that does not include "NA" substring
        non_na_values <- row[!grepl("NA", row)]
        # ---> there are 2 options: all NA or no NA. If only some elements are NA,
        # there's a problem with the groups
        if (ncol(non_na_values)){
          warning(paste0("There are NAs in the group_sizes since transfer 0. Something wrong with the groups or the simulations. (", simuls_folder, ")"))
        }
        last_non_na <- tail(non_na_values, 1)
        if (length(last_non_na) == 0) NA else last_non_na
      }))
    metadata[metadata$sample==sa &              
                 metadata$dilfactor==d, "group_sizes"] <- (df$result %>% paste(collapse = ";"))

    # group_success
    if (!is.null(pcgtable)) {
      metadata[metadata$sample==sa &
                 metadata$dilfactor==d, "group_success"] <- successes %>% paste(collapse = ";")
    }
    # success ==> if there are multiple groups, it happens when all the groups
    #             reach fixation (in percN of simulations)
    metadata[metadata$sample==sa &
               metadata$dilfactor==d, "success"] <- max(successes)

    # absolute abundances (initial)
    abs_ab <- data.table::fread(selec[selec$transfer==0,]$filename[[1]],
                                header = T,
                                drop = 1,
                                nrows = 1) %>% as_tibble()
    rel_ab <- abs_ab/selec[selec$transfer==0,]$size # initial size

    metadata[metadata$sample==sa & metadata$dilfactor==d,]$evenness     <- pielou(abs_ab)
    metadata[metadata$sample==sa & metadata$dilfactor==d,]$shannon  <- vegan::diversity(abs_ab)
    metadata[metadata$sample==sa & metadata$dilfactor==d,]$gini         <- Gini(abs_ab)
  }
}

metadata$transfer <- NULL

write.csv(x = metadata %>% apply(.,2,as.character), file = paste0(output_folder, "/", output_name, ".csv"), row.names = F)
