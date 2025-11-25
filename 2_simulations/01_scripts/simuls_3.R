# THIS SCRIPT
# we work with a community of <a reduced number of> members,
# <a reduced number of> functional groups, and a custom interactions table
#
# we print and plot the effects on growth probability in each dil-growth cycle
# + a comparison of the final relative abundances with and without interactions
# 
# NEW : because it checks for fixation in EACH group + doesn't simulate the
# no_interactions simulations
# --> adapted to 0g__plot_fixation_fixationPERgroups_gr_PREPARSED.R
# --> adapted to 0g__plot_extinction_extinctionPERgroups_gr_PREPARSED.R
# also: NEW does NOT INCLUDE PLOTS; plus includes data (success, extinctions and
# shannon) for EACH group separately 
# -->  %>% _per_group is a new function (similar to check_for_fixation)
# 
# OLD one creates ugly plots + the Total_success_with_NO_interactions
# which is unnecessary
# ==============================================================================
# ALGORITHM INFO
# in my_interactions, rows == affecting population and cols == affectED population
#
# This is the formula used to apply the interaction forces:
#   prob + (prob %*% as.matrix(my_interactions))
#
# Thus, the resulting vector is a horizontal vector which is the sum of all effects.
# 
# ==============================================================================
# libraries
# ==============================================================================
library("dilgrowth")
library("tidyverse")
library("parallel")
library("optparse")

# ==============================================================================
## params
# ==============================================================================
option_list <- list(
  # input params
  make_option(c("-d", "--dilution"), type = "numeric",
              default = NULL,
              help = "Dilution factor used"),
  make_option(c("--cores"), type = "integer",
              default = 8,
              help = "Number of cores to use in parallelization processes (mclapply). Default: 4.",
              metavar = "integer")
)

opt <- OptionParser(option_list = option_list) %>% parse_args

# data options
SAMPLE_ID <- 3
COMM_TYPE <- "EvenGroups_N3"
#COMM_TYPE <- "SkewedGroups_N10"
ABUN_TABLE <- "datos/lognorm_100sp_size10000.tsv"
PCG_TABLE <- paste0("datos/", COMM_TYPE, "_100sp.csv")
SUBSET <- 100
SIMULS_NAME <- COMM_TYPE # NOTA: aquí solía poner onlypos y cosas así
OUTPUT_FOLDER <- paste0("results_01_", COMM_TYPE, "_100sp")
# OUTPUT_FOLDER <- "results_01"
if (!file.exists(OUTPUT_FOLDER)) {
  system(paste0("mkdir -p ", OUTPUT_FOLDER, "/stats"))
}

# interaction options
INTER_VALUES <- c(1, 0.1, 0.01, 0.001, 0.0001)
INTER_FREQS <- c(1, 0.1, 0.01)
INTER_SIGNS  <- c(+1, -1, 2) # 2 == both pos and neg
ALLOW_MUTUAL <- FALSE
ALLOW_WITHIN_GROUP <- FALSE

# simulation options
#> 50% + 95%; 30% + 90%
FIXATION_THRESHOLD <- 0.5 # relative abundance an OTU needs to have to be considered "fixed"
SUCCESS_PERC <- 0.95 # percentage of iterations where fixation has to have happened to consider it a success
# 
DILUTIONS <- opt$dilution
if (is.null(DILUTIONS)) {
  stop("Dilution factor must be provided with --dilution or -d")
}
DILUTIONS <- c(DILUTIONS)  # convert to vector for loop compatibility
NO_OF_DIL <- 200
GROWTH_STEP <- 0.01
IS_GROWTH_STEP_A_PERC <- T
ALLOW_GROUP_EXTINCTIONS <- T
# 
NO_OF_SIMULS <- 100
CORES <- opt$cores

# ==============================================================================
# functions
# ==============================================================================
create_interact_table <- function(my_sample,
                                  INTER_VALUE,
                                  INTER_FREQ,
                                  ALLOW_MUTUAL,
                                  ALLOW_WITHIN_GROUP,
                                  carrying_capacities=NULL) {
  # my_sample : one-column dataframe
  # INTER_VALUE : magnitude of the interactions. Only a single possible value.
  # INTER_FREQ : percentage (0 - 1) of possible pairs that will have an interaction
  # ALLOW_MUTUAL : whether to allow mutual interactions or not (T or F)
  # ALLOW_WITHIN_GROUP : whether to allow interactions between members of the same group
  
  if (!ALLOW_WITHIN_GROUP) {
    if (is.null(carrying_capacities)) {
      error("carrying_capacities is needed for create_interact_table() when ALLOW_WITHIN_GROUP is FALSE")
    }
  }
  
  # 0) create empty interaction matrix
  my_interactions <- data.frame(matrix(0,
                                       nrow = nrow(my_sample),
                                       ncol = nrow(my_sample)),
                                row.names = rownames(my_sample))
  colnames(my_interactions) <- rownames(my_sample)
  
  # 1) identify possible interacting pairs (non within the same group, include mutual or not)
  if (ALLOW_MUTUAL) {
    pairs <- expand.grid(rownames(my_sample), rownames(my_sample)) %>% as.data.frame() 
  } else {
    pairs <- combn(rownames(my_sample), 2) %>% t %>%
      # !!!! If we leave it at that, from A B C we get only AB BC and AC. Never FROM C to A or to B. not totally random !!!!
      # !!!! So, we randomly swap elements with a probability of 50% for every element, making this list actually symmetrical !!!!
      apply(1, function(x) sample(x, length(x))) %>% t %>% 
      as.data.frame()
    

  }
  colnames(pairs) <- c("Var1", "Var2")
  
  for (group in unique(names(carrying_capacities))) {
    group_members <- rownames(my_sample)[names(carrying_capacities)==group]
    if (!ALLOW_WITHIN_GROUP) { # need to check here in case carrying_capacities is given but this is TRUE
      pairs_to_remove <- expand.grid(group_members, group_members) %>% as.data.frame() # regardless of the direction of the interaction!
      pairs <- anti_join(pairs, pairs_to_remove, by = c("Var1", "Var2"))
    }
  }
  
  # 2) select pairs according to frequency --> make TRUE or FALSE
  # NOTE: if ALLOW_MUTUAL, the total possible number of interaction is twice as many
  pairs <- pairs[sample(1:nrow(pairs), round(INTER_FREQ * nrow(pairs))),,drop=F]
  
  for (p in 1:nrow(pairs)) {
    row_name <- pairs[p, 1]
    col_name <- pairs[p, 2]
    my_interactions[row_name, col_name] <- INTER_VALUE
  }
  
  return(my_interactions)
}

calculate_shannon <- function(vec) {
  if (sum(vec) == 0) return(0)  # Handle cases with all zeros
  proportions <- vec / sum(vec)  # Calculate proportions
  shannon <- -sum(proportions * log(proportions), na.rm = TRUE)  # Shannon formula
  return(shannon)
}

calculate_shannon_per_group <- function (vec,
                                         group_assignations) # we need names(carrying_capacities) !!
{
    df <- data.frame(groups = group_assignations, 
                     abundances = as.numeric(vec))
    df <- df %>% group_by(groups) %>% 
      mutate(relative_abundances = abundances / sum(abundances)) %>% 
      summarise(shannon = -sum(relative_abundances * log(relative_abundances), na.rm = TRUE)) %>% # Shannon formula
      ungroup()
    
  return(df$shannon)
}

# ==============================================================================
# Load + parse sample
# ==============================================================================
my_sample <- read.csv(ABUN_TABLE,
                      # sep = "\t", row.names = 1)[SAMPLE_ID]
                      sep = "\t", row.names = 1)[SAMPLE_ID][0:SUBSET,, drop=F] # debug
abun_total <- sum(my_sample)

# Load + parse PCG table
pcg_table <- read_tsv(PCG_TABLE)
pcg_table <- pcg_table[1:(nrow(pcg_table) - 1 ), ] # remove last row (general info, not core info)

carrying_capacities <- rep(0, nrow(my_sample))
for (group in 1:nrow(pcg_table)) {
  leaves <- strsplit(pcg_table$Leaves[group], ";")[[1]]
  names(carrying_capacities)[rownames(my_sample) %in% leaves] <- pcg_table$Core[group]
  carrying_capacities[rownames(my_sample) %in% leaves]        <- (pcg_table$Average[group] * abun_total) %>% round
}
message(paste0("Simulating growth for groups: ", paste0(unique(names(carrying_capacities)), collapse = ", "), "."))

# ==============================================================================
# create empty dataframe for the results
# ==============================================================================
results_df <- data.frame(matrix(nrow = 0, ncol = 12))
colnames(results_df) <- c(
  "dilfactor", #> old: Dilution
  "Interaction_presence",
  "Interaction_value",
  "Interaction_sign",
  "Mutual_interactions",
  "Number of iterations",
  "Number of dilution cycles",
  "Fixation_threshold",
  "group_success_rate",
  "success_rate", #> old: Total_success_with_interactions
  # "Total_success_with_NO_interactions",
  "shannons",
  # "Diversity_with_NO_interactions"
  "group_extinctions"
)

results_df$group_success_rate <- results_df$group_success_rate %>% as.character
results_df$group_extinctions <- results_df$group_extinctions %>% as.character
results_df$shannons <- results_df$shannons %>% as.character

# ==============================================================================
# start
# ==============================================================================
for (DILUTION in DILUTIONS) {
  for (INTER_VALUE in INTER_VALUES) {
    for (INTER_SIGN in INTER_SIGNS) {
      for (INTER_FREQ in INTER_FREQS) {
                              # if ((INTER_SIGN == -1) & (INTER_VALUE==1) & (INTER_FREQ==1)) {
                              #   # in this case, no growth will be possible.
                              #   next
                              # } else {
        DESCRIPTION_INTER <- paste0(SIMULS_NAME, c("", "_nomutual_")[!ALLOW_MUTUAL + 1], INTER_SIGN, "_", INTER_VALUE, "_freq", 100 * INTER_FREQ)
        # ==============================================================================
        # simulations (dilgrowth functions)
        # ==============================================================================
        my_simuls_inter <- mclapply(X = 1:NO_OF_SIMULS,
                                    FUN = function(iter) {
                                      tryCatch({
                                      # create the interaction table
                                      my_interactions <- create_interact_table(my_sample,
                                                                               TRUE,
                                                                               INTER_FREQ,
                                                                               ALLOW_MUTUAL,
                                                                               ALLOW_WITHIN_GROUP,
                                                                               carrying_capacities)
                                      if (INTER_SIGN != 2) {
                                        my_interactions[my_interactions  == TRUE] <- INTER_VALUE * INTER_SIGN
                                      } else {
                                        true_indices <- which(my_interactions == TRUE, arr.ind = TRUE)  # Get indices of TRUE cells
                                        for (index in 1:nrow(true_indices)) {
                                          row <- true_indices[index, 1]
                                          col <- true_indices[index, 2]
                                          # Randomly assign either INTER_VALUE * +1 or INTER_VALUE * -1 to each TRUE cell
                                          my_interactions[row, col] <- sample(c(INTER_VALUE * +1, INTER_VALUE * -1), 1)
                                        }
                                      }
                                      
                                      
                                      # do the simulations
                                      trajectory <- simulate_timeseries(my_sample,
                                                                        carrying_capacities=carrying_capacities,
                                                                        interactions=my_interactions,
                                                                        abun_total=abun_total,
                                                                        
                                                                        dilution=DILUTION,
                                                                        no_of_dil=NO_OF_DIL,
                                                                        fixation_at=FIXATION_THRESHOLD,
                                                                        growth_step=GROWTH_STEP,
                                                                        is_growth_step_a_perc=IS_GROWTH_STEP_A_PERC,
                                                                        allow_group_extinctions=ALLOW_GROUP_EXTINCTIONS,
                                                                        
                                                                        logistic=FALSE,
                                                                        keep_all_timesteps=FALSE,
                                                                        force_continue=FALSE)
                                      print(paste("Simulation", iter, "finished"))
                                      return(trajectory)
                                      }, error = function(e) {
                                        warning(sprintf("Simulation %d failed. Error: %s",
                                                        iter, e$message))
                                        return(NULL)
                                      })
                                    }, mc.cores = CORES)
        
        # Calcular el porcentaje de simulaciones fallidas
        num_failed <- sum(sapply(my_simuls_inter, is.null))
        percent_failed <- (num_failed / NO_OF_SIMULS) * 100
        
        # Si hay simulaciones fallidas, imprimir el warning y saltar la iteración
        if (num_failed > 0) {
          warning(sprintf("Skipping iteration. DILUTION=%f, INTER_VALUE=%f, INTER_SIGN=%d, INTER_FREQ=%f. %.2f%% of simulations failed.",
                          DILUTION, INTER_VALUE, INTER_SIGN, INTER_FREQ, percent_failed))
          next
        }
        # ==============================================================================
        # Success rate comparison + extinctions + Shannon
        # ==============================================================================
        success_inter <- lapply(my_simuls_inter,
                                        FUN = function(iter) {
                                          check_for_fixation(iter,
                                                             carrying_capacities = carrying_capacities,
                                                             fixation_at = FIXATION_THRESHOLD)
                                          }
                                        ) %>% bind_rows()
        
        group_success_inter <- (apply(success_inter, 2, FUN = sum) / nrow(success_inter)) %>% as.character %>% paste0(., collapse = ";")
        total_success_inter <- (apply(success_inter, 1, FUN = function(row) all(row == TRUE)) %>% sum) / nrow(success_inter)
        
        # check extinction by group; code from: g_plot_results_facet_by_dil.r
        group_extinctions_inter <- lapply(my_simuls_inter,
                                      FUN = function(iter) {
                                        df <- data.frame(groups = names(carrying_capacities), 
                                                         abundances = as.numeric(iter))
                                        total_abundances <- aggregate(abundances ~ groups, df,
                                                                      sum) 
                                        return(setNames(total_abundances$abundances == 0, total_abundances$groups))
                                      }) %>%
          unname %>%
          bind_rows() %>%
          colMeans() %>% 
          paste0(., collapse = ";")
        
        # obtain Shannon diversity for each iteration in the list
        group_diversity_inter <- sapply(my_simuls_inter, FUN = function(vec) calculate_shannon_per_group(vec, names(carrying_capacities))) %>%
          apply(., 1, mean) %>% paste0(., collapse = ";")
        # total_diversity_inter <- sapply(my_simuls_inter, calculate_shannon) %>% mean
        
        # collect the current iteration into the results df
        iteration_result <- data.frame(
          DILUTION,
          INTER_FREQ,
          INTER_VALUE,
          INTER_SIGN,
          ALLOW_MUTUAL,
          NO_OF_SIMULS,
          NO_OF_DIL,
          FIXATION_THRESHOLD,
          group_success_inter, # NEW
          total_success_inter,
          group_diversity_inter, # NEW
          group_extinctions_inter # NEW
        ) %>% setNames(colnames(results_df))
        
        results_df <- bind_rows(results_df, iteration_result)
        
        # ==============================================================================
        # # DEBUG
        # filename <- sprintf("image_%s_%s_%s_%s.RData", DILUTION, INTER_VALUE, INTER_SIGN, INTER_FREQ)
        # save.image(file = filename)
                                  # }
      }
    }
  }

  # FOR EACH DILUTION, results_df grows in every loop iteration, then we save the whole table
  write.csv(results_df, file = paste0(OUTPUT_FOLDER, "/stats/stats_", DILUTION, "_", length(unique(names(carrying_capacities))), "groups_", SUBSET, "otus_", NO_OF_SIMULS, "iter.csv"), row.names = FALSE)
}


# debug
# carrying_capacities[my_simuls_inter[[1]]>0] %>% names %>% table

