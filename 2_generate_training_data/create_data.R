library(stringr)
library(tidyverse)
library(gridExtra)
library(parallel)
library("optparse")
for (dilfact in c('0.5', '0.25', '0.1', '0.05', '0.04', '0.025', '0.01', '0.008', '0.005', '0.004', '0.0025', '0.001')) {
  #for (dilfact in c('0.0005', '0.0004', '0.00025')) {
  dilfact <- as.character(dilfact)
  
  option_list <- list(
    make_option(c("-d", "--dilfact"), type="character", default=NULL)
  )
  
  # library(vegan) #pielou custom function
  
  #pielou custom function
  pielou <- function(vect) {
    bci_sub_h <- vegan::diversity(vect) # shannon
    bci_sub_s <- vegan::specnumber(vect) # richness
    bci_sub_j <- bci_sub_h / log(bci_sub_s) # pielou
    return(bci_sub_j[[1]])
  }
  
  ## Data ########################################################################
  cores <- 16
  ## Some info about the simulation parameters
  fixation_threshold <- 0.5
  simuls_folder <- paste0("../results_null_model/simuls001",dilfact) # en home!
  # simuls_folder <- paste0("~/allsimuls") # en home!
  # simuls_folder <- paste0("~/tomate_TEST_simuls_0.08")
  # simuls_folder <- "~/sim0.08/"
  
  ## Load all data and metadata ##################################################
  source("read_simul_data.R")
  # load("simul_data.RData")
  
  message("read data!") # DEBUG
  
  output_folder <- paste0("processed_data_", fixation_threshold)
  if (!file.exists(output_folder)) {system(paste("mkdir -p", output_folder))}
  
  # sample_info + first95 list ###################################################
  ################################################################################
  metadata <- simul_data
  first95 <- list(); namesfirst95 <- list(); sample_info <- list(); alist <- list()
  
  ## (1) sample_info
  for (sa in unique(metadata$sample)) {
    # let's also save the original abundances/richness for each OTU
    
    # sample_info[[sa]] will have dilfactor*corenumber fields, aprox.
    selec <- metadata[metadata$sample==sa & metadata$transfer==0,][c("core", "filename")]
    sample_info[[sa]] <-  lapply((1:nrow(selec)),
                                 FUN=function(row) {
                                   n <- selec$filename[[row]]
                                   firstline <- data.table::fread(file = n,header = T,nrows = 1,drop = 1) %>%
                                     as_tibble() # only need the first line; all simuls are the same at time==0
                                   firstline <- (firstline/sum(firstline)) %>% replace(is.na(.), fixation_threshold)
                                   
                                   cbind(
                                     "core"= selec$core[[row]],
                                     sort(firstline, decreasing = TRUE)[1:min(5, ncol(firstline))],
                                     "richness"=ncol(firstline)
                                   )
                                 }
    )
  }
  
  message("lets go for perc_95...") # DEBUG
  ## (2) max_abund, perc_95
  ## $max_abunds
  metadata$max_abuns <- lapply(metadata$filename, FUN= function(n){
    t <- data.table::fread(n, header = T, drop = 1) %>%
      as_tibble()
    t <- (t/rowSums(t)) %>% replace(is.na(.), fixation_threshold)
    apply(t, MARGIN=1, max) %>% # if we have 500 simuls/trajectories,
      sort(decreasing = T)      # we will have 500 max abundances
  })
  ## $perc_95
  metadata$perc_95 <- lapply(metadata$max_abuns, FUN=function(ma) {
    ma <- head(ma, trunc(length(ma)*0.95)) # aplico 95% !
    sum(ma>=fixation_threshold)/length(ma)}
  )
  
  ## (3) first95 FIELD
  metadata$first95   <- NA
  metadata$filt_even <- NA
  metadata$raw_even <- NA
  metadata$popul_size <- NA
  metadata$richness <- NA
  
  save.image(paste0("temp", dilfact)) # DEBUG
  message(unique(metadata$sample)) # DEBUG
  for (sa in unique(metadata$sample)) {
    # let's save the point where fixation is reached *for each PCG*; FOR EACH DF
    message(sa) # DEBUG
    for (df in unique(metadata$dilfactor)) {
      message(df) # DEBUG
      for (c in unique(metadata$core)) {
        message(c) # DEBUG
        selec <- metadata[metadata$sample==sa & metadata$core==c & metadata$dilfactor==df,]
        if (!!nrow(selec)) {
          metadata[metadata$sample==sa & metadata$core==c & metadata$dilfactor==df,]$first95 <- which((selec$perc_95 %>% as.numeric)>0.9)[1]
          
          message("now for (4) evenness, finally...") #	DEBUG
          
          ## (4) evenness. filtered. FIELD
          # coretib[[sa]][[c]] <- metadata[metadata$sample==sa & metadata$core==c,]
          abs_ab <- data.table::fread(selec$filename[[1]],
                                      header = T,
                                      drop = 1,
                                      nrows = 1) %>% as_tibble()
          ind <- rowSums(is.na(abs_ab)) == ncol(abs_ab); if (!ind) { #DEBUG            
            rel_ab <- abs_ab/sum(abs_ab)
            
            # filtro: OTUs de abundancia rel mayor al df; de estas calculo la evenness
            filt   <- abs_ab[rel_ab>=df]
            # filtro: que tras diluir queden al menos 10 bichos de ese otu
            # filt <- filt[filt>10]
            if (!!length(filt)) {
              metadata[metadata$sample==sa & metadata$core==c & metadata$dilfactor==df,]$filt_even  <- pielou(filt)
              metadata[metadata$sample==sa & metadata$core==c & metadata$dilfactor==df,]$raw_even   <- pielou(abs_ab)
              metadata[metadata$sample==sa & metadata$core==c & metadata$dilfactor==df,]$popul_size <- sum(abs_ab)
              metadata[metadata$sample==sa & metadata$core==c & metadata$dilfactor==df,]$richness   <- length(rel_ab)
            }
          } else {
            message(paste0("ojo, estaba vacía la df para: sa ", sa, " dil ", df, " core ", c)) #DEBUG
          }
          
        } else {
          message(paste0("selec está vacío para: sa ", sa, " dil ", df, " core ", c)) # DEBUG
        }
      }
    }
  }
  
  message("saving...") # DEBUG
  df <- metadata[metadata$transfer==0, ]%>% as.data.frame()
  df$max_abuns <- NULL
  df <- apply(df,2,as.character)
  write.csv(x = df, file = paste0(output_folder, "/simul_data_tomate001_", dilfact, ".csv"), row.names = F)
  message("saved...") # DEBUG
  
  ## 0 
  ## metadata$perc_95
  
  ## 1
  ## metadata$first95[metadata$transfer==0]
  
  
  ## mclr conversion of abundance tables.
  # library("metaMint")
  # row <- 1
  # df <-  paste0(simuls_folder, "/", metadata$filename[[row]])
  # df<- data.table::fread(file = df,header = T,drop = 1) 
  # df2 <- mclr(df, eps=0.1) #pseudocount==0.1, asumimos que los 0 son por undersampling. Otras cifras bajas también se sustituyen por 0.1.
  
} # end dilfact loop
