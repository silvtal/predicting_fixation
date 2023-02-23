### Read filenames and create tibble for all the simul_data
read_simul_data <- function(simuls_folder) {
  filenames <- list.files(path = simuls_folder, 
                          recursive = TRUE,
                          pattern="*.csv",
                          full.names=FALSE)
  
  simul_data  <- mapply(filenames, FUN=function(name){
    split <- str_split(name, "_")[[1]]
    return(c("core" = split[5],
             "sample" = split[6],
             "transfer" = str_split(split[8], ".csv")[[1]][1],
             "dilfactor" = str_split(split[1], "/")[[1]][1],
             "filename" = paste0(simuls_folder, "/", name))
    )
  }) %>% t %>% as_tibble()
  
  simul_data$transfer <- as.numeric(simul_data$transfer)
  simul_data$dilfactor <- as.numeric(simul_data$dilfactor)
  simul_data <- arrange(simul_data, transfer)
  total_transfers <- max(simul_data$transfer)
  
  return(simul_data)
}



##### Calculate ################################################################
#pielou custom function
pielou <- function(vect) {
  bci_sub_h <- vegan::diversity(vect) # shannon
  bci_sub_s <- vegan::specnumber(vect) # richness
  bci_sub_j <- bci_sub_h / log(bci_sub_s) # pielou
  return(bci_sub_j[[1]])
}

#function for comparisons
elementwise.all.equal <- Vectorize(function(x, y) {isTRUE(all.equal(x, y))})



##### Graphical ################################################################
# function for X axis formatting
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

# theme for tables
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.8)),
  colhead = list(fg_params=list(cex = 1)),
  rowhead = list(fg_params=list(cex = 1))
)



##### Parsing ##################################################################
create_processed_data <- function (metadata, fixation_threshold, percN=0.95, cores=1) {
  ###---------------------------------------------------------------------------
  ### Read and parse files, include fields max_abunds + fixation_N + size
  ###---------------------------------------------------------------------------
  ## what's the most abundant OTU's relative abundance? (in each trajectory)
  # i will look at all transfers but I don't need to save more than one
  processed_metadata <- metadata
  
  ## Normalize to relative abundance, save into the tibble
  tmp <- mclapply(processed_metadata$filename,
                  FUN=function(n) {
                    data.table::fread(n, header = T, drop = 1) %>%
                      as_tibble()
                  }, mc.cores = cores
  )
  
  
  processed_metadata$size <- mclapply(tmp,
                                           FUN=function(n) {
                                             max(rowSums(n, na.rm = T))
                                           }, mc.cores = cores
  )
  
  
  processed_metadata$table <- mclapply(tmp,
                                       FUN=function(n) {
                                         # "replace": if some row is filled with "NA", that
                                         # means some OTU reached fixation in the
                                         # last transfer of that trajectory
                                         (n/rowSums(n)) %>% replace(is.na(.), fixation_threshold)
                                       }, mc.cores = cores
  )
  
  rm(tmp)
  
  ## $max_abunds (the most fixated OTU is the one with the highest abund)
  ## we also record who that OTU is
  processed_metadata$max_abuns <- mclapply(processed_metadata$table, FUN=function(t){
    mapply(t, FUN=max) %>% # if we have 500 simuls/trajectories,
      sort(decreasing = T)      # we will have 500 max abundances
  }, mc.cores = cores
  )
  
  ## Remove the full table; we're only going to save a much smaller table
  processed_metadata$table <- NULL
  processed_metadata$table <- mclapply(processed_metadata$filename,
                                       FUN=function(n) {
                                         abuntable <- data.table::fread(n, header = T, nrows = 1, drop = 1) %>%
                                           as_tibble() # only need the first line; all simuls are the same at time==0
                                         abuntable <- (abuntable/rowSums(abuntable)) %>% replace(is.na(.), fixation_threshold)
                                         return(abuntable)
                                       }, mc.cores = cores
  )
  
  ## $fixation_N: and what's the mean maximum abundance when we look at the top 
  ## N% simulations? That is, excluding the (100-N)% of simuls that are the
  ## farthest to fixation what's the mean progress towards fixation of one OTU? (<-- in that given Core, in that given Sample)
  processed_metadata$fixation_N <- mclapply(processed_metadata$max_abuns, FUN=function(t){
    head(t, trunc(percN*length(t))) %>%   #pick top N (default 0.95, 95%)
      mean()                          #get mean
  }, mc.cores = cores
  )
  
  # perc_N: how many simulations (looking at top N%) reached the fixation goal
  # i.e., how many entries in processed_metadata$fixation_N surpass fixation_threshold
  processed_metadata$perc_N <- mclapply(processed_metadata$max_abuns, FUN=function(ma) {
    ma <- head(ma, trunc(length(ma)*percN)) # aplico 95% !
    sum(ma>=fixation_threshold)/length(ma) %>% as.numeric
  }, mc.cores = cores
  )
  
  return(processed_metadata) # processed_data[[sa]]
}



create_sample_info <- function(processed_data, cores=1) {
  ###---------------------------------------------------------------------------
  ### Table with info (original abundance, richness) for each PCG
  ###---------------------------------------------------------------------------
  info_data <- processed_data[processed_data$transfer==0,][c("core", "table")]
  info_data <- distinct(info_data)
  n <- info_data$core
  info_data <- mclapply((1:nrow(info_data)),
                        FUN=function(row) {
                          cbind(
                            "core"= info_data$core[[row]],
                            sort(info_data$table[[row]], decreasing = TRUE)[1:min(5, ncol(info_data$table[[row]]))],
                            "richness"=ncol(info_data$table[[row]])
                          )
                        }, mc.cores=cores
  )
  names(info_data) <- n
  return(info_data) # sample_info[[sa]]
}



record_fixation <- function(processed_data) {
  reached_fixation_at <- c()
  corenames <- unique(processed_data$core)
  for (c in corenames) {
    reached_fixation_at <- c(reached_fixation_at, which((processed_data[processed_data$core==c,]$perc_N %>% as.numeric)>0.9)[1])
  } ; names(reached_fixation_at) <- corenames
  return(reached_fixation_at)
}


get_diffs_and_fixpoints <- function(sample_info, reached_fix) {
  ###---------------------------------------------------------------------------
  ### diffs: abundance difference between the two most abundant OTUs
  ### Y: transfer where fixation is reached for percN*100% of samples
  ###---------------------------------------------------------------------------
  diffs <- c()
  all_fixation_points <- c()
  names_all_points <- c()
  
  for (sa in names(sample_info)) {
    corenames <- reached_fix[[sa]] %>% names
    for (core in corenames) { # in the case there are multiple columns (like when
      # reached_fix is reached_fixation_at[sa])
      si <- sample_info[[sa]][[core]]
      diffs <- c(diffs,
                 ifelse(test = length(si)>3,
                        yes  = si[[2]] - si[[3]],
                        no   = si[[2]] - 0)
      ); names_all_points <- c(names_all_points, core)
    }
    all_fixation_points <- c(all_fixation_points, reached_fix[[sa]])
  }; names(all_fixation_points) <- names_all_points
  
  return(list(diffs, all_fixation_points))
}  


get_diffs_and_fixpoints_all_dils <- function(sample_info, reached_fix) {
  ###---------------------------------------------------------------------------
  ### diffs: abundance difference between the two most abundant OTUs
  ### Y: transfer where fixation is reached for percN*100% of samples
  ### --> this function is different to the previous one in that it's not only a
  ###     list of samples, but each sample is itself a list of dilutions
  ###---------------------------------------------------------------------------
  diffs <- c()
  all_fixation_points <- c()
  names_all_points <- c()
  
  for (sa in names(sample_info)) {
    for (c in 1:length(reached_fix[[sa]])) { # !!
      corenames <- reached_fix[[sa]][[c]] %>% names # !!
      for (core in corenames) { # in the case there are multiple columns (like when
        # reached_fix is reached_fixation_at[sa])
        si <- sample_info[[sa]][[core]]
        diffs <- c(diffs,
                   ifelse(test = length(si)>3,
                          yes  = si[[2]] - si[[3]],
                          no   = si[[2]] - 0)
        ); names_all_points <- c(names_all_points, core)
      }
      all_fixation_points <- c(all_fixation_points, reached_fix[[sa]][[c]]) # !!
    }; names(all_fixation_points) <- names_all_points
    
    return(list(diffs, all_fixation_points))
  }
}

##### Plots ####################################################################
plot_table <- function(info_data, theme) {
  ###---------------------------------------------------------------------------
  ### A plot of tables with info for all PCGs in a sample (uses data from
  ### sample_info)
  ###---------------------------------------------------------------------------
  a <- list()
  for (core in 1:length(info_data)) {
    a[[core]] <- tableGrob(info_data[[core]], theme = theme)
  }
  return(ggpubr::ggarrange(plotlist=a,
                           ncol = 1, nrow = length(info_data)))
}

plot_fixation_v_transferN <- function(processed_data, sa, percN,
                                      scale=NULL) {
  ###---------------------------------------------------------------------------
  # X: Transfer
  # Y: Mean fixation at 95% simuls; mean maximum abundance when we look at the
  # top 95% simulations. That is - excluding the 5% of simuls that are the furthest
  # to fixation, what's the mean progress towards fixation of one OTU? (<-- in that given Core, in that given Sample)
  ###---------------------------------------------------------------------------
  b <- ggplot(data = processed_data[[sa]],
              aes(x=as.numeric(transfer), y=as.numeric(fixation_N), color=core)) +
    geom_line() +
    ggtitle(sa, dil_factor) +
    xlab("Transfer") + ylab(paste0("Mean fixation at ", percN, "% simuls")) +
    ylim(0, 1) + theme(plot.title = element_text(hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5))
  if (!is.null(scale)) {
    b <- b + scale_color_manual(values=scale)        
  }
  return(b)
}

plot_fixation_perc_v_transferN <- function(processed_data, sa, percN,
                                           reached_fixation_at,
                                           scale=NULL) {
  ###---------------------------------------------------------------------------
  # X: transfer
  # Y: PERCENTAGE of simuls/iterations at each timestep/transfer that have reached
  #    fixation_threshold% fixation of one OTU
  ###---------------------------------------------------------------------------
  b <- ggplot(data = processed_data[[sa]],
              aes(x=as.numeric(transfer), y=as.numeric(perc_N), color=core)) +
    geom_line() +
    ggtitle(sa, dil_factor) +
    geom_vline(xintercept = reached_fixation_at[[sa]], size=0.5) +
    annotate("text", x=reached_fixation_at[[sa]], y=seq(from=0, length.out=length(reached_fixation_at[[sa]]), by=0.05), label=unique(processed_data[[sa]]$core), angle=0, size=3) +
    xlab("Transfer") + ylab("Perc simuls which reached fixation") +
    ylim(0, 1) + theme(plot.title = element_text(hjust = 0.5),
                       plot.subtitle = element_text(hjust = 0.5))
  
  if (!is.null(scale)) {
    b <- b + scale_color_manual(values=scale)        
  }
  return(b)
}


plot_abundiff_v_fixationpoint <- function(diffs, all_fixation_points,
                                          dil_factor,
                                          scale=NULL) {
  ###---------------------------------------------------------------------------
  ## X: abundance difference between the two most abundant OTUs (diffs)
  ## Y: transfer where fixation is reached for percN*100% of samples (all_fixation_points)
  ###---------------------------------------------------------------------------
  d <- ggplot(mapping=aes(x=diffs, y=all_fixation_points, color=names(all_fixation_points))) +
    geom_point() +
    ggtitle(paste0("Effect of abundance gap between top 2 OTUs on objective completion"),
            paste0("Dilution factor: ", dil_factor)) +
    xlab("Abundance gap between top 2 OTUs") +
    ylab("Perc simuls which reached fixation") +
    labs(color="Core")
  
  if (!is.null(scale)) {
    d <- d + scale_color_manual(values=scale)        
  }
  return(d)
}


plot_dilfact_v_fixation <- function(dflist, fixation_points,
                                    scale=NULL) {
  ###---------------------------------------------------------------------------
  ## X: vector of dilfactors
  ## Y: transfer where fixation is reached for percN*100% of samples
  ###---------------------------------------------------------------------------
  e <- ggplot(mapping=aes(x=dflist, y=fixation_points, color=names(fixation_points))) +
    geom_point() +
    ggtitle(paste0(sa, " - Effect of dilution factor on objective completion")) +
    xlab("Dilution factor (--> more diluted)") +
    ylab(paste0("Transfer where ",fixation_threshold, " fixation is reached (--> faster fixation)")) +
    labs(color="Core") +
    scale_x_continuous(trans=reverselog_trans(10), labels = function(x) sprintf("%.5f", x))
  
  if (!is.null(scale)) {
    e <- e + scale_color_manual(values=scale)        
  }
  return(e)
}



plot_initabund_v_fixation <- function(initabund, fix_points,
                                      title, ytitle,
                                      scale=NULL) {
  ###---------------------------------------------------------------------------
  ## X: vector of initial abundances
  ## Y: transfer where fixation is reached for percN*100% of samples
  ###---------------------------------------------------------------------------
  g <-ggplot(mapping=aes(x=initabund,
                         y=fix_points,
                         color=names(fix_points))) +
    geom_point() +
    geom_line() +
    ggtitle(paste0(title, " - Effect of initial abundance on objective completion")) +
    xlab("Initial abundance") +
    ylab(paste0("% simuls with ", ytitle, " fixation (--> faster fixation)")) +
    labs(color="Core") 
  
  if (!is.null(scale)) {
    g <- g + scale_color_manual(values=scale)        
  }
  return(g)
}


plot_fixation_vs_evenness <- function(evenness, fix_points,
                                      title, ytitle,
                                      scale){
  h <-ggplot(mapping=aes(x=evenness,
                         y=fix_points,
                         color=names(fix_points))) +
    geom_point() +
    ggtitle(paste0(title, " - Effect of initial evenness on objective completion")) +
    xlab("Evenness") +
    ylab(paste0("Transfer where ", ytitle, " fixation is reached (--> faster fixation)")) +
    labs(color="Core") +
    scale_x_continuous(labels = function(x) sprintf("%.5f", x))
  
  if (!is.null(scale)) {
    h <- h + scale_color_manual(values=scale)        
  }
  return(h)
}