## this script is called from all create_data_ scripts; create_data_simcomms.R
## uses create_processed_data() and pielou()
## =============================================================================

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
create_processed_data <- function (metadata, fixation_threshold,
                                   cores=1,
                                   pcg_info=NULL) {
  ###---------------------------------------------------------------------------
  ### Read and parse files, include fields fixated + perc + size
  ###---------------------------------------------------------------------------
  ## it's important to sort the dataset by transfer, if it's not sorted already
  processed_metadata <- metadata[order(as.numeric(metadata$transfer)),]

  ## Read each sample/replicate X otu table
  ## Rows: OTUs; Columns: samples
  tmp <- mclapply(processed_metadata$filename,
                  FUN=function(n) {
                    t <- data.table::fread(file = n, header = T, drop = 1,) %>%
                      transpose %>%
                      as.data.frame()
                    rownames(t) <- paste0("otu", 1:nrow(t))
                    return(t)
                  }, mc.cores = cores
  )

  ## If there's a PCG table,
  ## - we add extra variable $group_sizes
  ## - $fixated and $perc will be a vector of values
  if (!is.null(pcg_info)) {
    ## [$group_sizes]
    processed_metadata$group_sizes <- mclapply(tmp, function(t) {
      mclapply(pcg_info$Core, FUN=function(g) {
        ## For each group, take note of its OTUs
        otus <- unlist(
          strsplit(as.character(pcg_info[grep(g, pcg_info$Core), "Leaves"]), ";")
        )
        ## Take note of their abundance in every instance
        apply(t[otus, ], 2, sum)
      }, mc.cores = 1
      ) %>% do.call(cbind, .) %>% apply(., MARGIN = 1, function(x) toString(as.list(x)))
    }, mc.cores = cores
    )

    ## [$fixated]
    ## (1) Normalize abundances
    ## (2) Replace NAs
    ## (3) check when fixation happened at each simulation
    ## >> Here, we simply ask ourselves: "for each recorded transfer, which
    ##    of our simulated communities has reached fixation?"
    processed_metadata$fixated <- mclapply(1:length(tmp), function(i) {
      ## For each tmp data.frame we'll return a string of comma-separated values.
      ## Each value indicates if fixation has been reached by each group in that
      ## sample and transfer.

      ## First we normalize the abundances
      t <- tmp[[i]]
      t <- sweep(t,2,colSums(t),`/`)
      t[is.na(t)] <- fixation_threshold

      group_df <-
        mclapply(pcg_info$Core, FUN=function(g) {
          ## For each group, take note of its OTUs
          otus <- unlist(
            strsplit(as.character(pcg_info[grep(g, pcg_info$Core), "Leaves"]), ";")
          )
          ## Check if that group reached fixation
          apply(t[otus, ], 2, function(x) any(x >= fixation_threshold))
        }, mc.cores = 1
        ) %>% do.call(cbind, .)

      return(group_df %>% apply(., MARGIN = 1, function(x) x) %>% as_tibble() %>% transpose())
      }, mc.cores = cores
    )

    ## [$perc]
    ## % of simulations that reached the fixation goal at that transfer
    ## i.e., how many communities surpass fixation_threshold
    processed_metadata$perc <- mclapply(processed_metadata$fixated,
                                        function(group_df) {
                                          # $perc
                                          return(group_df %>% apply(., 2,
                                                                    FUN=function(df){
                                                                      sum(df)/length(df)
                                                                      }) %>% as_tibble() %>% transpose()
                                          )
                                        }, mc.cores = cores
    )

  ## if no PCG table (no groups), $fixated and $perc will be a single value
  } else {
    ## [$fixated]
    processed_metadata$fixated <- mclapply(tmp,
                                         FUN=function(t) {
                                           t <- sweep(t,2,colSums(t),`/`)
                                           t[is.na(t)] <- fixation_threshold
                                           apply(t, 2, function(x) any(x >= fixation_threshold))
                                           }, mc.cores = cores
    )
    ## [$perc]
    processed_metadata$perc <- mclapply(processed_metadata$fixated,
                                        FUN=function(t){
                                          sum(t)/length(t)
                                        }, mc.cores = cores
    )
  }
  rm(tmp)
  return(processed_metadata) # processed_data[[sa]]
}



## For tomato rhizosphere samples, where we had one-core samples that could be
## grouped into multi-core communities.
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
