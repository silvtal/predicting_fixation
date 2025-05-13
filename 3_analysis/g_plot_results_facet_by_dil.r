# este script PLOTEA y también PARSEA las simulaciones, haciendo cálculos de diversidad
# las tablas resultantes incluyen datos sobre los tamaños de grupo y el éxito de grupo
# también registra la lista de TRANSFERS que contienen el último not NA ( <- fijación sin extinción, a menos que permitamos la extinción!)

# ---- plots ----
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
library(vegan)
library(dilgrowth)

setwd("~/repos/predicting_fixation/1_datasets/simulation_results/preparsed/")
source_folder <- "new_2025_simcomms_WITH_GROUPS/"

FIXATION_THRESHOLD <- 0.5

output_figure_folder <- paste0("../../../figures_groups/preparsed_", FIXATION_THRESHOLD)
missing_file = "~/temp_missing"

CORES <- 16

# ==============================================================================
# functions (diversity metrics, Gini) # from DescTools
{ # collapse this:
  
  #pielou custom function
  pielou <- function(vect) {
    bci_sub_h <- vegan::diversity(vect) # shannon
    bci_sub_s <- vegan::specnumber(vect) # richness
    bci_sub_j <- bci_sub_h / log(bci_sub_s) # pielou
    return(bci_sub_j[[1]])
  }
  
  # gini stuff:
  .NormWeights <- function(x, weights, na.rm=FALSE, zero.rm=FALSE, normwt=FALSE) {
    
    # Idea Henrik Bengtsson
    # we remove values with zero (and negative) weight. 
    # This would:
    #  1) take care of the case when all weights are zero,
    #  2) it will most likely speed up the sorting.
    
    if (na.rm){
      
      keep <- !is.na(x) & !is.na(weights)
      
      if(zero.rm)
        # remove values with zero weights
        keep <- keep & (weights > 0)
      
      x <- x[keep]
      weights <- weights[keep]
    } 
    
    if(any(is.na(x)) | (!is.null(weights) & any(is.na(weights))))
      return(NA_real_)
    
    n <- length(x)
    
    if (length(weights) != n) 
      stop("length of 'weights' must equal the number of rows in 'x'")
    
    # x and weights have length=0
    if(length(x)==0)
      return(list(x = x, weights = x, wsum = NaN))
    
    if (any(weights< 0) || (s <- sum(weights)) == 0) 
      stop("weights must be non-negative and not all zero")
    
    
    # we could normalize the weights to sum up to 1
    if (normwt) 
      weights <- weights * n/s
    
    return(list(x=x, weights=as.double(weights), wsum=s))
    
  }
  
  
  Mean <- function (x, ...)
    UseMethod("Mean")
  
  
  #' @rdname Mean
  #' @export
  Mean.Freq <- function(x, breaks, ...)  {
    sum(head(MoveAvg(breaks, order=2, align="left"), -1) * x$perc)
  }
  
  
  #' @rdname Mean
  #' @export
  Mean.default <- function (x, weights = NULL, trim = 0, na.rm = FALSE, ...) {
    
    if(is.null(weights)) {
      # use mean here instead of mean.default in order to be able to handle
      # mean.Date, mean.POSIXct etc.
      mean(x, trim, na.rm, ...)
      
    } else {
      if(trim!=0)
        warning("trim can't be set together with weights, we fall back to trim=0!")
      
      # # verbatim from stats:::weighted.mean.default
      # 
      # if (length(weights) != length(x))
      #   stop("'x' and 'w' must have the same length")
      # weights <- as.double(weights)
      # if (na.rm) {
      #   i <- !is.na(x)
      #   weights <- weights[i]
      #   x <- x[i]
      # }
      # sum((x * weights)[weights != 0])/sum(weights)
      
      # use a standard treatment for weights
      z <- .NormWeights(x, weights, na.rm=na.rm, zero.rm=TRUE)
      
      # we get no 0-weights back here...
      sum(z$x * z$weights) / z$wsum
      
    }
    
  }
  
  
  Gini <- function (x, weights = NULL, unbiased = TRUE, conf.level = NA, 
                    R = 10000, type = "bca", na.rm = FALSE) 
  {
    x <- as.numeric(x)
    if (is.null(weights)) {
      weights <- rep(1, length(x))
    }
    if (na.rm) {
      na <- (is.na(x) | is.na(weights))
      x <- x[!na]
      weights <- weights[!na]
    }
    if (any(is.na(x)) || any(x < 0)) 
      return(NA_real_)
    i.gini <- function(x, w, unbiased = FALSE) {
      w <- w/sum(w)
      x <- x[id <- order(x)]
      w <- w[id]
      f.hat <- w/2 + c(0, head(cumsum(w), -1))
      wm <- Mean(x, w)
      res <- 2/wm * sum(w * (x - wm) * (f.hat - Mean(f.hat, 
                                                     w)))
      if (unbiased) 
        res <- res * 1/(1 - sum(w^2))
      return(res)
    }
    if (is.na(conf.level)) {
      res <- i.gini(x, weights, unbiased = unbiased)
    }
    else {
      boot.gini <- boot(data = x, statistic = function(z, i, 
                                                       u, unbiased) i.gini(x = z[i], w = u[i], unbiased = unbiased), 
                        R = R, u = weights, unbiased = unbiased)
      ci <- boot.ci(boot.gini, conf = conf.level, type = type)
      res <- c(gini = boot.gini$t0, lwr.ci = ci[[4]][4], upr.ci = ci[[4]][5])
    }
    return(res)
  }
}
# ==============================================================================
# cálculos, parsing

# (1) take datos reorganizados (¡¡muchos!!)
folders <- list.files(source_folder, full.names = FALSE)
# folders <- folders[stringr::str_detect(folders, "3_SkewedGroups_")] # debug
# folders <- folders[stringr::str_detect(folders, "_lognorm_10000_sp_100_")] # debug
# folders <- readLines("~/temp_folders_missing") # DEBUG

# PARA CADA CARPETA (combinacion 1000 características + 1:30 comms):
diversity_metrics <- data.frame(matrix(ncol=0, nrow=0))
all_information <- data.frame(matrix(ncol=0, nrow=0))

for (folder in folders)  { # cada combi de características x30
  #> después leo las tablas
  simulations <- list()  # --> cada uno de estos vectores finales será un elemento de una lista de 200 elementos
  files <- list.files(paste0(source_folder, "/", folder), full.names = FALSE)
  for (f in files) {
    wholefn <- paste0(source_folder, "/", folder, "/", f)
    if (!file.exists(wholefn) || file.info(wholefn)$size == 0) {
      writeLines(wholefn, "~/temp_MISSING")
    } else { # only proceed if file is there and has at least the header !
      csv <- read.csv(wholefn)
      csv <- csv[stringr::str_detect(colnames(csv), "otu")] # So that it selects only otu info and not dummy columns
      # check NA rows...
      total_rows <- nrow(csv)
      not_NA <- which(rowSums(is.na(csv)) == 0)
      if (sum(dim(csv)) > 1) { # DEBUG: 1 or 0
        if (length(not_NA) < tail(not_NA, 1)) {
          warning(paste("There are missing time points for", f))
          }
        # else { # debug
          # (2) guardar solo el último timestep/row de cada iteración que no esté vacío (pre-extinción o post-fijación)
          simulations[[f]] <- csv[tail(not_NA, 1), ]
          non_empty_f <- f # for parsing the PCG table later without errors!
          # } # debug
        } else {
          warning(paste0(folder, "/", f, " is empty"))
      }
    }
  }
  
  if (length(simulations) == 0) {
    warning(paste0("No results at all for ", folder))
    # if no simulations at all for this combination of features, go to the next one
  } else {
    
    # ==========================================================================
    # Cálculo de medidas de diversidad
    # ==========================================================================
    
    # pielou
    evennesses <- mclapply(simulations,
                           FUN = function(iter) {
                             pielou(iter)
                           },
                           mc.cores = CORES) %>% bind_rows() %>% as.numeric
    mean_evenness <- evennesses %>% mean()
    evennesses <- evennesses %>% paste0(., collapse=";")
    
    # shannon
    shannons <- mclapply(simulations,
                           FUN = function(iter) {
                             vegan::diversity(iter)
                           },
                           mc.cores = CORES) %>% bind_rows() %>% as.numeric
    mean_shannon <- shannons %>% mean()
    shannons <- shannons %>% paste0(., collapse=";")
    
    # gini
    ginis <- mclapply(simulations,
                         FUN = function(iter) {
                           Gini(iter)
                         },
                         mc.cores = CORES) %>% bind_rows() %>% as.numeric
    mean_gini <- ginis %>% mean()
    ginis <- ginis %>% paste0(., collapse=";")
    
    # ==========================================================================
    # Miro la fijación y el éxito
    # ==========================================================================
    #> obtengo el metadata
    groups  <- as.numeric(sub(".*_SIMS_(\\d+)_.*", "\\1", folder))
    species <- as.numeric(sub(".*_sp_(\\d+)_.*", "\\1", folder))
    # richness <- as.numeric(sub(".*_sp_(\\d+)_\\d+$", "\\1", folder))
    group_distr <- sub(".*_(Even|Skewed)Groups_.*", "\\1", folder)
    dilution_factor <- as.numeric(sub(".*Groups_(\\d+\\.\\d+)_.*", "\\1", folder))
    distr <- sub(".*_(\\w+)_\\d+_sp_.*", "\\1", folder)
    size <- as.numeric(sub(".*_(\\d+)_sp_.*", "\\1", folder))
    community_type <- paste(groups, group_distr, species, sep = "_")
    comm_number <- as.numeric(sub(".*_sp_\\d+_(\\d+)$", "\\1", folder))
    
    #> cargo la tabla de PCGs necesaria
    pcg_table <- read_tsv(
      paste0("~/repos/predicting_fixation/1_datasets/PCGtables/",
             group_distr,
             "Groups_N",
             groups,
             "_",
             species,
             "sp.csv")
    )
    pcg_table <- pcg_table[1:(nrow(pcg_table) - 1 ), ] # remove last row (general info, not core info)
    
    #> configurar las capacidades de carga en base a la PCG_table
    # (solo para la última comunidad de este tipo, porque todas [a menos que
    # estén vacías, pero esas se saltan] tienen las mismas OTUs en el mismo orden)
    carrying_capacities <- rep(0, ncol(simulations[[non_empty_f]]))
    for (group in 1:nrow(pcg_table)) {
      leaves <- strsplit(pcg_table$Leaves[group], ";")[[1]]
      names(carrying_capacities)[colnames(simulations[[non_empty_f]]) %in% leaves] <- pcg_table$Core[group]
      carrying_capacities[colnames(simulations[[non_empty_f]]) %in% leaves]        <- (pcg_table$Average[group] * species) %>% round
    }
    # (3) mirar si hay fijación o no (separando por grupos)
    success_rates_folder <- mclapply(simulations,
                                        FUN = function(iter) {
                                            check_for_fixation(iter,
                                                               carrying_capacities = carrying_capacities,
                                                               fixation_at = FIXATION_THRESHOLD)
                                          },
                                        mc.cores = CORES) %>% bind_rows()
    # FIX COLUMN ORDER: SAME ORDER AS PCG TABLE
    success_rates_folder <- success_rates_folder[, pcg_table$Core]
    
    # (4) mirar si hay EXTINCIÓN o no (separando por grupos)
    group_extinctions <- mclapply(simulations,
                                     FUN = function(iter) {
                                       df <- data.frame(groups = names(carrying_capacities), 
                                                        abundances = as.numeric(iter))
                                       total_abundances <- aggregate(abundances ~ groups, df,
                                                                     sum) 
                                       return(setNames(total_abundances$abundances == 0, total_abundances$groups))
                                       
                                       },
                                     mc.cores = CORES) %>% unname %>% bind_rows()
    # FIX COLUMN ORDER: SAME ORDER AS PCG TABLE
    group_extinctions <- group_extinctions[, pcg_table$Core]
    

    
    # (a) se obtiene un porcentaje de éxito POR GRUPO para cada CARPETA
    # (b) se obtiene un porcentaje de éxito TOTAL para cada CARPETA (<- NOTA: no obtengo el transfer en el que ocurre cada cosa!!!)
    # (c) se obtiene un porcentaje de extinción POR GRUPO para cada CARPETA
    success_rates_folder <- bind_cols(community=folder,
                                      group_success_rate=paste0(colMeans(success_rates_folder), collapse = ";"), # (a)
                                      success_rate=(apply(success_rates_folder,  # (b)
                                                          1, FUN = function(row) all(row == TRUE)) %>% sum) / length(simulations),
                                      group_extinctions=paste0(colMeans(group_extinctions), collapse = ";"), # (c)
                                      num_groups=groups,
                                      richness=species,
                                      group_distribution=group_distr,
                                      dilution_factor=dilution_factor,
                                      distribution=distr,
                                      community_type=community_type,
                                      comm_number=comm_number
                                      )
    
    # ==========================================================================
    # Poner todo junto
    # ==========================================================================
    all_information <- bind_rows(all_information, bind_cols(success_rates_folder,
                                                        # record how many simulations we have for each folder
                                                        length=length(simulations),
                                                        # record diversity measures
                                                        mean_gini=mean_gini,
                                                        mean_evenness=mean_evenness,
                                                        mean_shannon=mean_shannon,
                                                        ginis=ginis,
                                                        evennesses=evennesses,
                                                        shannons=shannons
                                                        ))
  }
}

save.image("debug_ALL_information.RData") # debug
write_csv(all_information, "all_information.csv")

# ==============================================================================
# plots

#> realmente la "main" var es el community type, todo el grupo en group vars
my_study_var  <- "comm_number"
# my_study_var  <- "Interaction_value"

#> variables que definen el community_type:
my_subtitle_vars <- c(
  # "success_rate",
  "num_groups", 
  "group_distribution",
  # "distribution",
  "richness"
)

usar_my_other_var <- TRUE
my_other_var <- "dilution_factor"

# others: "Tipo", "Total_success",

# convertir a categóricas
all_information[["dilution_factor"]] <- as.factor(all_information[["dilution_factor"]])
all_information[["comm_number"]]     <- as.factor(all_information[["comm_number"]])

# ==============================================================================
# Generar los gráficos
for (ctype in unique(all_information$community_type)) {
  success_rates_ctype <- all_information[(all_information$community_type == ctype), ]
  grafico <- ggplot(success_rates_ctype, aes(x = !!sym(my_study_var),
                                   y = success_rate,
                                   color = !!sym(my_study_var))) +
    # {if (usar_my_other_var) 
    #   geom_point(aes(shape = !!sym(my_other_var)), size = 3) else 
    #     geom_point(size = 3)
    # } +
    geom_point(size = 3, position = "identity") +
    # geom_line(aes(group = interaction(!!sym(my_study_var), !!sym(my_other_var)))) +
    geom_hline(yintercept = 0.9, color = "red") + # shows success threshold
    {if (usar_my_other_var) facet_wrap(as.formula(paste("~", my_other_var))) else NULL} +
    labs(
      title = paste("Relación entre", my_study_var, "y Success Rate"),
      subtitle = paste0(
        sapply(my_subtitle_vars, function(var) {
          paste0(var, ": ", unique(success_rates_ctype[[var]]), ", ")
        }),
        collapse = ""
      ),
      x = my_study_var,
      y = "Success",
      color = "Community"
    ) +
    theme_minimal() +
    ylim(0, 1) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  ggsave(
    filename = paste0(file.path(output_figure_folder, ctype), ".png"),
    plot = grafico,
    width = 14, height = 10
  )
  save.image(paste0("../../../3_analysis/debug/", ctype))
}


# code for checking missing FOLDERS
system(paste0("rm ", missing_file, "; touch ", missing_file))
for (cn in unique(all_information$comm_number)) {
  for (ct in unique(all_information$community_type)) {
    for (dilf in unique(all_information$dilution_factor)) {
      success_rates_check <- all_information[(all_information$comm_number == cn),]
      success_rates_check <- success_rates_check[(success_rates_check$community_type == ct),]
      success_rates_check <- success_rates_check[(success_rates_check$dilution_factor == dilf),]
              if (nrow(success_rates_check) == 0) {
                splitct <- stringr::str_split(ct, "_")[[1]]
                for (i in 1:200) {
                  cat(paste0("new_SIMCOMM_SIMS_", splitct[1], "_", splitct[2], "Groups_", dilf, "_lognorm_10000_sp_", splitct[3], "_", cn,
                             "/community_", i, ".csv"),
                      file=missing_file, append = T)
                  cat("\n",
                      file=missing_file, append = T)
                }
        # cat(paste0("new_SIMCOMM_SIMS_", splitct[1], "_", splitct[2], "Groups_", dilf, "_lognorm_10000_sp_", splitct[3], "_", cn),
        #     file="~/temp_missing", append = T)
        # cat("\n",
        #     file="~/temp_missing", append = T)
      }
    }
  }
}
