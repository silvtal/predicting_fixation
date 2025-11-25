start_time <- Sys.time()
## =============================================================================
##                              ÚLTIMA VERSIÓN
## =============================================================================
## En la versión anterior (v3), había un wrapper cuyo input era un fichero con
## información de cada PCG (sus hojas y sus abundancias relativas finales), y
## este script era ejecutado por separado para cada PCG.
##
## En esta nueva versión, las simulaciones ocurren de forma paralela para todos
## los PCGs. Ese fichero va a ser pues el input directo de este script. Durante
## las simulaciones, en cada iteración, se va a tener en cuenta el crecimiento
## proporcional que debe tener cada organismo
##
## minor changes: tax/pcgc variables were unneccessary and removed; abuntable
## is now formated with character colnames + no taxonomy column.
##
## =============================================================================
##                                SIMULACIÓN
## =============================================================================
## El UNTB consiste en una serie de comunidades locales neutrales que reciben
## inmigración de una metacomunidad también neutral. El modelo original
## sustituye, en cada generación/timestep, DE MEDIA una vez a cada organismo.
## Esto significa que en cada timestep habrá N iteraciones, siendo N el número
## total de bichos, en las que se escogerá al azar un individuo que morirá y
## será sustituido por un nuevo individuo de una especie ya presente (división
## celular) o bien de una no presente en la comunidad local pero sí en la
## metacomunidad (inmigración).
##
## En mi modelo hay un par de cambios: en primer lugar se trata de comunidades
## aisladas, no hay metacomunidad ni inmigración (son aisladas porque estoy
## replicando experimentos de laboratorio); en segundo lugar, cada vez que se
## alcanza cierta abundancia total hago una dilución, simulando los transfers
## del experimento del wet-lab. Esto quiere decir que no hago N iteraciones o
## sustituciones, sino tantas como haga falta para alcanzar la misma abundancia
## total que había antes de la dilución. Hay 12 transfers, contando el inóculo
## inicial.
##
## Las abundancias relativas en cada timestep se guardan en una matriz y
## posteriormente se hace la media para cada OTU en cada momento y se calculan
## la richness, evenness y diversidad, medidas de las cuales también se calcula
## una media.
##
## =============================================================================
## disclaimer
## =============================================================================
## Hacemos unas asunciones que hay que entender bien:
##
##   1) Asumimos que tras cada transfer se ha alcanzado ya la proporción final
##   conocida de cada PCG. Obviamente en la vida real no tiene por qué ser así,
##   por ejemplo: que tras 7 transfers haya 75% de Enterobacteriaceae y 25% de
##   Pseudomonadaceae pero que tras la 3ª haya todavía 55% y 45%.
##   [En Talavera et al. 2022 se tuvo en cuenta simulando a partir de una
##   comunidad inicial ya estable]
##
##   2) Los bichos de cada PCG crecen a la par, al mismo ritmo (pero proporcio-
##   nalmente). Esto quiere decir que no se va a dar el caso en el que (dentro
##   de uno de los ciclos de crecimiento) primero crezca PCG1 y luego PCG2. Sin
##   embargo, puede que en la vida real pase justo eso: que crezca primero el
##   más abundante y luego el otro, consumiendo sus productos de desecho, sea el
##   que crece. Asumimos que este orden no tiene mucho impacto, ni siquiera al
##   introducir las interacciones durante el sampleo, porque las proporciones
##   son fijas.
##
##   3) Asumimos que la abundancia total de cada PCGs es la misma para todas las
##   simulaciones, sin ruido
## =============================================================================

# ===========
# mis datos [glc]
# ===========
library("tidyverse", exclude = c("purrr::transpose"))
library("dilgrowth") #> includes data.table::transpose
library("gsubfn")
library("parallel")
library("optparse")
library("dplyr") # para bind_rows tras paralelización

option_list <- list(
  # input params
  make_option(c("-a", "--abuntable"), type = "character",
              default = NULL,
              help = "Abundance table"),
  make_option(c("-p", "--pcgtable"), type = "character",
              default = NULL,
              help = "PCG table; table with information with each PCG, output from BacterialCore.py; required fields: Core (PCG name), Average (relative abundance), Leaves."),
  make_option(c("--skip_lines_abuntable"), type = "integer",
              default = 0,
              help = "How many lines of abuntable should we skip (read.csv()'s 'skip'). Should be 1 for Qiime 1 output, for instance"),
  make_option(c("-s", "--sample"), type = "character",
              default = NULL,
              help = "Sample name in the abundance table"),
  make_option(c("-i", "--interactions"), type = "character",
              default = NULL,
              help = "Path to a .csv interactions table"),
  make_option(c("--selected_species"), type = "character",
              default = NULL,
              help = "semicolon-separated list of species to select from the abundance table"),
  make_option(c("--dilution"), type = "character",
              default = "8 * 10 ** (-3)",
              help = "Dilution before each simulated transfer. Can be an equation (e.g. 10**(-3))"),
  make_option(c("--no_of_dil"), type = "integer",
              default = 12,
              help = "Number of dilutions/transfers"),
  make_option(c("--fixation_at"), type = "double",
              default = 1,
              help = "The simulations stop when only one bug is fixed at 100% [fixation_at=1], since it's meaningless to keep going when only one species is left. Nevertheless, we can consider an OTU is fixed when it has reached a lower relative abundance (e.g. fixation_at=0.95) and stop earlier."),
  make_option(c("--save_all"), type = "logical",
              default = FALSE,
              help = "save all intermediate states of the simulations? default FALSE"),
  make_option(c("--growth_step"), type = "double",
              default = 1,
              help = "How many bugs are born on each iteration. 1 by default"),
  make_option(c("--is_growth_step_a_perc"), type = "logical",
              default = FALSE,
              help = "If FALSE, growth_step is taken as a fixed value, so the step will always be the same. If TRUE, it's read as a percentage - the step will be changed proportionally to the community size. If growth_step is 0.02, 2% of the members in the community will grow the next iteration. FALSE by default"),
  make_option(c("--logistic"), type = "logical",
              default = FALSE,
              help = "¿Should growth be logistic? If TRUE, the growth_step/growth rate of each species will change over the course of the simulation, in a logistic growth manner (see growth_log)."),
  make_option(c("--no_of_simulations"), type = "integer",
              default = 1,
              help = "Number of simulations"),
  make_option(c("--allow_group_extinctions"), type = "logical",
              default = TRUE,
              help = "If TRUE, simulations will continue even if one or more groups go extinct, and the function will try to reach fixation in the rest of the groups. Only applicable when carrying_capacities is not NULL (when there are multiple functional groups) Also, this being FALSE does NOT affect groups that were not in the community from the start (if there are missing groups from the start, there will be a warning). default TRUE"),
  make_option(c("--cores"), type = "integer",
              default = 1,
              help = "Number of cores to use in parallelization processes (mclapply). Default: 4.",
              metavar = "integer"),
  # output params
  make_option(c("--outputname"), type = "character",
              default = "out",
              help = "name for output .csv file"),
  make_option(c("--outdir"), type = "character",
              default = "results",
              help = "output directory")
)

parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

abuntable  <- opt$abuntable  # e.g. "./original_100percArbol/Tree/0.99/table.from_biom_0.99.txt"
skip_lines_abuntable <- opt$skip_lines_abuntable
pcgtable <- opt$pcgtable
s <- opt$sample      # e.g. "X2", "sa2"...
interactions <- opt$interactions
selected_species <- opt$selected_species

outdir <- opt$outdir # e.g. "my_neutral_model_v2_test_16_simuls_8_cores"
outputname <- opt$outputname # e.g. "X2_rep4"
cores <- opt$cores  # e.g. 16

save_all  <- opt$save_all
growth_step <- opt$growth_step
is_growth_step_a_perc <- opt$is_growth_step_a_perc
allow_group_extinctions <- opt$allow_group_extinctions

if (is_growth_step_a_perc && (growth_step == 1)) {
  message("WARNING: growth_step has been defined as a percentage of value 1. All organisms will duplicate every iteration. If you meant for growth_step to be a fixed value, set is_growth_step_a_perc to FALSE instead.")
}

logistic <- opt$logistic

dilution <-  eval(parse(text=opt$dilution)) # e.g. 8 * 10 ** (-3)
no_of_dil <- opt$no_of_dil
fixation_at <- opt$fixation_at

no_of_simulations <- opt$no_of_simulations

# read abundance table
exp <- read.csv(
  abuntable,
  sep = "\t",
  skip = skip_lines_abuntable,
  row.names = 1,
  check.names = FALSE
)
species_are_rows <- TRUE; if (!species_are_rows) {exp <- .my_transpose(exp)}
exp <- exp[colnames(exp) != "taxonomy"] # remove taxonomy column if present
colnames(exp) <- as.character(colnames(exp)) # in case sample names are numbers

# read PCG table
if (!(is.null(pcgtable))) {
  pcg_table <- read.csv(pcgtable, sep="\t")
  pcg_table <- pcg_table[1:(nrow(pcg_table)-1),] # remove last row (general info, not core info)
  pcg_table <- pcg_table[c("Core", "Average", "Leaves")]
  pcg_table$Average <- as.numeric(pcg_table$Average)
}

# read interactions table
if (!(is.null(interactions))) {
  interactions <- read.csv(file = interactions, row.names = 1, check.names = F)
}

# create output dir
system(paste("mkdir -p", outdir))

# =============
# create counts
# =============
# s is the original sample, exp[s] is a vector of abundances
# selected_species is a vector of species from that vector
if (is.null(selected_species)) {
  counts <- exp[s]
} else {
  selected_species <- strsplit(selected_species, ";")[[1]]
  counts <- exp[selected_species, s, drop=F]
}

# ==========
# simulation
# ==========
## results will be stored here:
final_abund <- list()

# initial/final total abundance
abun_total <- sum(counts)

# check first if there's anything to simulate
if (abun_total == 0) {
  message(paste0("There are no detectable OTUs in initial sample ",
                 s, ". Moving to next PCG..."))
} else if (abun_total * dilution < 3) {
  stop("EXIT: 3 or less bugs will be left after diluting! Consider changing your dilution factor.")

} else {
  if (!(is.null(pcgtable))) {
    # parse PCG table if everything's OK
    abun_others <- abun_total * (1 - sum(pcg_table$Average))
    carrying_capacities <- rep(round(abun_others), nrow(counts))
    names(carrying_capacities) <- rep("others", nrow(counts))
    for (group in 1:nrow(pcg_table)) {
      leaves <- strsplit(pcg_table$Leaves[group], ";")[[1]]
      names(carrying_capacities)[rownames(counts) %in% leaves] <- pcg_table$Core[group]
      carrying_capacities[rownames(counts) %in% leaves]        <- (pcg_table$Average[group] * abun_total) %>% round
    }
    message(paste0("Simulating growth for groups: ", paste0(unique(names(carrying_capacities)), collapse = ", "), "."))
  } else {
    carrying_capacities=NULL
    message("No PCG table provided; simulating growth without groups")
  }
  # start simulations
  abund_temp <- mclapply(X = 1:no_of_simulations,
                         FUN = function(iter) {
                           trajectory <- simulate_timeseries(counts,
                                                             carrying_capacities = carrying_capacities,
                                                             interactions = interactions,
                                                             logistic = logistic,
                                                             dilution = dilution,
                                                             no_of_dil = no_of_dil,
                                                             fixation_at = fixation_at,
                                                             abun_total = abun_total,
                                                             growth_step = growth_step,
                                                             is_growth_step_a_perc = is_growth_step_a_perc,
                                                             keep_all_timesteps = save_all,
                                                             allow_group_extinctions = allow_group_extinctions,
                                                             force_continue = FALSE)
                           print(paste("Simulation", iter, "finished for", s))
                           return(trajectory)
                         }, mc.cores = cores)

  # save data
  message(paste0("Saving data for sample ", s, "..."))
  if (save_all == T) {
    for (timepoint in 1:(no_of_dil + 1)) {
      # pick all rows number "timepoint" from all the lists in all_abund
      # "temp" refers to a single timepoint
      temp <- lapply(abund_temp, FUN = function(traj) {
        return(traj[timepoint, , drop = F] %>% as.data.frame)}) %>%
        bind_rows %>% as_tibble # one file per timepoint

      write.csv(temp,
                file = paste0(outdir, "/simul_", outputname, "_t_", timepoint - 1, ".csv"))
    }
  } else {
    final_abund <- as.data.frame(bind_rows(abund_temp))
    write.csv(final_abund,
              file = paste0(outdir,"/simul_", outputname, ".csv"))
  }
}

end_time <- Sys.time()
print(end_time - start_time)

