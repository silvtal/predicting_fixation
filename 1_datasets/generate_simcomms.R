a=Sys.time()

## =============================================================================
## FUNCTIONS
## =============================================================================
my_transpose <- function(df) {
  library(data.table, quietly = 1)
  t_df = data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  return(t_df)
}


simulate_uniform_community <- function(num_species=10,
                                       comm_size=10**4,
                                       name="otu") {
  table(
    sample(x = paste0(name, 1:num_species),
           replace = TRUE,
           size = comm_size)
  )
}

simulate_log_normal_community <- function(num_species=10,
                                          comm_size=10**4,
                                          name="otu",
                                          meanlog = 0,
                                          sdlog = 1) {
  # https://www.nature.com/scitable/knowledge/library/explaining-general-patterns-in-species-abundance-and-23162842/
  table(
    sample(x = paste0(name, 1:num_species),
           replace = TRUE,
           size = comm_size,
           prob = stats::rlnorm(n = num_species, meanlog = 0, sdlog = 1))
  )
}

simulate_log_normal_community2 <- function (num_species = 10,
                                            comm_size = 10000,
                                            name="otu",
                                            meanlog = -.5,
                                            sdlog = 1.4,
                                            plot = FALSE,
                                            ...) {
  # (This function is more consistent, meaning the shape of the distribution is
  # always similar. On the contrary, simulate_log_normal_community() would
  # sometimes return distributions that are not that skewed)
  # num_species ==> number of species
  # comm_size   ==> community size
  x_dlnorm <- seq(1, num_species, by = 1)   # Specify x-values for dlnorm function; this is 100 species

  y_dlnorm <- dlnorm(x_dlnorm, meanlog = meanlog, sdlog = sdlog)   # Apply dlnorm function
  y_dlnorm <- y_dlnorm * comm_size / sum(y_dlnorm) # adjust to final abun

  if (plot) {
    barplot(y_dlnorm, ...)
  }
  return(y_dlnorm)
}

# Take a look:
# ------------------------------------------------------------------------------
# # pdf("test_lognormal.pdf")
# for (ns in c(10, 100, 1000)) {
#   for (cs in c(10000, 1000000)) {
#     simulate_log_normal_community2(ns, cs, 0, 1, plot = TRUE, col="cyan")
#       barplot(sort(simulate_log_normal_community(ns, cs, "otu", 0, 1), decreasing = T), col="coral")
#   }
# }
# # dev.off()

simulate_skewed_community <- function(num_species=10,
                                      comm_size=10**4,
                                      name="otu",
                                      num_winners=2,
                                      fold=2) {
  # num_winners --> number of exceptionally abundant OTUs
  # fold ---------> how exceptionally abundant they are
  table(
    sample(x = paste0(name, 1:num_species),
           replace = TRUE,
           size = comm_size,
           prob = c(rep(2, num_winners), rep(1, num_species - num_winners)))
  )
}


## =============================================================================
## Parameters
## =============================================================================
# Experimental design
total_samples_u  <- 30 # going to simulate N original communities for each
# variant / abundance distribution
total_samples_ln <- 30

# 10, 100, 1000 especies
num_species <- c(10, 100, 1000)
# tamaños: 10⁴, 10⁶, 10⁵, 10³
comm_size <- c(10**4)
# skewed
# total_samples_s  <- 10
# num_winners <- c(1, 2, 3, 5, 8, 10)
# fold <- c(1.5, 2, 3, 5, 10)

outputdir <- "~/AAA/2023-02-01_initial_community_simulation/simcomms_2"
if (!file.exists(outputdir)) {system(paste("mkdir -p", outputdir))}

## =============================================================================
## Uniform
## =============================================================================
simcomms <- c()
for (ns in num_species) {
  for (cs in comm_size) {
    s <- 0
    while (s < total_samples_u) {
      s <- s + 1
      simcomms <- rbind(simcomms, # NOTE : column names might not match for R < 3.5.2
                        # not super important
                        simulate_uniform_community(num_species = ns,
                                                   comm_size = cs,
                                                   name = "otu"))
      message(paste0("Created uniform simcomm #", s))
    }
    write.table(my_transpose(as.data.frame(simcomms)),
                sep = "\t",
                col.names = NA,
                paste0(outputdir, "/uniform_", ns,"sp_", "size", cs, ".tsv")
    )
    simcomms <- c()
  }
}


## =============================================================================
## Skewed
## =============================================================================
# simcomms <- c()
# for (ns in num_species) {
#   for (cs in comm_size) {
#     for (nw in num_winners) {
#       if (nw <= (comm_size/2)) { # no point on more than half the members being winners
#         for (f in fold) {
#           s <- 0
#           while (s < total_samples_s) {
#             s <- s + 1
#             simcomms <- rbind(simcomms, # NOTE : column names might not match for R < 3.5.2
#                               # not super important
#                               simulate_skewed_community(num_species = ns,
#                                                         comm_size = cs,
#                                                         name = "otu",
#                                                         num_winners = nw,
#                                                         fold = f))
#             message(paste0("Created skewed simcomm #", s))
#           }
#           write.table(my_transpose(as.data.frame(simcomms)),
#                       sep = "\t",
#                       col.names = NA,
#                       paste0(outputdir, "/skewed_", ns,"sp_", "size", cs, "_", nw, "winners_", f, "fold", ".tsv")
#           )
#           simcomms <- c()
#         }
#       }
#     }
#   }
# }

## =============================================================================
## Log-normal
## =============================================================================
simcomms <- c()
for (ns in num_species) {
  for (cs in comm_size) {
    s <- 0
    while (s < total_samples_u) {
      s <- s + 1
      simcomms <- rbind(simcomms, # NOTE : column names might not match for R < 3.5.2
                        # not super important
                        simulate_log_normal_community2(num_species = 100,
                                                       comm_size = 1000000,
                                                       name = "otu",
                                                       meanlog = 0,
                                                       sdlog = 1,
                                                       plot=TRUE)
                        )
      message(paste0("Created log-normal simcomm #", s))
    }
    write.table(my_transpose(as.data.frame(simcomms)),
                sep = "\t",
                col.names = NA,
                paste0(outputdir, "/lognorm_", ns,"sp_", "size", cs, ".tsv")
    )
    simcomms <- c()
  }
}

print(Sys.time()-a)
