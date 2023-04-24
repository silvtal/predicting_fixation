## Read filenames and create tibble for all the simul_data
filenames <- list.files(path = simuls_folder, 
                        recursive = TRUE,
                        pattern=glob2rx(paste0("*csv")),
                        full.names=FALSE)

simul_data  <- mapply(filenames, FUN=function(name){
  split <- str_split(name, "_")[[1]]
  return(c("core" = split[5],
           "sample" = split[6],
           "transfer" = str_split(split[8], ".csv")[[1]][1],
           "dilfactor" = str_split(split[1], "/")[[1]][1],
           "filename" = name)
  )
}) %>% t %>% as_tibble()

# ## Load files into a tibble
# simul_data$absabund <- lapply(simul_data$filename,
#                             FUN=function(n) {
#                               abuntable <- data.table::fread(paste0(simuls_folder, "/", n), header = T, drop = 1) %>%
#                                 as_tibble()
#                               # if some row is filled with "NA", that means some OTU reached fixation in the
#                               # last transfer of that trajectory
#                               return(abuntable)
#                             }
# )
# 
# simul_data$initabund <- lapply(simul_data$absabund,
#                              FUN=function(n) {
#                                return(sum(n[1,])) ## pillo una simul cualquiera, la 1; nos interesa solo a tiempo 0 y en ese caso todas las filas son iguales
#                              }
# )
# 
# ## Normalize to relative abundance, save into the tibble
# simul_data$table <- lapply(simul_data$filename,
#                          FUN=function(n) {
#                            abuntable <- data.table::fread(paste0(simuls_folder, "/", n), header = T, drop = 1) %>%
#                              as_tibble()
#                            # if some row is filled with "NA", that means some OTU reached fixation in the
#                            # last transfer of that trajectory
#                            # TODO with this trick we wont know *which one* is fixed !
#                            abuntable <- (abuntable/rowSums(abuntable)) %>% replace(is.na(.), fixation_threshold)
#                            return(abuntable)
#                          }
# )
# 
# 
# ## $max_abunds
#   simul_data$max_abuns <- lapply(simul_data$table, FUN= function(t){
#     apply(t, MARGIN=1, max) %>% # if we have 500 simuls/trajectories,
#       sort(decreasing = T)      # we will have 500 max abundances
#   })                            
#   
# ## $perc_95
# simul_data$perc_95 <- lapply(simul_data$max_abuns, FUN=function(ma) {
#   ma <- head(ma, trunc(length(ma)*0.95)) # aplico 95% !
#   sum(ma>=fixation_threshold)/length(ma)})


simul_data$transfer <- as.numeric(simul_data$transfer)
simul_data$dilfactor <- as.numeric(simul_data$dilfactor)

simul_data <- arrange(simul_data, transfer)

total_transfers <- max(simul_data$transfer)
