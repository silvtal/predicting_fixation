library(dplyr)
library(flexplot)

out_folder = "../figures/0__feature_effects"; if (!file.exists(out_folder)) {system(paste("mkdir -p", out_folder))}
for (threshold in c(0.5, 0.9)) {
  # OPTIONS -----------------------------------------------------------------
  my_file = paste0("../1_datasets/simulation_results/processed_data_simcomms_", threshold, "_full_jun")
  prefix = paste0(threshold*100, "_")
  if (threshold == 0.9) {
    maxdilfactor <- 0.1
    mindilfactor <- 0
  } else {
    maxdilfactor <- 0.25
    mindilfactor <- 0
  }
  my_family <- Gamma(link = "log")
  
  # Load the input data
  csv <- read.csv(my_file)
  csv <- csv[!(colnames(csv) %in%  c("final_size", "filename", "sample"))]
  
  csv["distrib"] <- ifelse(csv$distrib == "uniform", 1, 0)
  csv["logdf"]   <- log(csv$dilfactor)
  print(paste("Muestras procesadas:", nrow(csv)))
  
  csv <- csv %>%
    na.omit()
  print(paste("Eliminando los NA:", nrow(csv)))
  
  csv <- csv[csv$dilfactor < maxdilfactor & csv$dilfactor > mindilfactor, ]
  print(paste("Filtrando por dilution factor (>", mindilfactor, "; <", maxdilfactor, "):", nrow(csv)))
 
  for (var in colnames(csv)[colnames(csv) != "success"]) {
    png(paste0(out_folder, "/", prefix, var, "_effect.png"),
        width = 1000, height = 600)
    csv$var <- csv[[var]]
    f <- flexplot(success~var, data=csv)
    plot(f)
    csv$var <- NULL
    dev.off()
  }
}
