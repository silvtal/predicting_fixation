library(party)
library(dplyr)
library(pscl)
library(MASS)

# visualization packages
library(flexplot) #devtools::install_github("dustinfife/flexplot")
    ## yaxis ~ xaxis + colors | column_panels + row_panels <- this is the format
    ### "|" symbol:
    ### Instead of assuming a single relationship between success and shannon
    ### across all levels of dilfactor and size, the model will estimate separate
    ### relationships for each combination. This allows us to account for potential
    ### interactions or differences in the effect of shannon on success based
    ### on the levels of dilfactor and size.
library(ggplot2)
library(DHARMa)
library(sjPlot)

for (threshold in c(0.5, 0.9)) {
  
  # OPTIONS -----------------------------------------------------------------
  out_folder = "../figures/RF_success"; if (!file.exists(out_folder)) {system(paste("mkdir -p", out_folder))}
  my_file = paste0("../1_datasets/simcomms/processed_data_simcomms_", threshold, "_full_jun")
  prefix = paste0(threshold*100, "_")
  maxdilfactor <- 0.01
  mindilfactor <- 0.00025
  my_family <- Gamma(link = "log")
  
  # Load the input data
  csv <- read.csv(my_file)
  csv <- csv[!(colnames(csv) %in%  c("final_size", "filename", "sample"))]
  
  csv["distrib"]   <- ifelse(csv$distrib == "uniform", 1, 0)
  csv["dilfactor"] <- log(csv$dilfactor)
  print(paste("Muestras procesadas:", nrow(csv)))
  
  csv <- csv %>%
    na.omit()
  print(paste("Eliminando los NA:", nrow(csv)))
  
  csv <- csv[csv$dilfactor < maxdilfactor & csv$dilfactor > mindilfactor, ]
  print(paste("Filtrando por dilution factor (>", mindilfactor, "; <", maxdilfactor, "):", nrow(csv)))
  
  
  
  # RANDOM FOREST (variable selection) --------------------------------------
  # is not really necessary to transform / even out the data but I'm
  # using a boxcox transformation. That's because if we even out our data, the
  # decision trees will be as accurate for values from one extreme as for
  # the rest of the values. See: https://stats.stackexchange.com/a/448153
  png(paste0(out_folder, "/", prefix, "before_boxcox.png"),
      width = 1000, height = 600)
  flexplot(success~1, data=csv)
  dev.off()
  
  bc <- boxcox(success~1, data=csv)
  # csv$success_transformed <- csv$success^(bc$x[which(bc$y==max(bc$y))])
  csv$success <- csv$success^(bc$x[which(bc$y==max(bc$y))])
  png(paste0(out_folder, "/", prefix, "after_boxcox.png"),
      width = 1000, height = 600)
  flexplot(success~1, data=csv)
  dev.off()
  
  # make a model with all variables
  rf_model <- cforest(success~.,
                      controls = cforest_control(ntree = 1000),
                      data=csv)
  
  rf_estimates <- estimates(rf_model)
  
  print(rf_estimates$rsq)
  # R² is high
  
  print(rf_estimates$oob)
  # OOB performance: abs(predicted - actual)
  # there's one calculation for each sample, and then the results get sorted
  # 50%: on average, the OOB performance is this, meaning the model only fails
  # by <this value> in each prediction
  # 100%: the maximum deviation that happens, between what's predicted and reality
  # > not that good
  
  print(rf_estimates$importance)
  # Variable importance: uses Mean Decrease Accuracy (MDA). This method
  # assesses the importance of a variable by evaluating the decrease in model
  # accuracy when the variable is randomly permuted
  
  # The ones I cannot remove are size and dilfactor, but Shannon/evenness/Gini look
  # relatively important, as does richness. If I understand the method correctly,
  # however, this does not account for redundancy.
  
  png(paste0(out_folder, "/", prefix, "feature_importance.png"),
      width = 1000, height = 600)
  barplot(rf_estimates$importance,
          col = c("#FFB3D1", "#FFC3A0", "#FFD8BF", "#FFE8A0", "#FFF6B3",
                  "#D1FFB3", "#A0FFC3", "#BFFFD8", "#A0FFE8", "#B3FFF6")[1:length(rf_estimates$importance)]
  )
  dev.off()
  
  # choosing a diversity metric / should I exclude richness? -----------------
  ## Made 3 models, each with a different diversity metric, and asked ChatGPT
  ## about it (reached the same conclusion as me). In all three cases, the
  ## diversity metric is more important on its own than with the others.
  rf_small1 <- cforest(success~size + dilfactor + evenness,
                       data=csv)
  
  rf_small2 <- cforest(success~size + dilfactor + shannon,
                       data=csv)
  
  rf_small3 <- cforest(success~size + dilfactor + gini,
                       data=csv)
  
  
  rich_rf_small1 <- cforest(success~size + dilfactor + evenness + richness,
                            data=csv)
  
  rich_rf_small2 <- cforest(success~size + dilfactor + shannon + richness,
                            data=csv)
  
  rich_rf_small3 <- cforest(success~size + dilfactor + gini + richness,
                            data=csv)
  
  ## Comparison:
  diversities <- c()
  all_estimates <- data_frame()
  all_r2 <- data_frame()
  all_oob <- data_frame()
  
  # even, shannon, gini, even, shannon, gini
  models <- c(rf_small1, rf_small2, rf_small3, rich_rf_small1, rich_rf_small2, rich_rf_small3)
  for (m in 1:length(models)) {
    e <- estimates(models[[m]])
    all_estimates <- bind_rows(all_estimates, e$importance)
    all_r2 <- bind_rows(all_r2, e$rsq)
    all_oob <- bind_rows(all_oob, e$oob)
    diversities[m] = e[["importance"]][names(e[["importance"]])[names(e[["importance"]]) %in% c("gini", "evenness", "shannon")]]
    png(paste0(out_folder, "/", prefix, "divers_check", m, "_feature_importance.png"),
        width = 1000, height = 600)
    barplot(e$importance,
            col = c("#FFB3D1", "#FFC3A0", "#FFD8BF", "#FFE8A0", "#FFF6B3",
                    "#D1FFB3", "#A0FFC3", "#BFFFD8", "#A0FFE8", "#B3FFF6")[1:length(rf_estimates$importance)]
    )
    dev.off()
  }
  
  ## With richness, the Shannon diversity metric loses importance. however, the others
  ## don't lose that much
  # size dilfactor evenness shannon  gini richness
  # <dbl>     <dbl>    <dbl>       <dbl> <dbl>    <dbl>
  #   1  334.      256.     69.3        NA    NA       NA  
  # 2  319.      259.     NA          84.2  NA       NA  
  # 3  321.      251.     NA          NA    69.9     NA  
  # 4  337.      263.     68.2        NA    NA       60.0
  # 5  318.      259.     NA          64.1  NA       66.5
  # 6  337.      263.     NA          NA    71.5     59.1
  
  ## R²: richness improves the model, but without richness the best model is the
  ## Shannon one:
  # success
  # 1 0.9189764
  # 2 0.9398986
  # 3 0.9220420
  # 4 0.9423304
  # 5 0.9441217
  # 6 0.9415435
  
  
  
  # improving the model: account for interactions ---------------------------
  
  ## First, some visualization
  ## the diagonal red lines here are more inclined when that feature is relevant
  visualize(rf_small2)
  
  predictions <- compare.fits(success~dilfactor + shannon + size,
                              data = csv,
                              model1 = rf_small2,
                              return.preds = T)
  predictionsR <- compare.fits(success~dilfactor + shannon + size + richness,
                               data = csv,
                               model1 = rich_rf_small2,
                               return.preds = T)
  head(predictions)
  # 1. shannon depends on dilfactor and size
  rf_slash1 <- cforest(success~ shannon | I(log(dilfactor)) + size,
                       data = csv)
  visualize(rf_slash1)
  
  # flexplot(formula = success~ shannon | dilfactor + size + richness,
  #          data = csv,
  #          method = "glm",
  #          prediction = predictionsR,
  #          plot.type = "line", ghost.reference = "size"
  #          )
  flexplot(formula = success~  log(dilfactor) + shannon | size + richness,
           data = csv,
           method = "glm",
           prediction = predictionsR,
           plot.type = "line", ghost.reference = "size"
  )
  
  flexplot(formula = success~  log(dilfactor) + shannon +size,
           data = csv,
           method = "glm",
           prediction = predictionsR,
           plot.type = "line", ghost.reference = "size"
  )
  
}

# just save the image -----------------------------------------------------
save.image(paste0(out_folder, "/RF_success.RData"))
