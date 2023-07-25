
# libraries ---------------------------------------------------------------
library(party)
library(flexplot) #devtools::install_github("dustinfife/flexplot")
library(dplyr)

# visualization packages
library(ggplot2)
library(DHARMa)
library(sjPlot)

results <- list()
out_folder = "../figures/2__GLM_failure/"; if (!file.exists(out_folder)) {system(paste("mkdir -p", out_folder))}
for (threshold in c(0.5, 0.9)) {
  # OPTIONS -----------------------------------------------------------------
  my_file = paste0("../1_datasets/simulation_results/processed_data_simcomms_", threshold, "_full_jul")
  prefix = paste0(threshold*100, "_")

  # Load the input data
  fcsv <- read.csv(my_file)
  fcsv <- fcsv[!(colnames(fcsv) %in%  c("final_size", "filename", "sample"))]
  fcsv["distrib"]   <- ifelse(fcsv$distrib == "uniform", 1, 0)
  fcsv["dilfactor"] <- log(fcsv$dilfactor)
  print(paste("Muestras procesadas:", nrow(fcsv)))
  
  # fcsv <- fcsv[fcsv$dilfactor < maxdilfactor & fcsv$dilfactor > mindilfactor, ]
  # print(paste("Filtrando por dilution factor (>", mindilfactor, "; <", maxdilfactor, "):", nrow(fcsv)))
  
  # family has to be logistic
  my_family <- binomial(link='logit') # A binomial logistic regression attempts to
                                      # predict the probability that an observation
                                      # falls into one of two categories of a
                                      # dichotomous dependent variable based on one
                                      # or more independent variables that can be
                                      # either continuous or categorical.

  for (l in c(10, 25, 50, 100, 200, 400, 800, 1000)) {
    results[[as.character(threshold)]][[as.character(l)]] <- list()
    
    # Create "failure" variable
    # if fixation takes more than <limit> cycles to happen, that's a failure too
    fcsv["failure"] <- is.na(fcsv["success"]) | fcsv["success"] > l
    
    prefix = paste0(threshold*100, "_", l, "__")
    
    ### Change diversity metric again
    failuremodel <- glm("failure ~ size:dilfactor +  shannon  + size  + dilfactor", data = fcsv, family = my_family)
    
    # assess model --------------------------------------------------------------
    results[[as.character(threshold)]][[as.character(l)]]$summary <- summary(failuremodel)
    results[[as.character(threshold)]][[as.character(l)]]$model   <- failuremodel
    
    ## first we check the residuals
      # Deviance Residuals: 
      #   Min        1Q    Median        3Q       Max  
      # -2.70390  -0.13151  -0.00003  -0.00001   2.04691  
    ## > 1Q and 3Q abs values should be similar.
    ## > median should be close to 0
    ## min and max should be under 3 (good)
    ## and also similar to each other (good)
    
    ## we can see how the levine test goes well here, stable variance between feature values
    ## unlike:
    #### successmodel <- glm("success ~ size:dilfactor +  shannon  + size  + dilfactor", data = csv, family = my_family)
    #### simulatedOutputSuccess <- simulateResiduals(successmodel, integerResponse = F)
    simulatedOutput <- simulateResiduals(failuremodel, integerResponse = T)
    png(paste0(out_folder, "/", prefix, "dharma_simulateResiduals_FAILURE.png"), width = 1000, height = 600)
    plot(simulatedOutput)
    dev.off()
    
    
    
    tab_model(failuremodel)
    
    
    
    # Get the z-values from the model summary
    #### from summary-glm:
    #### This third column is labelled t ratio if the dispersion is estimated,
    #### and z ratio if the dispersion is known (or fixed by the family)
    #### Z-values are a proxy for the importance of each variable.
    z_values <- summary(failuremodel)$coefficients[, "z value"]
      png(paste0(out_folder, "/", prefix, "z-values_FAILURE.png"), width = 1000, height = 600)
    # Create a bar plot of the x-values
    barplot(z_values, 
            names.arg = names(z_values), 
            xlab = "Variables", 
            ylab = "z-values",
            main = "Z-values of Coefficients",
            border = "black"
    )
    dev.off()
    
    
    
    
    # Now we can run the anova() function on the model to analyze the table of deviance
    results[[as.character(threshold)]][[as.character(l)]]$anova <- anova(failuremodel, test="Chisq")
    # Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
    # NULL                                 1844     6125.7              
    # evenness             1      1.8      1843     6123.9 4.766e-07 ***
    #   size                 1   4818.8      1842     1305.2 < 2.2e-16 ***
    #   log(dilfactor)       1   1111.7      1841      193.5 < 2.2e-16 ***
    #   size:log(dilfactor)  1     21.2      1840      172.3 < 2.2e-16 ***
    # The difference between the null deviance and the residual deviance shows how our
    # model is doing against the null model (a model with only the intercept). The
    # wider this gap, the better.
    
    # While no exact equivalent to the R2 of linear regression exists, the McFadden
    # R2 index can be used to assess the model fit.
    results[[as.character(threshold)]][[as.character(l)]]$mcfadden <- pscl::pR2(failuremodel)
    
    
    
    # barplots ----------------------------------------------------------------
    
    # feature importance (RF !!!!)
    rf_estimates <- estimates(cforest(failure ~ size:dilfactor +  shannon  + size  + dilfactor,
                                      controls = cforest_control(ntree = 400),
                                      data=fcsv))
    feature_importance <- rf_estimates$importance
    results[[as.character(threshold)]][[as.character(l)]]$FI <- feature_importance
    
    # bar plot for feature importance
    feature_df <- data.frame(Feature = names(feature_importance),
                             Importance = feature_importance)
    
    feature_plot <- ggplot(data = feature_df, aes(x = Feature, y = Importance)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      labs(x = "Feature", y = "Importance") +
      ggtitle("Feature Importance")
    
    # display
      png(paste0(out_folder, "/", prefix, "failure_FI.png"), width = 600, height = 400)
    print(feature_plot)
    dev.off()
    
    # coefs themselves (GLM)
    coefficients <- coef(failuremodel)
    
    # bar plot for coefficients
    coef_df <- data.frame(Coefficient = names(coefficients),
                          Value = coefficients)
    
    coef_plot <- ggplot(data = coef_df, aes(x = Coefficient, y = Value)) +
      geom_bar(stat = "identity", fill = "lightpink") +
      labs(x = "Coefficient", y = "Value") +
      ggtitle("Coefficients") + coord_flip()
    
    # display
      png(paste0(out_folder, "/", prefix, "failure_coefs.png"), width = 600, height = 400)
    print(coef_plot)
    dev.off()
  }
}

print(results)

# just save the image -----------------------------------------------------
save.image(paste0(out_folder, "/GLM_failure.RData"))
#for (l in c(10, 25, 50, 100, 200, 400, 800, 1000)) {print(c(as.integer(l), results[[as.character(threshold)]][[as.character(l)]][["mcfadden"]]))}
# for (n in c("25", "50", "100", "200")) {print(results[["0.9"]][[n]]$summary)}