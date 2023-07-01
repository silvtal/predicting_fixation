
# libraries ---------------------------------------------------------------
library(party)
library(flexplot) #devtools::install_github("dustinfife/flexplot")
library(dplyr)

# visualization packages
library(ggplot2)
library(DHARMa)
library(sjPlot)

out_folder = "../figures/GLM_failure/"
for (threshold in c(0.5, 0.9)) {
  # OPTIONS -----------------------------------------------------------------
  
  my_file = paste0("../1_datasets/simcomms/processed_data_simcomms_", threshold, "_full_jun")

  # family has to be logistic
  my_family <- binomial(link='logit') # A binomial logistic regression attempts to
                                      # predict the probability that an observation
                                      # falls into one of two categories of a
                                      # dichotomous dependent variable based on one
                                      # or more independent variables that can be
                                      # either continuous or categorical.
  results <- list()
  for (l in c(10, 25, 50, 100, 200, 400, 800, 1000)) {
  # for (l in c(100)) {
    results[[as.character(l)]] <- list()
    
    limit <- l # if fixation takes more than <limit> cycles to happen, that's a failure too
    
    prefix = paste0(threshold*100, "_", l, "__")
    
    # Load the input data
    fcsv <- read.csv(my_file) %>%
      select(-reached_fixation_at, -final_size, -initial_size, -transfer, -filt_shannon, -filt_even, -filename, -sample) 
    fcsv["distrib"]   <- ifelse(fcsv$distrib == "uniform", 1, 0)
    fcsv["dilfactor"] <- log(fcsv$dilfactor)
    print(paste("Muestras procesadas:", nrow(fcsv)))
    
    # fcsv <- fcsv[fcsv$dilfactor < maxdilfactor & fcsv$dilfactor > mindilfactor, ]
    # print(paste("Filtrando por dilution factor (>", mindilfactor, "; <", maxdilfactor, "):", nrow(fcsv)))
    
    # Create "failure" variable
    fcsv["failure"] <- is.na(fcsv["success"]) | fcsv["success"] > limit
    
    ### Change diversity metric again
    failuremodel <- glm("failure ~ size:log(dilfactor) +  raw_even  + size  + log(dilfactor)", data = fcsv, family = my_family)
    
    # assess model --------------------------------------------------------------
    results[[as.character(l)]]$summary <- summary(failuremodel)
    results[[as.character(l)]]$model   <- failuremodel
    
    ## first we check the residuals
      # Deviance Residuals: 
      #   Min        1Q    Median        3Q       Max  
      # -2.70390  -0.13151  -0.00003  -0.00001   2.04691  
    ## > 1Q and 3Q abs values should be similar.
    ## > median should be close to 0
    ## min and max should be under 3 (good)
    ## and also similar to each other (good)
    
    ## we can see how the levine test goes well here, stable variance between feature values
    simulatedOutput <- simulateResiduals(failuremodel, integerResponse = T)
    png(paste0(out_folder, "/", prefix, "dharma_simulateResiduals_FAILURE.png"), width = 1000, height = 600)
    plot(simulatedOutput)
    dev.off()
    
    
    
    tab_model(failuremodel)
    
    
    
    # Get the z-values from the model summary
    z_values <- summary(failuremodel)$coefficients[, "z value"]
      png(paste0(out_folder, "/", prefix, "z-values_FAILURE.png"), width = 1000, height = 600)
    # Create a bar plot of the x-values
    barplot(z_values, 
            names.arg = names(z_values), 
            xlab = "Variables", 
            ylab = "z-values",
            main = "Z-values of Coefficients",
            border = "black",
            col = z_values
    )
    dev.off()
    
    
    
    
    # Now we can run the anova() function on the model to analyze the table of deviance
    anova(failuremodel, test="Chisq")
    # Df Deviance Resid. Df Resid. Dev  Pr(>Chi)    
    # NULL                                 1844     6125.7              
    # raw_even             1      1.8      1843     6123.9 4.766e-07 ***
    #   size                 1   4818.8      1842     1305.2 < 2.2e-16 ***
    #   log(dilfactor)       1   1111.7      1841      193.5 < 2.2e-16 ***
    #   size:log(dilfactor)  1     21.2      1840      172.3 < 2.2e-16 ***
    # The difference between the null deviance and the residual deviance shows how our
    # model is doing against the null model (a model with only the intercept). The
    # wider this gap, the better.
    
    # While no exact equivalent to the R2 of linear regression exists, the McFadden
    # R2 index can be used to assess the model fit.
    results[[as.character(l)]]$mcfadden <- pscl::pR2(failuremodel)
    
    
    
    # barplots ----------------------------------------------------------------
    
    # feature importance (RF !!!!)
    rf_estimates <- estimates(cforest(failure ~ size:log(dilfactor) +  raw_even  + size  + log(dilfactor),
                                      controls = cforest_control(ntree = 400),
                                      data=fcsv))
    feature_importance <- rf_estimates$importance
    results[[as.character(l)]]$FI <- feature_importance
    
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