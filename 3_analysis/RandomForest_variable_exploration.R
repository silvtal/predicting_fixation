## Script viejo en el que uso GLM para success, pero al final he decidido no
## hacerlo porque la distribución es rara
## -----------------------------------------------------------------------------

### https://youtu.be/QctVIuLTKtY?list=PL8F480DgtpW9W-PEX0f2gHl8SnQ7PtKBv
# 1)  Usar RF para sacar las variables más relevantes
# 2)  Tomar esos predictores e incluirlos en mi GLM (interacciones, términos no lineales...)

library(party)
library(dplyr)
library(pscl)


# visualization packages
library(flexplot) #devtools::install_github("dustinfife/flexplot")
## yaxis ~ xaxis + colors | column_panels + row_panels <- this is the format
### "|" symbol:
### Instead of assuming a single relationship between success and raw_shannon
### across all levels of dilfactor and size, the model will estimate separate
### relationships for each combination. This allows us to account for potential
### interactions or differences in the effect of raw_shannon on success based
### on the levels of dilfactor and size.
library(ggplot2)
library(DHARMa)
library(sjPlot)


# OPTIONS -----------------------------------------------------------------
# input
threshold <- 50
my_file = paste0("/home/silvia/repos/predicting_fixation/1_datasets/simcomms/processed_data_simcomms_", threshold/10, "_full_jun16")

# output
out_folder = paste0("/home/silvia/repos/predicting_fixation/3_analysis/RF_", threshold, "%/")
if (!file.exists(out_folder)) {system(paste0("mkdir -p ", out_folder))}

# filtrar por dilfactor
maxdilfactor <- 0.01
mindilfactor <- 0.00025

# model
my_family <- Gamma(link = "log")


# load data ---------------------------------------------------------------
csv <- read.csv(my_file) %>%
  select(-reached_fixation_at, -final_size, -initial_size, -transfer, -filt_shannon, -filt_even, -filename, -sample) 
csv["distrib"] <- ifelse(csv$distrib == "uniform", 1, 0)
csv["log(dilfactor)"] <- log(csv$dilfactor)
print(paste("Muestras procesadas:", nrow(csv)))
csv <- csv %>%
  na.omit()
print(paste("Eliminando los NA:", nrow(csv)))
csv <- csv[csv$dilfactor < maxdilfactor & csv$dilfactor > mindilfactor, ]
print(paste("Filtrando por dilution factor (>", mindilfactor, "; <", maxdilfactor, "):", nrow(csv)))


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


# choosing a diversity metric / should I exclude richness? -----------------
## Made 3 models, each with a different diversity metric, and asked ChatGPT
## about it (reached the same conclusion as me). In all three cases, the
## diversity metric is more important on its own than with the others.
rf_small1 <- cforest(success~size + dilfactor + raw_even,
                     data=csv)

rf_small2 <- cforest(success~size + dilfactor + raw_shannon,
                     data=csv)

rf_small3 <- cforest(success~size + dilfactor + gini,
                     data=csv)


rich_rf_small1 <- cforest(success~size + dilfactor + raw_even + richness,
                          data=csv)

rich_rf_small2 <- cforest(success~size + dilfactor + raw_shannon + richness,
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
  diversities[m] = e[["importance"]][names(e[["importance"]])[names(e[["importance"]]) %in% c("gini", "raw_even", "raw_shannon")]]
}

## With richness, the Shannon diversity metric loses importance. however, the others
## don't lose that much
# size dilfactor raw_even raw_shannon  gini richness
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

predictions <- compare.fits(success~dilfactor + raw_shannon + size,
                            data = csv,
                            model1 = rf_small2,
                            return.preds = T)
predictionsR <- compare.fits(success~dilfactor + raw_shannon + size + richness,
                             data = csv,
                             model1 = rich_rf_small2,
                             return.preds = T)
head(predictions)
# 1. shannon depends on dilfactor and size
rf_slash1 <- cforest(success~ raw_shannon | I(log(dilfactor)) + size,
                     data = csv)
visualize(rf_slash1)

# flexplot(formula = success~ raw_shannon | dilfactor + size + richness,
#          data = csv,
#          method = "glm",
#          prediction = predictionsR,
#          plot.type = "line", ghost.reference = "size"
#          )
flexplot(formula = success~  log(dilfactor) + raw_shannon | size + richness,
         data = csv,
         method = "glm",
         prediction = predictionsR,
         plot.type = "line", ghost.reference = "size"
)

flexplot(formula = success~  log(dilfactor) + raw_shannon +size,
         data = csv,
         method = "glm",
         prediction = predictionsR,
         plot.type = "line", ghost.reference = "size"
)




# GLM: test different models ----------------------------------------------
## Now we'll take the selected variables and combine them into different GLMs
# Residual deviance: The residual deviance represents the deviance of the fitted
# model after accounting for the predictors. It measures the lack of fit of the
# model after considering the predictors. In this case, the residual deviance is
# 180.82, indicating that the model has significantly reduced the deviance
# compared to the null model, suggesting a good fit.
#
# intercept: the success when everything is 0.

# BETTER NOT TO -- Normalize the variables in X for GLM
# X <- as.data.frame(sapply(X, function(x) (x - min(x)) / (max(x) - min(x))),

test_glm <- function(variable_names) {
  # separate_with -->  character that separates the interacting features from
  #                    the others in the formula
  glm_models <- list()
  formulas <- list()
  for (i in 1:length(variable_names)) {
    combos <- combn(variable_names, i, simplify = FALSE)
    for (j in 1:length(combos)) {
      if (i >= 2) {
        # Create the formula string
        formula_str <- paste("success ~", paste(combos[[j]], collapse = "+"))
        
        # Store the formula in the list
        formulas[[length(formulas) + 1]] <- formula_str
        
        # Create the formula object
        formula_obj <- formula(formula_str)
        
        # Fit the GLM
        glm_model <- glm(formula_obj, data = csv, family = my_family)
        
        # Store the model in the list
        glm_models[[length(glm_models) + 1]] <- glm_model
        
        # Generate formulas with interactions
        interactions <- combn(combos[[j]], 2, simplify = FALSE)
        for (k in 1:length(interactions)) {
          # Create the formula string with interactions
          formula_str_int <- paste("success ~", paste(interactions[[k]], collapse = ":"), "+", paste(combos[[j]], collapse = "+"))
          
          # Store the formula in the list
          formulas[[length(formulas) + 1]] <- formula_str_int
          
          # Create the formula object with interactions
          formula_obj_int <- formula(formula_str_int)
          
          # Create the GLM model with interactions
          glm_model_int <- glm(formula_obj_int, data = csv, family = my_family)
          
          # Store the model in the list
          glm_models[[length(glm_models) + 1]] <- glm_model_int
        }
      }
    }
  }
  
  # Access the formulas and models for further analysis
  results <- data_frame()
  for (f in 1:length(formulas)) {
    model <- glm_models[[f]]
    
    summary(model)
    coefficients <- coef(model)
    
    # Regression formula
    regression_formula <- "y = "
    for (v in seq_along(variable_names)) {
      coefficient <- coefficients[v]
      variable_name <- names(coefficients)[v]
      regression_formula <- paste0(regression_formula, sprintf("(%0.4f * %s) + ", coefficient, variable_name))
    }
    # Lower deviance values indicate a better fit
    results <- bind_rows(results, data_frame(formulas[[f]], model$deviance, substr(regression_formula, 1, nchar(regression_formula) - 3)))
  }
  r <- results[order(results$`model$deviance`), ]
  View(r)
  return(r)
}

test1 <- test_glm(c("size", "dilfactor", "raw_shannon", "richness", "log(dilfactor)"))
### Some models are terrible and some others are good:

### when not including size, all models are very bad (5000 deviance)
###       --> success ~ dilfactor:raw_shannon + dilfactor+raw_shannon+richness

### with size but without dilfactor, models are still very bad (1000 deviance)

### models that use size and dilfactor with deviance of ~300 don't use "log(dilfactor)". See:
# > print(r[(r$`model$deviance` >300) & (r$`model$deviance` <500), c("formulas[[f]]", "model$deviance")])
# # A tibble: 17 × 2
# `formulas[[f]]`                                                       `model$deviance`
# <chr>                                                                            <dbl>
#   1 success ~ size:dilfactor + size+dilfactor+raw_shannon+richness                    328.
# 2 success ~ size:dilfactor + size+dilfactor+raw_shannon                             328.
# 3 success ~ size:dilfactor + size+dilfactor+richness                                338.
# 4 success ~ size:dilfactor + size+dilfactor                                         343.
# 5 success ~ dilfactor:raw_shannon + size+dilfactor+raw_shannon+richness             357.
# 6 success ~ dilfactor:raw_shannon + size+dilfactor+raw_shannon                      358.
# 7 success ~ raw_shannon:richness + size+dilfactor+raw_shannon+richness              361.
# 8 success ~ size:raw_shannon + size+dilfactor+raw_shannon+richness                  368.
# 9 success ~ size:raw_shannon + size+dilfactor+raw_shannon                           368.
# 10 success ~ dilfactor:richness + size+dilfactor+raw_shannon+richness                368.
# 11 success ~ size:richness + size+dilfactor+raw_shannon+richness                     372.
# 12 success ~ size+dilfactor+raw_shannon+richness                                     374.
# 13 success ~ size+dilfactor+raw_shannon                                              375.
# 14 success ~ size:richness + size+dilfactor+richness                                 381.
# 15 success ~ dilfactor:richness + size+dilfactor+richness                            381.
# 16 success ~ size+dilfactor+richness                                                 382.
# 17 success ~ size+dilfactor                                                          384.
### finally, the top models (<300, there are like 20) include  size and log(dilfactor) but not necessarily the others.
### richness is not necessary: models are exactly as good or bad without it
# 1 success ~ size:log(dilfactor) + size+dilfactor+raw_shannon+richness+log(dilfactor)                   181.
# 2 success ~ size:log(dilfactor) + size+dilfactor+raw_shannon+log(dilfactor)                            181.

# 3 success ~ size:log(dilfactor) + size+raw_shannon+richness+log(dilfactor)                             182.
# 4 success ~ size:log(dilfactor) + size+raw_shannon+log(dilfactor)                                      182.
### and apparently dilfactor (independent) is not useful either, compare 1-2 above with 3-4.

### now I'm going to check a couple more models
# original : "success ~ size:log(dilfactor) + size+raw_shannon+log(dilfactor)"

summary(glm("success ~ size:log(dilfactor) + raw_shannon+log(dilfactor)", data = csv, family = my_family))
# removing the independent "size" variable increases the deviance to 477.41

summary(glm("success ~ size:log(dilfactor) + size:raw_shannon + log(dilfactor)", data = csv, family = my_family))
# 484.53

summary(glm("success ~ size:log(dilfactor) + raw_shannon + size + raw_shannon:log(dilfactor)", data = csv, family = my_family))
# 314.29

summary(glm("success ~ size:log(dilfactor) + size + raw_shannon:log(dilfactor)", data = csv, family = my_family))
# removing raw_shannon on its own is hurtful to the model as well. 700

summary(glm("success ~ size:log(dilfactor):raw_shannon  + size + log(dilfactor)", data = csv, family = my_family))
# 220.59

### Change diversity metric again
myformula <- "success ~ size:log(dilfactor) +  raw_even  + size  + log(dilfactor)"
mymodel <- glm(myformula, data = csv, family = my_family)
summary(mymodel)$deviance
# 172.27 <- choose this model

summary(glm("success ~ size:log(dilfactor) +  gini  + size  + log(dilfactor)", data = csv, family = my_family))
# 190

summary(glm("success ~ size:log(dilfactor) +  raw_even*size  + log(dilfactor)", data = csv, family = my_family))
# 168.29 (mucho complicar quiza)

summary(glm("success ~ size * log(dilfactor) * raw_even", data = csv, family = my_family))
# 167.35 (mucho complicar...)
# (Intercept)                   2.223e+00  7.010e-01   3.171 0.001545 ** 
# size                          6.590e-06  1.104e-06   5.968 2.88e-09 ***
# log(dilfactor)                2.865e-01  1.091e-01   2.627 0.008696 ** 
# raw_even                      4.173e+00  7.468e-01   5.588 2.64e-08 ***
# size:log(dilfactor)           6.135e-07  1.684e-07   3.643 0.000277 ***
# size:raw_even                -1.123e-06  1.184e-06  -0.949 0.342923    
# log(dilfactor):raw_even       4.317e-01  1.160e-01   3.720 0.000205 ***
# size:log(dilfactor):raw_even -4.071e-07  1.801e-07  -2.261 0.023887 * 
# just choose mymodel





# assess model --------------------------------------------------------------
myformula <- "success ~ log(dilfactor):raw_even + raw_even  + size  + log(dilfactor)"
mymodel <- glm(myformula, data = csv, family = my_family)
summary(mymodel)
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -1.19969  -0.14187   0.02204   0.16801   0.52744  
## > 1Q and 3Q abs values should be similar.
## > median should be close to 0
## min and max should be under 3 (good)
## and also similar to each other (not good). HERE IT FAILS. wrong distribution?
tab_model(mymodel)



# Get the t-values from the model summary
t_values <- summary(mymodel)$coefficients[, "t value"]
# Create a bar plot of the t-values
barplot(t_values, 
        names.arg = names(t_values), 
        xlab = "Variables", 
        ylab = "t-values",
        main = "T-values of Coefficients",
        border = "black",
        col = t_values
)




# Now we can run the anova() function on the model to analyze the table of deviance
# of course, anova can also be used for comparing two models
anova(mymodel, test="Chisq")
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
pscl::pR2(mymodel)




## we can see how the levine test fails here, different variance between feature values
simulatedOutput <- simulateResiduals(mymodel, )
png(paste0(out_folder, "/dharma_simulateResiduals.png"), width = 1000, height = 600)
plot(simulatedOutput)
dev.off()
## from https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
# If you see this pattern, note that a common reason for underdispersion is
# overfitting, i.e. your model is too complex.
## mymodel2 <- glm("success ~ raw_even  + size  + log(dilfactor)", data = csv, family = my_family)
##  --> dharma's the same actually, just more deviance




testDispersion(mymodel)
testDispersion(simulatedOutput)
# In the normal DHARMa residual plots, zero-inflation will look pretty much like overdispersion
# The reason is that the model will usually try to find a compromise between
# the zeros, and the other values, which will lead to excess variance in the residuals.
testZeroInflation(simulatedOutput)







# we should make two models -----------------------------------------------
# Finally, let's visualize a model with richness and shannon, but algo with logged dilfactor
csv$logdilfactor <- log(csv$dilfactor)
rich_rf2_small_log <- cforest(success~size + logdilfactor + raw_shannon + richness,
                              data=csv)

predictionsL <- compare.fits(formula = success~size + logdilfactor + raw_shannon + richness,
                             data = csv,
                             model1 = rich_rf2_small_log,
                             return.preds = T)

png(paste0(out_folder, "/success_by_dilfactorxsize.png"), width = 1200, height = 500)
flexplot(formula = success~ raw_shannon | logdilfactor + size,
         data = csv, #[csv$distrib==1000000",],
         method = "glm",
         prediction = predictionsL,
         plot.type = "line", spread = "quartiles",
         alpha = .4, bins = 7
)
dev.off()

## This one shows us that maybe we should separate the data by size,
## maybe do two models

mymodel1 <- glm("success ~ raw_even  + log(dilfactor) + richness",
                data = csv[csv$size == 10000, ],
                family = my_family)
summary(mymodel1)

# evenness is better for this one...
# richness is not that much of an improvement for this one.
mymodel2 <- glm("success ~ raw_even  + log(dilfactor) + richness",
                data = csv[csv$size == 1000000, ],
                family = my_family)
summary(mymodel2)


# barplots ----------------------------------------------------------------
# barplots from "GLM short.ipynb" are easier to understand, with and without normalization

# feature importance (RF !!!!)
rf_estimates <- estimates(cforest(success ~ size:log(dilfactor) +  raw_even  + size  + log(dilfactor),
                                  controls = cforest_control(ntree = 1000),
                                  data=csv))
feature_importance <- rf_estimates$importance

# bar plot for feature importance
feature_df <- data.frame(Feature = names(feature_importance),
                         Importance = feature_importance)

feature_plot <- ggplot(data = feature_df, aes(x = Feature, y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Feature", y = "Importance") +
  ggtitle("Feature Importance")

# display
png(paste0(out_folder, "/success_FI.png"), width = 600, height = 400)
print(feature_plot)
dev.off()

# coefs themselves (GLM)
coefficients <- coef(mymodel)

# bar plot for coefficients
coef_df <- data.frame(Coefficient = names(coefficients),
                      Value = coefficients)

coef_plot <- ggplot(data = coef_df, aes(x = Coefficient, y = Value)) +
  geom_bar(stat = "identity", fill = "lightpink") +
  labs(x = "Coefficient", y = "Value") +
  ggtitle("Coefficients") + coord_flip()

# display
png(paste0(out_folder, "/success_coefs.png"), width = 600, height = 400)
print(coef_plot)
dev.off()
