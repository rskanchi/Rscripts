# function to perform logit regression

logitRegression <- function(data, outcome, continuos_predictors, ordinal_predictors, 
                            nominal_predictors, output_dir = NULL, ...){
  library(finalfit)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(broom)
  library(patchwork)

  if (is.null(output_dir)) output_dir <- paste("logit_analysis_output", outcome, sep = "/")
  createDir(output_dir)
  
  # outcome variable
  getPie(data = data, traits = outcome, folder = paste(output_dir, "Piecharts", "outcome", sep = "/"))
  
  # continuous predictors
  if(!is.null(continuos_predictors)){
    # histograms
    getHistogram(data = data, traits = continuos_predictors, folder = paste(output_dir, "histograms", sep = "/"))
    # AOV boxplots for the continous predictors
    getAOVBoxplot(data = data, xNames = outcome, yNames = continuos_predictors, folder = paste(output_dir, "AOVboxplots", sep = "/"))
  }
  
  # categorical predictors
  # chisq association of outcome with categorical (ordinal and nominal) predictors; barplots
  # pie charts 
  
  getPie(data = data, traits = c(ordinal_predictors, nominal_predictors), 
         folder = paste(output_dir, "Piecharts", "predictors", sep = "/"))
  
  getChiSq.association(data = data, xNames = c(ordinal_predictors, nominal_predictors),
                       yNames = outcome, output_folder = paste(output_dir, "chisq_association", sep = "/"))
  
  # logistic regression
  predictors <- c(continuos_predictors, ordinal_predictors, nominal_predictors)
  formula <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  model <- glm(formula, data = data, family = binomial)
  
  # Save model summary
  sink(file = paste(output_dir, "model_summary.txt", sep = "/"), append = TRUE)
  cat("Data Summary \n")
  print(data%>%
          summary_factorlist(dependent = outcome, explanatory = predictors, p = TRUE, na_include = TRUE, add_dependent_label = TRUE))
  cat("############### \n\n")
  cat("Logistic Regression Model Summary \n\n")
  print(summary(model))
  sink()
  
  # Coefficients as a csv/table
  coefficients_df <- tidy(model)
  write.csv(coefficients_df, paste(output_dir, "coefficients.csv", sep = "/"), row.names = FALSE)
} # end of function logitRegression


