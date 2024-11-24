# function to perform survial analysis using ML methods
# data: the outcome (status, time) and predictors
performSurvivalML <- function(data, # provide one of these two
                              #data_type = NULL, # "a tag for the data type, example "clinical", "proteome"
                              status, time, # specify column names in data that correspond to these
                              predictors, # specify column names in data that correspond to these
                              nRepeat = 10, 
                              nFold = 3, # glmnet
                              pTrain = 0.7, # 
                              pImpFreq = 0.7, # frequency of times the feature is in a model without NaN, 0
                              output_folder_tag = NULL,
                              seed = 1234567
                              
                              
){
  
  library(survival)
  library(ranger)
  library(ggplot2)
  library(dplyr)
  library(ggfortify)
  library(survminer)
  library(glmnet)
  library(caret)
  library(ggplot2)
  library(tidyr)
  
  output_folder <- paste(status, time, paste0("nRep", nRepeat), 
                         paste0("nFold", nFold), paste0("pTrain", pTrain), 
                         paste0("pImpFreq", pImpFreq), sep = "_")
  
  if (!is.null(output_folder_tag)) output_folder <- paste(output_folder_tag, output_folder, sep = "_")
  
  dir.create(output_folder, recursive = T)
  
  cat("outputs in folder:", output_folder)
  
  log_file <- file(paste0(output_folder, "/log_file.txt"), open = "wt")
  sink(log_file, append = TRUE, split = TRUE)
  sink(log_file, append = TRUE, type = "message")

  
  # data preprocessing
  data.filt <- data[, c(status, time, predictors)] #-na.indices
  colnames(data.filt) <- c("status", "time", predictors)
  data.filt[1:5,1:5]
  na.indices <- which(apply(data, 1, FUN = function(x) any(is.na(x)))) # none
  if (length(na.indices) > 0) {
    data.filt <- data.filt[-na.indices,]
    cat("Removed", length(na.indices), "rows due to NA values")
    cat("Number of rows of data", nrow(data.filt))
  }
  
  model_equation <- paste("Surv(time, status)", " ~ ", paste(predictors, collapse = "+"))
  
  # Result objects - feature importance/coefficient
  feature_importance <- list() # full matrix for all iterations and models
  ranger_importance <- list() # individual models, for ease of computing the method stats related to pImpFreq
  lasso_coef <- list()
  enet_coef <- list()
  # concordance index on test data
  ranger_cindex <- numeric()
  lasso_cindex <- numeric()
  enet_cindex <- numeric()
  
  # repeated usage
  models <- c("ranger", "lasso", "enet")
  nModels <- length(models)

  #  nRepeat loop
  for (i in 1:nRepeat) {
    cat("\n#### cross validation iteration:", i)
    set.seed(seed = seed + i)  # seed for each repeat to ensure reproducibility
    cat("\nSeed:", seed + i)
    
    # Stratify train/test split based on status for each iteration
    train_indices <- createDataPartition(data.filt$status, p = pTrain, list = FALSE, times = 1)
    cat("\ntrain data indices:", train_indices[,1])
    train_data <- data.filt[train_indices, ]
    test_data <- data.filt[-train_indices, ]
    cat("\ndimensions of train data:", dim(train_data))
    cat("\ndimensions of test data:", dim(test_data))
    
    # scale the predictor data
    train_data[, predictors] <- apply(train_data[ ,predictors], 2, scale)
    test_data[,predictors] <- apply(test_data[ ,predictors], 2, scale)
    
    # predictor and response for train and test
    x_train <- model.matrix(as.formula(paste("~", paste(predictors, collapse = "+"), "- 1")), train_data)
    x_test <- model.matrix(as.formula(paste("~", paste(predictors, collapse = "+"), "- 1")), test_data)
    y_train <- survival::Surv(time = train_data$time, event = train_data$status)
    y_test <- survival::Surv(time = test_data$time, event = test_data$status)
    
    # feature importance
    feature_importance[[i]] <- data.frame(matrix(nr=length(predictors), nc = 5)) 
    colnames(feature_importance[[i]]) <- c("iteration", "feature", models)
    feature_importance[[i]][, 1] <- rep(i, length(predictors)); feature_importance[[i]][, 2] <- predictors
      
    # 1. Ranger 
    ranger_model <- ranger(formula = as.formula(model_equation),
      data = train_data,
      num.trees = 1000, mtry = NULL,
      splitrule = "logrank",
      min.node.size = 3, # 
      importance = "impurity_corrected" #"impurity" # "permutation" gave NaNs in variable importance
    )
    ranger_importance[[i]] <- ranger_model$variable.importance
    feature_importance[[i]]$ranger <- ranger_model$variable.importance

    # Evaluate ranger model on test data; compute C index
    ranger_preds <- predict(ranger_model, data = test_data[, predictors])
    ranger_risk_scores <- ranger_preds$chf[, ncol(ranger_preds$chf)] # cumulative hazard at the last time point
    ranger_cindex[i] <- Cindex(ranger_risk_scores, y_test)

    
    # 2. LASSO Model
   lasso_model <- cv.glmnet(x = x_train, y = y_train,
                          family = "cox", alpha = 1, # alpha = 1 for lasso
                          nfolds = nFold, type.measure = "C"
    )
    
    lasso_coef[[i]] <- coef(lasso_model, s = "lambda.min")
    feature_importance[[i]]$lasso <- coef(lasso_model, s = "lambda.min")[,1]
    lasso_pred <- predict(lasso_model, newx = x_test, s = "lambda.min", type = "response")
    lasso_cindex[i] <- Cindex(lasso_pred, y_test)

    
    # 3. Elastic Net Model
    enet_model <- cv.glmnet(x = x_train, y = y_train,
      family = "cox", alpha = 0.5,
      nfolds = nFold, type.measure = "C"
    )
    enet_coef[[i]] <- coef(enet_model, s = "lambda.min")
    feature_importance[[i]]$enet <- coef(enet_model, s = "lambda.min")[,1]
    enet_pred <- predict(enet_model, newx = x_test, s = "lambda.min", type = "response")
    enet_cindex[i] <- Cindex(enet_pred, y_test)

  } # end of for in nRepeat
  
  
  # C index violin plot
  # Combine the C-index values into a single data frame
  cindex_data <- data.frame(
    Cindex = c(ranger_cindex, lasso_cindex, enet_cindex),
    Model = rep(models,
                times = c(length(ranger_cindex), length(lasso_cindex), length(enet_cindex)))
  )
  write.csv(cindex_data, paste0(output_folder, "/Cindex_nReps", nRepeat, ".csv"))
  
  # Calculate median C-index values for each model
  cindex_median_values <- cindex_data %>%
    group_by(Model) %>%
    summarize(MedianCindex = median(Cindex))
  write.csv(cindex_median_values, paste0(output_folder, "/Cindex_median_values_nReps", nRepeat, ".csv"))
  
  cindex_plot <- ggplot(cindex_data, aes(x = Model, y = Cindex)) +
    geom_violin(trim = FALSE, alpha = 0.6, width = 0.4) +  # Create the violin plot
    geom_jitter(width = 0.1, alpha = 0.4, color = "black") +  # Add points for individual C-index values
    stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "blue", color = "blue") +
    coord_cartesian(ylim = c(min(cindex_data$Cindex) - 0.2, NA)) +
    geom_text(data = cindex_median_values, aes(x = Model, y = min(cindex_data$Cindex) - 0.1, label = round(MedianCindex, 3)),
              vjust = -0.5, color = "blue", size = 3.5) + # Add text for median values
    labs(title = "Concordance index", subtitle = paste("number of iterations:", nRepeat),
         y = "C index", x = "Model") +
    theme_minimal(15) +
    theme(panel.grid.minor = element_blank())
  
  pdf(paste0(output_folder, "/Concordance_Index.pdf"), height = 5, width = 5)
  print(cindex_plot)
  dev.off()

  # list of variable importance and coefficients to matrices
  feature_importance <- do.call(rbind, feature_importance) # all iterations and models 
  write.csv(feature_importance, paste0(output_folder, "/feature_importance_wide_format_nReps", nRepeat, ".csv"))

  ranger_importance_matrix <- do.call(cbind, ranger_importance)
  lasso_coef_matrix <- do.call(cbind, lasso_coef)
  enet_coef_matrix <- do.call(cbind, enet_coef)

  # Calculate the percentage of NaN coefficients/importance across repetitions for each feature
  features_nonZero <- data.frame(matrix(data = cbind(rowSums(ranger_importance_matrix > 0) / nRepeat, 
                              rowSums(lasso_coef_matrix != 0) / nRepeat,
                              rowSums(enet_coef_matrix != 0) / nRepeat),
                              ncol = 3,
                              dimnames = list(predictors, models)
                              ))
  features_nonZero$features <- rownames(features_nonZero)
  write.csv(features_nonZero, paste0(output_folder, "/pImpFreq_all_wide_format_nReps", nRepeat, ".csv"))
  
  features_nonZero <- features_nonZero %>% pivot_longer(cols = -features, names_to = "method", values_to = "pImpFreq")
  features_pImpFreq <- features_nonZero[features_nonZero$pImpFreq >= pImpFreq, ]
  
  # Create a bar plot
  pImp_nFeatures <- data.frame(matrix(nr = nModels, nc = 1), row.names = models)
  colnames(pImp_nFeatures) <- "nFeatures"
  
  pImpFreq_plots <-   lapply(models, FUN = function(m){
    method_features <- features_pImpFreq$features[which(features_pImpFreq$method == m)]
    pImp_nFeatures[m,1] <- length(method_features)
    pImpFreq_plot_data <- feature_importance[feature_importance$feature %in% method_features, c("iteration", "feature", m)]
    colnames(pImpFreq_plot_data)[3] <- "value"
    pImpFreq_median_importance_data <- pImpFreq_plot_data %>% group_by(feature) %>%
      summarize(median = median(value))
    pImpFreq_median_importance_data <- pImpFreq_median_importance_data[order(pImpFreq_median_importance_data$median, decreasing = T),]
    write.csv(pImpFreq_median_importance_data, paste0(output_folder, "/pImpFreq_median_Imp_Coef_data_", m, ".csv"), row.names = T)
    colnames(pImpFreq_median_importance_data)[1] <- "Variable"
    
    if (m == "ranger"){
      ylab <- "Importance"
      title <- paste("RF, C index = ", round(cindex_median_values[cindex_median_values$Model == m,2], 3)) 
    } else {
      ylab <- "Coefficient"
      title <- paste(m, ", C index = ", round(cindex_median_values[cindex_median_values$Model == m,2], 3)) 
    }
    ggplot(pImpFreq_median_importance_data, aes(x = reorder(Variable, median), y = median)) +
      geom_col(width = 0.4) + coord_flip() +
      labs(x = NULL, y = ylab) + ggtitle(title) + theme_bw(10)
  }) # end of plot list
  names(pImpFreq_plots) <- models
  

  # Plot the variable importance from both and the c-index
  # arrange rows, cols by c-index
  library(gridExtra)
  library(grid)

  pdf(paste0(output_folder, "/variable_importance.pdf"), height = 4, width = 9)
  grid.arrange(pImpFreq_plots$lasso, pImpFreq_plots$enet, nrow = 1)
  print(pImpFreq_plots$ranger)
  print(pImpFreq_plots$lasso)
  print(pImpFreq_plots$enet)
  dev.off()

  sink()  # stop redirecting standard output
  sink(type = "message")  # stop redirecting message output
  
} # end of function performSurvivalML


