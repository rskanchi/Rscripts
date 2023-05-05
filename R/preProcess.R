# script to read the input trait and expr files and merge them based on the common samples
# expr data with expression rows and sample columns. If samples are rows, use the option samplesinRows = TRUE

preProcess <- function(traitDataFile, exprDataFile,
                       merged = FALSE, # the trait and expr data files are separate inputs
                       samplesinRows = FALSE, # expr data file is assumed to have samples in columns
                       prop_samples_expressed = 0.8,
                       na.list = NULL){
  
  # this function 
  func <- "preProcess"
  # NAs list before reading data
  na.list <- union(na.list, c(".", "", "Don't know","NA"))
  
  filters <- data.frame(filter = c("trait_data_samples", "expr_data_samples", "expr_data_features",
                                   "common_samples", 
                                   "after_removing_all_zero_expr", "after_removing_no_variation_expr",
                                   "prop_samples_expressed", 
                                   "expr_after_prop_filter"
  ), value = c(rep(0, 6), prop_samples_expressed, 0))
  
  # read the csv files
  traitData <- read.csv(file = traitDataFile, row.names = 1, na.strings = na.list)
  if (samplesinRows) exprData <- read.csv(file = exprDataFile, row.names = 1, na.strings = na.list) else
    exprData <- data.frame(t(read.csv(file = exprDataFile, row.names = 1, na.strings = na.list)),
                           check.names = FALSE, check.rows = FALSE) # read file and transpose so that samples are in rows
  
  filters$value[filters$filter %in% c("trait_data_samples", "expr_data_samples", 
                                      "expr_data_features")] <- c(nrow(traitData), nrow(exprData), ncol(exprData))
  
  getFunctionMessage(text = paste("Number of samples in trait file:", nrow(traitData)), func = func)
  getFunctionMessage(text = paste("Number of samples in expr file:", nrow(exprData)), func = func)
  getFunctionMessage(text = paste("Number of features in expr file:", ncol(exprData)), func = func)
  
  # merge trait and expr data; samples in both the data sets are kept
  getFunctionMessage(text = "Merging files ---------->", func = func)
  samplesCommon <- intersect(rownames(traitData), rownames(exprData))
  traitData <- traitData[samplesCommon, , drop = FALSE]
  exprData <- exprData[samplesCommon, , drop = FALSE]
  
  filters$value[filters$filter %in% c("common_samples")] <- c(length(samplesCommon))
  getFunctionMessage(text = paste("Number of samples common to trait and expr data:", length(samplesCommon)), func = func)
  
  # remove expr that are all 0s
  getFunctionMessage(text = "Removing features that are all 0s ---------->", func = func)
  expr.notZero <- colSums(exprData) != 0
  #exprZero <- colSums(exprData) == 0
  exprData <- exprData[, expr.notZero, drop = FALSE]
  filters$value[filters$filter %in% c("after_removing_all_zero_expr")] <- c(ncol(exprData))
  getFunctionMessage(text = paste("Number of features after removing all-zero features:", ncol(exprData)), func = func)
  
  # remove expr that have 0 sd (no variation)
  getFunctionMessage(text = "Removing features that have sd=0 i.e., no variation ---------->", func = func)
  expr.SDnot0 <- apply(exprData, 2, sd) != 0
  #exprSD0 <- apply(exprData, 2, sd) == 0
  exprData <- exprData[, expr.SDnot0, drop = FALSE]
  filters$value[filters$filter %in% c("after_removing_no_variation_expr")] <- c(ncol(exprData))
  getFunctionMessage(text = paste("Number of features after removing sd = 0 (no variation) features:", ncol(exprData)), func = func)
  
  # keep features expressed in more than prop_samples_expressed samples
  # Example: if prop_samples_expressed = 0.85, and 
  # if there 100 samples, keep only those features that are expressed in >=85 samples
  getFunctionMessage(text = "Checking feature expression proportion ---------->", func = func)
  expr.prop <- apply(exprData, 2, function(x) sum(x!=0 & !is.na(x))/length(x)) >= prop_samples_expressed
  exprData <- exprData[, expr.prop, drop = FALSE]
  filters$value[filters$filter %in% c("expr_after_prop_filter")] <- c(ncol(exprData))
  getFunctionMessage(text = paste("Number of features after expr_after_prop_filter:", ncol(exprData)), func = func)
  
  filters$value[c(1:6,8)] <- round(filters$value[c(1:6,8)], 0)
  # return the data files
  return(list("traitData" = traitData, "exprData" = exprData, "filters" = filters))
  
} # end of function preProcess

