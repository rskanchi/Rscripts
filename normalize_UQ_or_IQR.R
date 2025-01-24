# function to perform upper-quartile (UQ) or interquartile range (IQR) normalization

normalize_UQ_or_IQR <- function(d, method = c("UQ", "IQR")) {
  # get method
  method <- match.arg(method)
  
  # Ensure the input is numeric
  if (!is.numeric(d)) {
    stop("The input data must be numeric.")
  }
  
  # normalization based on the chosen method
  if (method == "UQ") {
    # UQ normalization
    uq_values <- apply(d, 2, function(column) quantile(column, 0.75, na.rm = TRUE))
    uq_values[uq_values == 0] <- 1e-6  # Avoid division by zero
    normalized_data <- sweep(d, 2, uq_values, FUN = "/")
  } else if (method == "IQR") {
    # IQR normalization
    quartiles <- apply(d, 2, function(column) quantile(column, probs = c(0.25, 0.75), na.rm = TRUE))
    iqr_values <- quartiles[2, ] - quartiles[1, ]
    iqr_values[iqr_values == 0] <- 1e-6  # to avoid division by zero
    normalized_data <- sweep(d, 2, iqr_values, FUN = "/")
  }
  
  # Return the normalized matrix
  return(normalized_data)
}
