# function to check if numeric traits in a dataframe were read as characters and to convert them back to numeric

numericFlag <- function(d){
  flag <- suppressWarnings(sapply(d, function(x){!any(is.na(as.numeric(x)))}))
  d[,flag] <- sapply(d[,flag], as.numeric)
  return(d)
}
