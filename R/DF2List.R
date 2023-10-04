# convert data frame columns to list elements and remove the NAs

DF2List <- function(df){
  lapply(df, FUN = function(x) x[!is.na(x)])
}