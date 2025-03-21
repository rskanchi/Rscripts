# function to create a dataframe from a list of vectors (of different lengths)
# shorter vectors will have trailing NAs in the dataframe

list2DF <- function(vList) {
  # maximum length of the vectors
  max_length <- max(sapply(vList, length))
  
  # add NAs to shorter vectors
  vList_sameLen <- lapply(vList, function(vector) {
    c(vector, rep(NA, max_length - length(vector)))
  })
  
  # dataframe 
  df <- as.data.frame(do.call(cbind, vList_sameLen))
  return(df)
} # end of function list2DF
