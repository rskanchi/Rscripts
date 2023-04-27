#  function to check if an object is a positive integer

is.positive.integer <- function(x){
  if (is.integer(x) & (x > 0)) flag <- TRUE else flag <- FALSE
  if (!flag) print(paste(x, "is not a positive integer | error from function is.positive.integer()"))
  if (!flag) getFunctionMessage(paste("Running function is.positive.integer() |", x, "is not a positive integer")) 
  return(flag)
} # end of function is.positive.integer


