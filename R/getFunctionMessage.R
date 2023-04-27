# messaging when in a function

getFunctionMessage <- function(text, func){
  cat(format(Sys.time(),usetz = TRUE), 
      " | Running function: " , func, "() \n", 
      text, "\n", sep = "")
} # end of function getFunctionMessage 


