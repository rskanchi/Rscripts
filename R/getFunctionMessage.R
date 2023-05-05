# messaging when in a function

getFunctionMessage <- function(text, func){
  cat(format(Sys.time(),usetz = TRUE), 
      " | Function: " , func, "() | ", 
      text, "\n", sep = "")
} # end of function getFunctionMessage 


# usage: getFunctionMessage(text = , func = func)
