# function to read R scripts from a specified folder

readRscripts <- function(path){
  rFiles <- list.files(path = path, pattern = "*.R")
  sapply(rFiles, FUN = function(rF) source(paste(path, rF, sep = "/")))
  
  message(paste("R scripts sourced from the folder", path))
}
