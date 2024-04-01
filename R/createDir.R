# function to create a folder recursively, if it doesn't exist

createDir <- function(folder = NULL) {
  if(is.null(folder)) dir.create("output") else
  if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
}