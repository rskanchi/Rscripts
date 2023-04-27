# function to create a folder recursively, if it doesn't exist

createDir <- function(folder) {
  if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}       
}