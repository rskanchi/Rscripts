# setup workflow
setupWorkflow <- function(project_name){
  dir.create(project_name)
  dir.create(paste(project_name, "raw_data", sep = "/"), recursive = TRUE)
  dir.create(paste(project_name, "processed_data", sep = "/"), recursive = TRUE)
  dir.create(paste(project_name, output, sep = "/"), recursive = TRUE)
  source("/Users/rupakanchi/01projects/github/Rscripts/R/readRscripts.R")
  readRscripts(path = "/Users/rupakanchi/01projects/github/Rscripts/R")
} # end of function setupWorkflow
