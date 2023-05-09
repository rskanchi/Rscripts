# function to install any packages not already installed; and load all the required packages
# input: a vector of R packages

loadRequiredPackages <- function(required.packages){
  # what packages are not installed? 
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  # install if any not installed already
  if(length(packages.not.installed)>0) {install.packages(packages.not.installed, repos = "http://cran.us.r-project.org")}
  lapply(required.packages, require, character.only = TRUE)
} # end of function loadRequiredPackages
