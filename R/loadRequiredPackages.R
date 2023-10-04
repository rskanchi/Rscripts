# function to install any packages not already installed; and load all the required packages
# input: a vector of R packages

loadRequiredPackages <- function(required.packages){
  # what packages are not installed? 
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
  
  
  # install if any not installed already
  if(length(packages.not.installed) > 0) {
    # if installation through biocManager
    biocPackages <- packages.not.installed[packages.not.installed %in% BiocManager::available()]
    if (length(biocPackages) > 0) {
      BiocManager::install(biocPackages)
      packages.not.installed <- setdiff(packages.not.installed, biocPackages)
    } # end of if in length(biocPackages) > 0
    
    # after installing biocPackages, if any
    install.packages(packages.not.installed, repos = "http://cran.us.r-project.org")
    } # end of if in length(packages.not.installed)>0
  lapply(required.packages, require, character.only = TRUE)
} # end of function loadRequiredPackages
