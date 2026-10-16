# Function to convert a gmx file to gmt file
# Mitocarta 3.0 pathway file is in gmx format
library(data.table)
gmx2gmt <- function(gmx.file){
  # Read the GMX file (one pathway per column)
  gmx <- fread(gmx.file, header = FALSE, sep = "\t", fill = TRUE,
               colClasses = "character")
  gmt <- lapply(gmx, function(pathway) {
    pathway <- pathway[!is.na(pathway) & pathway != ""] # drop padding/NA
    c(pathway[1:2], unique(pathway[-c(1, 2)]))         # name, description, unique genes
  })
  
  # Write to a gmt file (swap only the extension)
  gmt_file <- sub("\\.gmx$", ".gmt", gmx.file)
  fileConn <- file(gmt_file, "w")
  
  for (i in seq_along(gmt)) {
    writeLines(paste(gmt[[i]], collapse = "\t"), fileConn)
  }
  close(fileConn)
  
  invisible(gmt_file)  # return the path so you can pipe it into gmtPathways()
} # end of gmx2gmt