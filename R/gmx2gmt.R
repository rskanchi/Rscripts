# Function to convert a gmx file to gmt file
# Mitocarta 3.0 pathway file in in gmx format

gmx2gmt <- function(gmx.file){
  # Read the GMX file
  gmx <- fread(gmx.file, header = FALSE) # is read as a list
  gmt <- lapply(gmx, function(pathway) {
    pathway <- pathway[pathway != ""] # remove trailing empty strings ("") from each list element
    pathway # in gmt format: [pathway name, description, Gene1, Gene2, ...]
  })
  
  # Write to a gmt file
  gmt_file <- gsub("gmx", "gmt", gmx.file)
  fileConn <- file(gmt_file, "w") # connection to write the GMT file
  
  for (i in 1:length(gmt)) {
    # pathway name, description, and tab-separated genes
    gmt_line <- paste(c(gmt[[i]][1], gmt[[i]][2], gmt[[i]][-c(1,2)]), collapse = "\t")
    writeLines(gmt_line, fileConn) # write the line to the file
  }
  close(fileConn)
} # end of gmx2gmt