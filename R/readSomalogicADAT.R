# function to read a file in Somascan's adat format and output the expression data & feature info in csv format
# REF: https://github.com/SomaLogic/SomaDataIO

readSomalogicADAT <- function(adat_file, output_prefix = "raw", output_folder = "data") {
  #install.packages("SomaDataIO")
  library(SomaDataIO)
  
  adat_data <- read_adat(adat_file) # read ADAT file
  adat_data <- adat_data[adat_data$SampleType == "Sample",]
  rownames(adat_data) <- adat_data$SampleId
  expression_data <- adat_data[, getAnalyteInfo(adat_data)$AptName, drop = FALSE] # extract expression data: sample rows, feature columns
  feature_data <- getAnalyteInfo(adat_data) # extract feature information
  
  # save to csv
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = T)
  full_adat_csv <- paste0(output_folder, "/", output_prefix, "_full_adat_data.csv")
  expression_csv <- paste0(output_folder, "/", output_prefix, "_expression_data.csv")
  feature_csv <- paste0(output_folder, "/", output_prefix, "_feature_data.csv")
  write.csv(adat_data, file = full_adat_csv)
  write.csv(expression_data, file = expression_csv, row.names = TRUE)
  write.csv(feature_data, file = feature_csv, row.names = FALSE)
  
  message("Files saved: \n", expression_csv, " \n", feature_csv)
} # end of function readSomalogicADAT

# How to use the function
#readSomalogicADAT(adat_file = adat_file_path, output_prefix = project_or_PI_name, output_folder = output_folder_path)
