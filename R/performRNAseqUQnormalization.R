# function to perform upper-quartile (UQ) normalization for RNA-seq data
# input: RNA-seq counts in csv file 
# the functions filters lowly expressed genes and performs UQ normalization
library(edgeR)
performRNAseqUQnormalization <- function(data_file, min_CPM = 1, min_samples = 2) {
  # Step 1: Load the data
  raw_data <- read.csv(data_file)  
  gene_info <- raw_data[, 1:3]  # First three columns are gene information
  count_data <- raw_data[, -(1:3)]  # Remaining columns are sample counts
  
  # Step 2: Filter lowly expressed genes
  keep <- rowSums(cpm(count_data) > min_CPM) >= min_samples
  filtered_count_data <- count_data[keep, ]
  filtered_gene_info <- gene_info[keep, ]  # Subset gene info to match filtered rows
  
  # Step 3: Perform upper-quartile (UQ) normalization
  dge <- DGEList(counts = filtered_count_data)
  dge <- calcNormFactors(dge, method = "upperquartile")
  normalized_data <- cpm(dge, normalized.lib.sizes = TRUE)
  
  # Combine normalized data with gene information for output
  output_data <- cbind(filtered_gene_info, normalized_data)
  
  # Return the combined data
  return(output_data)
}
