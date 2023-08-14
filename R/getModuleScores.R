# Function to generate enrichment scores using gene expression data based on an input gene list
# Input 
  # specify either gene expression data or the csv file name; expression data format = samples x genes (row x col)

  # Output is an excel file with several worksheets (change to csv files at a later time)
  # overlap genes in each module/cluster and frequency of overlap genes
  # ucellScores: 
  # matrix of enrichment scores for each sample; sample x geneset scores (col x row)
  # Method: Ucell scores, the scores are normalized Mann-Whitney U statistics (0 to 1) and 
  # depend on the gene expression ranks of individual sample
  # Rank based, robust across datasets irrespective of dataset composition
  
  # AddGenesetScores: matrix with geneset scores added


getModuleScores <- function(geneExprData = NULL,
                            geneExprFile = NULL,
                            geneList,
                            outFolder = NULL
                            
){
  
  library(UCell)
  if (is.null(geneExprData) & !is.null(geneExprFile)){geneExprData <- read.csv(file = geneExprFile, row.names = 1)}
  geneOverlapFreq <- data.frame(matrix(ncol = 3, nrow = length(geneList))) # output
  colnames(geneOverlapFreq) <- c("module", "nModuleGenes", "nOverlap")
  overlapGenes <- list()
  
 
  AddGenesetScores <- matrix(nrow = nrow(geneExprData), ncol = length(geneList))
  rownames(AddGenesetScores) <- rownames(geneExprData)
  if (is.null(names(geneList))) geneListnames <- 1:length(geneList) else geneListnames <- names(geneList)
  colnames(AddGenesetScores) <- geneListnames
  
  for (m in 1:length(geneListnames)){
    m.genes <- geneList[[geneListnames[m]]]
    
    overlapGenes[[geneListnames[m]]] <- intersect(colnames(geneExprData), m.genes)
    geneOverlapFreq[m,] <- c(geneListnames[m], length(m.genes), length(overlapGenes[[geneListnames[m]]]))
    AddGenesetScores[, m] <- rowSums(geneExprData[, overlapGenes[[geneListnames[m]]]])
  } # end of for loop in m
  
  ucellScores <- ScoreSignatures_UCell(t(geneExprData), features = overlapGenes)
  overlapGenesMatrix <- data.frame(lapply(overlapGenes, function(g) {
    g <- unlist(g)
    length(g) <- max(lengths(overlapGenes))
    return(g)
  }))
  colnames(overlapGenesMatrix) <- geneListnames
  
  if (!is.null(outFolder)) {
    createDir(outFolder)
    library(openxlsx)
    write.xlsx(list("geneOverlapFreq" = geneOverlapFreq, "overlapGenes" = overlapGenesMatrix, 
                    "ucellScores" = ucellScores, "AddGenesetScores" = AddGenesetScores),
               file = paste(outFolder, "geneModuleScores.xlsx", sep = "/"),
               rowNames = c(FALSE, FALSE, TRUE, TRUE), overwrite = TRUE)
  } # end of writing output 
  
  return(list("geneOverlapFreq" = geneOverlapFreq,
              "overlapGenes" = overlapGenesMatrix, "ucellScores" = ucellScores, "AddGenesetScores" = AddGenesetScores))
  
} # end of function getModuleScores

