# Function to generate pathway scores using gene expression data
# Input 
    # specify either gene expression data or the csv file name; expression data format = samples x genes (row x col)
    # information: 
        # species: default human; use the function msigdbr_species() to see the available species in msigdb
        # category: default hallmark "H"; use the function msigdbr_collections() 
        # subcategory: this is optional, some categories (like C2) have sub-categories (such as KEGG, REACTOME..)and some don't (like Hallmark)
        # so check the reference link https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp to specify category and sub-category

# Output is an excel file with several worksheets (change to csv files at a later time)
    # overlap genes in each module/cluster and frequency of overlap genes
    # ucellScores: 
        # matrix of signature enrichment scores for each sample; sample x pathway scores (col x row)
        # Method: Ucell scores, the scores are normalized Mann-Whitney U statistics (0 to 1) and 
        # depend on the gene expression ranks of individual sample
        # Rank based, robust across datasets irrespective of dataset composition

    # AddGenesetScores: matrix with overlapping geneset scores added

getPathwayScores.msigdb <- function(geneExprData = NULL,
                             geneExprFile = NULL,
                             species = "Homo sapiens", 
                             category = "H",
                             subcategory = NULL,
                             outFolder = NULL
                             ){
  
  library(UCell)
  library(msigdbr)
  if (is.null(geneExprData) & !is.null(geneExprFile)){geneExprData <- read.csv(file = geneExprFile, row.names = 1)}

  # msigdbr_species()
  # collections <- msigdbr_collections()
  msigdb.df <- msigdbr(species = species, category = category, subcategory = subcategory) # category subcategory coding later
  collection.pathways <- unique(msigdb.df$gs_name)
  
  geneOverlapFreq <- data.frame(matrix(ncol = 3, nrow = length(collection.pathways))) # output
  colnames(geneOverlapFreq) <- c("pathway", "nPathwayGenes", "nOverlap")
  overlapGenes <- list()

  AddGenesetScores <- matrix(nrow = nrow(geneExprData), ncol = length(collection.pathways))
  rownames(AddGenesetScores) <- rownames(geneExprData); colnames(AddGenesetScores) <- collection.pathways
  
  for (p in 1:length(collection.pathways)){
    p.genes <- msigdb.df$gene_symbol[msigdb.df$gs_name == collection.pathways[p]]
    
    overlapGenes[[collection.pathways[p]]] <- intersect(colnames(geneExprData), p.genes)
    geneOverlapFreq[p,] <- c(collection.pathways[p], length(p.genes), length(overlapGenes[[collection.pathways[p]]]))
    AddGenesetScores[, p] <- rowSums(as.matrix(geneExprData[, overlapGenes[[collection.pathways[p]]]]))
  } # end of for loop in p
  
  ucellScores <- ScoreSignatures_UCell(t(geneExprData), features = overlapGenes)
  overlapGenesMatrix <- data.frame(lapply(overlapGenes, function(g) {
    g <- unlist(g)
    length(g) <- max(lengths(overlapGenes))
    return(g)
  }))
  
  if (!is.null(outFolder)) {
    createDir(outFolder)
    library(openxlsx)
    write.xlsx(list("msigDB" = msigdb.df, "geneOverlapFreq" = geneOverlapFreq,
                    "overlapGenes" = overlapGenesMatrix, "ucellScores" = ucellScores, "AddGenesetScores" = AddGenesetScores),
               file = paste(outFolder, "msigDB_PathwayScores.xlsx", sep = "/"),
               rowNames = c(FALSE, FALSE, FALSE, TRUE, TRUE), overwrite = TRUE)
  } # end of writing output 
  

  return(list("species" = species, category = "H", "msigDB" = msigdb.df, "geneOverlapFreq" = geneOverlapFreq,
         "overlapGenes" = overlapGenesMatrix, "ucellScores" = ucellScores, "AddGenesetScores" = AddGenesetScores))

} # end of function getPathwayScores.msigdb

