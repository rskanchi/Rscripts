# function to get consensus clusters of the columns of a dataframe

getConsensusClusters <- function(d, folder, dist.method = "spearman", clusterAlg = "hc", 
                                 innerLinkage = "ward.D", finalLinkage = "ward.D",
                                 maxK = 8, nReps = 500, pItem = 0.9, pFeature = 1, corUse = "pairwise.complete.obs", 
                                 seed = 7654321){
  library(ConsensusClusterPlus)
  library(clustree)
  library(ggplot2)
  
  if (!dir.exists(folder)){ dir.create(folder, recursive = TRUE) }
  
  if (dist.method == "spearman" | dist.method == "pearson")
    distMat <- as.dist(1-cor(d, use = corUse, method = dist.method)) else # cor() returns cor between cols of dataframe
      distMat <- as.dist(dist(t(d), method = dist.method)) # dist() returns dist between rows of the dataframe so t(d)
    
    write.csv(as.matrix(distMat), file = paste(folder, paste0("distMat_", dist.method, ".csv"), sep = "/"), row.names = TRUE)
    
    consensusRes <- ConsensusClusterPlus(distMat, maxK = maxK, reps = nReps, pItem = pItem, pFeature = pFeature,
                                         innerLinkage = innerLinkage, finalLinkage = finalLinkage, clusterAlg = clusterAlg, 
                                         title = folder, seed = seed, plot = "png")
    # cluster tree for expression/features
    df_clusters <- data.frame(row.names = colnames(d))
    df_clusters[, 1] <- rep(1, ncol(d))
    for (k in 2:maxK){ df_clusters[, k] <- consensusRes[[k]][["consensusClass"]] }
    colnames(df_clusters) <- paste("K", 1:maxK, sep = "")
    write.csv(df_clusters, paste(folder, paste0("clustree_", dist.method, ".csv"), sep = "/"), row.names = TRUE)
    
    #df_clusters <- df_clusters[, c(1:12, 18, 24)]
    clusterTree <- clustree(df_clusters, prefix = "K", node_label = "size", prop_filter = 0, node_text_size = 5, node_label_size = 5,
                            node_label_nudge = -0.4) + 
      guides(fill="none", x = "none", y = "none") + labs(x = NULL, y = NULL) +
      theme_classic(15) # use_core_edges = FALSE, edge_arrow = TRUE, edge_arrow_ends = "both"
    
    pdf(paste(folder, paste0("clustree_", dist.method, ".pdf"), sep = "/"), height = 10, width = 10)
    print(clusterTree) 
    dev.off()
    
    return(consensusRes)
    
} # end of function getConsensusClusters
