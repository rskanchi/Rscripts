
# function to perform the Weighted Gene Correlation Network Analysis

runWGCNA <- function(expr_file,
                     trait_file,
                     na.list = NULL,
                     folder,
                     sample_names_file = NULL,
                     feature_names_file = NULL,
                     corMethod = "pearson",
                     padjMethod = "BH",
                     scaling = T,
                     traitScaling = F,
                     traitVarFilter = 0.1,
                     feature_selection = NULL,
                     softPower = NULL,
                     networkType = "signed", 
                     clustMethod = "average",
                     min_module_size = c(5, 10, 25, 50, 100, 250, 500, 1000),
                     softCutoff = 0.9, # adjust this based on the sft$fitIndices table; corresponds to using an R^2 cut-off
                     MEDissThres = 0.1, # height cut of 0.1, corresponding to correlation of 0.9
                     # module-trait correlation heatmap
                     MTCHlabSize = 0.5,
                     MTCHxSize = 0.6,
                     
                     # intra-trait correlation heatmap
                     ITCHlabSize = 2,
                     ITCHtextSize = 9,
                     ITCHtitleSize = 10,
                     barplotYLabSize = 8
                     ) {
  
  library(ComplexHeatmap)
  
  # Creating separate output directory to store data in
  if (!dir.exists(paste(folder,"data",sep = "/"))){dir.create(paste(folder,"data",sep = "/"), recursive = TRUE)}
  
  sink(paste0(folder,"/log.Rout"), type = c("output", "message"))
  
  # Loading required packages and installing ones not present
  list.of.packages <- c("WGCNA", "ggplot2", "ggrepel", "gridExtra", "dplyr", "openxlsx","ggcorrplot", "circlize")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0) {install.packages(new.packages)} else {lapply(list.of.packages, require, character.only = TRUE)}
  
  # Adding any user specified NA strings to the default na.list
  na.list <- union(na.list,c("","Don't know","NA"))
  
  # Reading in the expression and clinical data files
  rawExprData <- data.frame(t(read.csv(expr_file, na.strings = na.list, row.names = 1)), check.names = F, check.rows = F)
  rawTraitData <- read.csv(trait_file, na.strings = na.list, row.names = 1)
  
  cat(getTime(), " Number of samples (raw): ", nrow(rawExprData),"\n", sep = "")
  cat(getTime(), " Number of features (raw): ", ncol(rawExprData),"\n", sep = "")
  cat(getTime(), " Number of traits (raw): ", ncol(rawTraitData),"\n", sep = "")
  
  # dropping traits with more than half the values missing
  keep <- !(apply(rawTraitData, 2, function(x){sum(is.na(x))}) > (ncol(rawTraitData)) * 0.5)
  rawTraitData <- rawTraitData[, keep, drop=F]
  
  # Restricting expression data to samples and features specified by the user
  if(!(is.null(sample_names_file))) {
    keep_samps <- unname(unlist(read.csv(sample_names_file)))
    rawExprData <- rawExprData[keep_samps, ]
  }
  
  if(!(is.null(feature_names_file))) {
    keep_ftrs <- unname(unlist(read.csv(feature_names_file)))
    rawExprData <- rawExprData[, keep_ftrs]
  }
  
  # Removing zero-variance genes and samples with missing entries (NAs)
  gsg <- goodSamplesGenes(rawExprData, verbose = 0)
  rawExprData <- rawExprData[gsg$goodSamples, gsg$goodGenes]
  
  # Feature selection by variance
  if(!(is.null(feature_selection))) {
    if(ncol(rawExprData) < feature_selection) {stop(cat(getTime(), "You have asked to select more features (", feature_selection,") than available in the dataset (", nrow(rawExprData),"). Please try again with a lower number of features to select.", "\n", sep = ""))}
    rawExprData <- rawExprData[ , order(apply(rawExprData, 2, var), decreasing = T)]
    rawExprData <- rawExprData[ , 1:feature_selection]
  }
  
  # Restricting minModuleSize to maximum size of n/2, where n is the number of features in the expression data
  min_module_size <- min_module_size[min_module_size < (ncol(rawExprData)/2)]
  cat(getTime(), "The minModuleSize range is: ", min_module_size, "\n")
  
  # Restricting clinical data to samples included in the expression data and matching the sample order between the datasets
  rawTraitData <- rawTraitData[rownames(rawTraitData) %in% rownames(rawExprData), , drop=F]
  rawTraitData <- rawTraitData[match(rownames(rawExprData), rownames(rawTraitData)), , drop=F]
  if (!(identical(rownames(rawTraitData), rownames(rawExprData)))) {stop("Sample names in the expression and clinical file are not matching up, please check. \n")}
  
  # Ensuring no numerical traits were read in as strings
  flag <- suppressWarnings(sapply(rawTraitData, function(x){!any(is.na(as.numeric(x)))}))
  rawTraitData[,flag] <- sapply(rawTraitData[,flag], as.numeric)
  
  # Dropping multi class traits
  ## Selecting traits with more than 2 classes but less than 6, these are multi class traits and need to be dropped
  keep <- ((sapply(rawTraitData, function(x){length(unique(x[!(is.na(x))]))}) == 2) | (sapply(rawTraitData, function(x){length(unique(x[!(is.na(x))]))}) > 6))
  cat(getTime(),"Dropping following columns for being discrete or having variance of 0: \n", colnames(rawTraitData)[!keep], "\n")
  rawTraitData <- rawTraitData[,keep,drop=F]

  cat(getTime(), " Number of samples (filtered): ",nrow(rawExprData), "\n", sep = "")
  cat(getTime(), " Number of features (filtered): ",ncol(rawExprData), "\n", sep = "")
  cat(getTime(), " Number of traits (filtered): ",ncol(rawTraitData), "\n", sep = "")
  
  # Saving the raw traits file
  write.csv(rawTraitData, paste(folder,"data/traits_raw.csv",sep = "/"), row.names = T)
  
  # Scaling numerical traits if set to True
  if (traitScaling) {
    varsToScale <- sapply(rawTraitData, function(x){length(unique(x[!(is.na(x))]))}) > 2
    rawTraitData[,varsToScale] <- scale(rawTraitData[, varsToScale, drop=F]) 
  }
  
  # Converting binary traits to 0s and 1s
  flag <- sapply(rawTraitData, function(x){length(unique(x[!(is.na(x))]))}) == 2
  traitData <- rawTraitData
  traitData[,flag] <- sapply(traitData[,flag, drop=F], function(x){as.numeric(as.factor(x))})
  traitData[,flag] <- apply(traitData[,flag, drop=F], 2, function(x){ifelse(x==1, 0, 1)})
  
  # Saving the processed traits file
  write.csv(traitData, paste(folder,"data/traits_fixed.csv",sep = "/"), row.names = T)
  
  # Saving the raw data
  write.csv(rawExprData, paste(folder,"data/exprs_raw.csv",sep = "/"), row.names = T)
  
  # Performing scaling if option set to T otherwise ignoring it
  if (isTRUE(scaling)) {
    d <- data.frame(scale(rawExprData), check.names = F)
    # Saving the normalized data
    write.csv(d, paste(folder,"data/exprs_scaled.csv",sep = "/"), row.names = T)
  } else {
    d <- rawExprData
  }
  
  nSamples <- nrow(d)
  nClin <- ncol(traitData)
  clinNames <- colnames(traitData)
  
  # WGCNA to explore co-expressed expression modules and their association with clinical variables
  # Module/network construction
  
  enableWGCNAThreads()
  # Step 0: Soft thresholding - choose a set of soft-thresholding powers #--------------------
  
  powers <- c(c(1:10), seq(from = 12, to=20, by=2))
  sft <- pickSoftThreshold(d, powerVector = powers, verbose = 0)
  # sft$fitIndices
  # sft$powerEstimate
  # Plot the indices
  sft_plot1 <- ggplot(data = sft$fitIndices, aes(x = Power, y = -sign(slope)*SFT.R.sq)) +
    geom_point(size = 0.5, color = "red") +
    geom_text_repel(label = sft$fitIndices$Power, size = 4, color = "red") +
    xlab("Soft threshold (power)") +
    ylab("Scale-free topology (signed R^2)") +
    labs(title = "Scale independence") +
    geom_hline(yintercept = softCutoff, linetype = "dashed", color = "red") +
    theme_bw(15)
  
  sft_plot2 <- ggplot(data = sft$fitIndices, aes(x = Power, y = mean.k.)) +
    geom_point(size = 0.5, color = "red") +
    xlab("Soft Threshold (power)") +
    ylab("Mean Connectivity") +
    labs(title = "Mean connectivity") +
    geom_text_repel(label = sft$fitIndices$Power, size = 4, color = "red", direction = "y") +
    theme_bw(15)
  
  ## results in a folder
  
  if (!dir.exists(folder)){ dir.create(folder, recursive = TRUE) }
  write.csv(sft$fitIndices, file = paste(folder, "soft_thresholding.csv", sep = "/"))
  pdf(paste(folder, "plots_soft_thresholding.pdf", sep = "/"), width = 12)
  grid.arrange(sft_plot1, sft_plot2, ncol = 2)
  dev.off()
  

  if(is.null(softPower)) {
    softPower <- sft$powerEstimate 
  }
  
  if(is.na(softPower)) {
    powerdf <- sft$fitIndices
    softPower <- powerdf$Power[powerdf$SFT.R.sq==max(powerdf$SFT.R.sq)]
    cat(getTime(), "Power is set to ", softPower,"\n\n")
  } else {cat(getTime(), "Power is set to ", softPower,"\n\n")}
  
  # Step 1: inter trait correlation (Only for datasets with more than 2 clinical traits) #---------------
  
  if (ncol(traitData)>2) {
    
    # dropping traits with more than half the values missing
    temp.traitData <- traitData[,!(apply(traitData, 2, function(x){sum(is.na(x))}) > (ncol(traitData))/2)]
    # dropping traits with variance < traitVarFilter
    temp.traitData <- temp.traitData[,!(apply(temp.traitData, 2, function(x){var(x[!is.na(x)])}) < traitVarFilter)]
    
    if(ncol(temp.traitData)>2) {
      cormat <- cor(temp.traitData, use = "p", method = corMethod)
      # pmat <- cor_pmat(temp.traitData)
      pmat <- data.frame(matrix(data = 0, nrow = ncol(temp.traitData), ncol = ncol(temp.traitData), dimnames = list(colnames(temp.traitData), colnames(temp.traitData))), check.names = F, check.rows = F)
      for (x in 1:ncol(temp.traitData)) {
        for (y in 1:ncol(temp.traitData)) {
          pmat[x,y] <- cor.test(temp.traitData[,x],temp.traitData[,y], method = corMethod)[["p.value"]]
        }
      }
      ncols <- ncol(pmat); n <- ncols*ncols - ncols - (ncols*(ncols-1)/2) # the number of comparisons
      
      padjmat <- apply(pmat, c(1,2), FUN = function(x) p.adjust(x, method = padjMethod, n = n))
      
      corTitle <- paste("Inter Trait Correlation, n =", nrow(temp.traitData))
      
      pdf(paste(folder, "inter_trait_corr_heatmap.pdf", sep = "/"))
      
      print(
        ggcorrplot(cormat, hc.order = TRUE, type = "lower", lab = TRUE, lab_size = ITCHlabSize, tl.cex = ITCHtextSize, 
                   title = corTitle) + 
          theme(plot.title = element_text(size = ITCHtitleSize))
      )
      print(
        ggcorrplot(cormat, hc.order = TRUE, type = "lower", lab = TRUE, lab_size = ITCHlabSize, tl.cex = ITCHtextSize, 
                   insig = "blank", p.mat = pmat, 
                   title = paste(corTitle, ", p < 0.05 (p-adjust method = none)")) + 
          theme(plot.title = element_text(size = ITCHtitleSize))
      )
      print(
        ggcorrplot(cormat, hc.order = TRUE, type = "lower", lab = TRUE, lab_size = ITCHlabSize, tl.cex = ITCHtextSize,
                   insig = "blank", p.mat = padjmat, sig.level = 0.25, 
                   title = paste0(corTitle, ", p < 0.25 (p-adjust method = ", padjMethod,")")) + 
          theme(plot.title = element_text(size = ITCHtitleSize))
      )
      print(
        ggcorrplot(cormat, hc.order = TRUE, type = "lower", lab = TRUE, lab_size = ITCHlabSize, tl.cex = ITCHtextSize, 
                   insig = "blank", p.mat = padjmat, sig.level = 0.05, 
                   title = paste0(corTitle, ", p < 0.05 (p-adjust method = ", padjMethod, ")")) + 
          theme(plot.title = element_text(size = ITCHtitleSize))
      )
      
      dev.off()
      
      list_of_datasets <- list("correlation" = cormat, "pvalues" = pmat, "FDR" = padjmat)
      write.xlsx(list_of_datasets, file = paste(folder, "inter_trait_corr_heatmap.xlsx", sep = "/"), overwrite = TRUE, rowNames = TRUE)
      #Also save the long form of these data
      To_long_form <- function(matrix, type){
        pairs <-t(combn(colnames(matrix),2))
        pairwise_corr <- data.frame(pairs, name =matrix[pairs])
        names(pairwise_corr)[3] <- type
        return(pairwise_corr)
      }
      
      list_of_datasets_long <- list("correlation"=To_long_form(cormat, "correlation"),
                                    "pvalues"= To_long_form(pmat, "pvalues"),
                                    "FDR"= To_long_form(padjmat,"FDR"))
      write.xlsx(list_of_datasets_long, file = paste(folder, "inter_trait_corr_heatmap_long.xlsx", sep = "/"), overwrite = TRUE, rowNames = FALSE)
    } else {cat(getTime(), "Either 1 trait or none passing the variance filter or having less than half the values missing so inter trait correlation not possible \n\n")}
    
  } else {cat(getTime(), "Either 1 trait or none remaining so inter trait correlation not possible \n\n")}
  
  
  # Step 2: Build network #---------------
  
  # Similarity/distance between each pair of expr computed based on shared neighbors (network/topological overlap)
  # Topological overlap measures interconnectedness
  # To define the adjacency matrix, a soft threshold power `r softPower` was chosen, the lowest power (`r softCutoff`) for
  # which the scale-free topology fit index curve flattens
  # for network type https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
  # signed or signed hybrid: https://peterlangfelder.com/2018/11/25/__trashed/

  adjacency <- adjacency(d, corFnc = "cor", corOptions = list(use = 'p', method = corMethod),
                         type = networkType, power = softPower)
  # to minimize effects of noise and spurous associations, transform the adjacency into topological overlap matrix (TOM)
  #calculate corresponding dissimilarity
  #TOM <- TOMsimilarity(adjacency); dissTOM <- 1-TOM
  dissTOM <- TOMdist(adjacency, verbose = 0, TOMType = networkType)
  geneTree <- hclust(as.dist(dissTOM), method = clustMethod)
  
  for (minModuleSize in min_module_size) {
    tryCatch({
      cat(getTime(), "Performing WGCNA for minimum module size of ",minModuleSize,"\n")
      
      dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 3, pamRespectsDendro = FALSE,
                                   minClusterSize = minModuleSize, verbose = 0)
      #table(dynamicMods) # Label 0 (grey) is for unassigned
      nMods <- length(table(dynamicMods))-1
      dynamicColors <- labels2colors(dynamicMods) # Convert numeric labels into colors
      #table(dynamicColors)
      
      # Plot the dendrogram and colors underneath
      clust_plot <- plotDendroAndColors(geneTree, dynamicColors, "Modules",
                                        dendroLabels = FALSE, hang = 0.03,
                                        addGuide = TRUE, guideHang = 0.05,
                                        main = "clusters and module colors")
      
      # merging modules whose expression profiles are similar
      # To quantify co-expression similarity of entire modules,
      # calculate their eigengenes and cluster them on their correlation
      MEList <- moduleEigengenes(d, colors = dynamicColors, softPower = softPower, scale = FALSE) # already scaled data so turning off the scale option
      MEs <- MEList$eigengenes
      MEDiss <- 1-cor(MEs, use = "p", method = corMethod) # Calculate dissimilarity of module eigengenes
      
      # Spearman ERROR: Error in r + t(r) - diag(diag(r)) : non-conformable arrays
      
      METree <- hclust(as.dist(MEDiss), method = clustMethod) # Cluster module eigengenes
      
      moduleColors <- dynamicColors # Save module colors and labels for use in subsequent parts
      colorOrder <- c("grey", standardColors(50)) # Construct numerical labels corresponding to the colors
      moduleLabels <- match(moduleColors, colorOrder)-1
      
      # Smaller modules with similar expression profiles merged?
      
      merge <- mergeCloseModules(d, dynamicColors, cutHeight = MEDissThres, verbose = 0)
      mergedColors <- merge$colors # merged module colors
      mergedMEs <- merge$newMEs # eigengenes of the new merged modules
      # visualize result of merging
      # plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Initial clusters", "Merged clusters"),
      #                     dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
      
      # the same as initial in this case
      # Rename to moduleColors
      moduleColors <- mergedColors
      # Construct numerical labels corresponding to the colors
      colorOrder <- c("grey", standardColors(50));
      moduleLabels <- match(moduleColors, colorOrder)-1;
      MEs <- mergedMEs
      nMods <- ncol(MEs)
      
      # Step 3: module-trait association #---------------
      
      
      # For each module, an eigengene is computed which is the first principal component based on expr of members in that module.
      # Eigengenes can be thought of as weighted average expression for the co-expression modules.
      # Rows in the figure below correspond to eigengenes of the modules and the columns to the clinical variables.
      # Modules significantly associated with the clinical variales can be further studied.
      
      # Recalculate MEs with color labels
      MEs0 <- moduleEigengenes(d, colors = moduleColors, softPower = softPower, scale = FALSE)$eigengenes
      MEs <- orderMEs(MEs0)
      # If you want traits in a specific order, specify it here - like cluster them?
      # not required for this data since the traits are in a predefined order
      
      # ---------------------- code below for a trait order in module-trait cor heatmap ----------------------
      #traitTree <- hclust(dist(t(traitData)), method = clustMethod) #, method = "canberra", "average"
      #plot(traitTree, main = "Trait clustering to order them", sub="", xlab="", cex.lab = 1.8,
      #     cex.axis = 1.5, cex.main = 1.8) # rev, sort to get an order
      # ----------------------------------------------------------------------------------------
      
      moduleTraitCor <- cor(MEs, traitData[, clinNames, drop=F], use = "p", method = corMethod); # clinNames[traitTree$order]
      moduleTraitCor[is.na(moduleTraitCor)] <- 0
      moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
      textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
      dim(textMatrix) <- dim(moduleTraitCor)
      yLabels <- substring(rownames(moduleTraitCor),3)
      ySymbols <- table(mergedColors)[yLabels]
      ySymbols <- paste(yLabels, "(" , ySymbols, ")" , sep = "")
      moduleMembership <- data.frame(exprId = colnames(d), moduleColor = moduleColors)
      
      # Creating separate output directory for each module size parameter
      if (!dir.exists(paste(folder,paste0("minModuleSize",minModuleSize),sep = "/"))){
        dir.create(paste(folder,paste0("minModuleSize",minModuleSize),sep = "/"), recursive = TRUE)
      }
      
      pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,".pdf"), sep = "/"), width = 14, height = 10)
      # Display the correlation values within a heatmap plot
      mar.def <- par()$mar
      par(mar = mar.def + c(0,2,0,0))
      labeledHeatmap(Matrix = moduleTraitCor,
                     xLabels = colnames(moduleTraitCor),
                     yLabels = rownames(moduleTraitCor),
                     ySymbols = ySymbols,
                     colorLabels = FALSE,
                     colors = blueWhiteRed(50),
                     textMatrix = textMatrix,
                     setStdMargins = FALSE,
                     cex.text = MTCHlabSize,
                     cex.lab.x = MTCHxSize,
                     cex.lab.y = 0.6,
                     zlim = c(-1,1),
                     main = paste("Module-trait correlations"))
      dev.off()
      par(mar = mar.def)
      write.csv(moduleMembership, file = paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("moduleMembership_MMS",minModuleSize,".csv"), sep = "/"), row.names = FALSE)
      
      # Saving excel file of data used to create the module trait correlation heatmap
      moduleTraitCombined <- (function(df.corr, df.pval){
        rownames(df.corr) <- gsub("ME","", rownames(df.corr))
        rownames(df.pval) <- gsub("ME","", rownames(df.pval))
        colnames(df.corr) <- paste0(colnames(df.corr),"_CORR")
        colnames(df.pval) <- paste0(colnames(df.pval),"_PVAL")
        finaldf <- c()
        for (idx in 1:ncol(df.corr)) {
          finaldf <- cbind(finaldf, df.corr[,idx], df.pval[,idx])
          colnames(finaldf)[c(2*idx-1, 2*idx)] <- c(colnames(df.corr[,idx, drop=F]), colnames(df.pval[,idx, drop=F]))
        }
        return(finaldf)
      })(moduleTraitCor, moduleTraitPvalue)
      write.xlsx(list(Combined=moduleTraitCombined, CorrCoeff=moduleTraitCor, R2=moduleTraitCor^2, Pvalues=moduleTraitPvalue), paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,".xlsx"), sep = "/"), row.names = T, overwrite = T, headerStyle = createStyle(textRotation = 45), colWidths = "auto")
      
      # scatter plot showing the module-trait correlations
      # separate file for these scatter plots (too many)
      pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("scatterplots_ME_vs_trait_MMS",minModuleSize,".pdf"), sep = "/"))
      dScatter <- data.frame(MEs, traitData)
      for (m in colnames(MEs)){
        for (trt in clinNames){ # clinNames[traitTree$order]
          print(ggplot(data = dScatter,
                       aes(x = dScatter[,m], y = dScatter[,trt])) +
                  labs(title = "module eigengene (ME) x trait scatter plot",
                       subtitle = paste("cor = ", round(moduleTraitCor[m,trt],3), "; p-val = ", round(moduleTraitPvalue[m,trt], 2))) +
                  geom_jitter(shape = 21, size = 2, colour = "black", fill = substring(m,3), width = 0.05, height = 0.2) +
                  geom_smooth(formula = y~x, method=lm, se=FALSE, linetype="dashed", color="black") +
                  xlab(m) + ylab(trt) +
                  theme_bw(15) +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5))
          )
        }
      }
      dev.off()
      
      # Filtering data to only include module trait correlations with pval < 0.05
      moduleTraitCor_f <- moduleTraitCor
      moduleTraitCor_f[moduleTraitPvalue > 0.05] <- 0
      if (all(moduleTraitCor_f == 0)) {
        cat(getTime(), "No significant module trait correlations for minimum module size of",minModuleSize,"at pval < 0.05 \n")
      } else {
        flag_row <- rowSums(moduleTraitCor_f) != 0
        flag_col <- colSums(moduleTraitCor_f) != 0
        moduleTraitCor_f <- moduleTraitCor_f[flag_row, flag_col, drop=F]
        moduleTraitPvalue_f <- moduleTraitPvalue[flag_row, flag_col, drop=F]
        textMatrix <- paste(signif(moduleTraitCor_f, 2), "\n(", signif(moduleTraitPvalue_f, 1), ")", sep = "")
        dim(textMatrix) <- dim(moduleTraitCor_f)
        yLabels <- substring(rownames(moduleTraitCor_f),3)
        ySymbols <- table(mergedColors)[yLabels]
        ySymbols <- paste(yLabels, "(" , ySymbols, ")" , sep = "")
        
        pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,"-pval_0.05.pdf"), sep = "/"), width = 14, height = 10)
        # Display the correlation values within a heatmap plot
        mar.def <- par()$mar
        par(mar = mar.def + c(0,2,0,0))
        labeledHeatmap(Matrix = moduleTraitCor_f,
                       xLabels = colnames(moduleTraitCor_f),
                       yLabels = rownames(moduleTraitCor_f),
                       ySymbols = ySymbols,
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = MTCHlabSize,
                       cex.lab.x = MTCHxSize,
                       cex.lab.y = 0.6,
                       zlim = c(-1,1),
                       main = paste("Module-trait correlations (pval < 0.05)"))
        dev.off()
        par(mar = mar.def)
        
        # Saving excel file of data used to create the module trait correlation heatmap
        moduleTraitCombined <- (function(df.corr, df.pval){
          rownames(df.corr) <- gsub("ME","", rownames(df.corr))
          rownames(df.pval) <- gsub("ME","", rownames(df.pval))
          colnames(df.corr) <- paste0(colnames(df.corr),"_CORR")
          colnames(df.pval) <- paste0(colnames(df.pval),"_PVAL")
          finaldf <- c()
          for (idx in 1:ncol(df.corr)) {
            finaldf <- cbind(finaldf, df.corr[,idx], df.pval[,idx])
            colnames(finaldf)[c(2*idx-1, 2*idx)] <- c(colnames(df.corr[,idx, drop=F]), colnames(df.pval[,idx, drop=F]))
          }
          return(finaldf)
        })(moduleTraitCor_f, moduleTraitPvalue_f)
        write.xlsx(list(Combined=moduleTraitCombined, CorrCoeff=moduleTraitCor_f, R2=moduleTraitCor_f^2, Pvalues=moduleTraitPvalue_f), paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,"-pval_0.05.xlsx"), sep = "/"), row.names = T, overwrite = T, headerStyle = createStyle(textRotation = 45), colWidths = "auto") 
      }
      
      # Filtering data to only include module trait correlations with pval < 0.05 and R2 > 5%
      moduleTraitCor_f <- moduleTraitCor
      moduleTraitCor_f[moduleTraitCor^2 < 0.05] <- 0
      moduleTraitCor_f[moduleTraitPvalue > 0.05] <- 0
      if (all(moduleTraitCor_f == 0)) {
        cat(getTime(), "No significant module trait correlations for minimum module size of",minModuleSize,"at pval < 0.05 and R^2 > 5% \n")
      } else {
        flag_row <- rowSums(moduleTraitCor_f) != 0
        flag_col <- colSums(moduleTraitCor_f) != 0
        moduleTraitCor_f <- moduleTraitCor_f[flag_row, flag_col, drop=F]
        moduleTraitPvalue_f <- moduleTraitPvalue[flag_row, flag_col, drop=F]
        textMatrix <- paste(signif(moduleTraitCor_f, 2), "\n(", signif(moduleTraitPvalue_f, 1), ")", sep = "")
        dim(textMatrix) <- dim(moduleTraitCor_f)
        yLabels <- substring(rownames(moduleTraitCor_f),3)
        ySymbols <- table(mergedColors)[yLabels]
        ySymbols <- paste(yLabels, "(" , ySymbols, ")" , sep = "")
        
        pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,"-pval_0.05_R2_0.05.pdf"), sep = "/"), width = 14, height = 10)
        # Display the correlation values within a heatmap plot
        mar.def <- par()$mar
        par(mar = mar.def + c(0,2,0,0))
        labeledHeatmap(Matrix = moduleTraitCor_f,
                       xLabels = colnames(moduleTraitCor_f),
                       yLabels = rownames(moduleTraitCor_f),
                       ySymbols = ySymbols,
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = MTCHlabSize,
                       cex.lab.x = MTCHxSize,
                       cex.lab.y = 0.6,
                       zlim = c(-1,1),
                       main = paste("Module-trait correlations (pval < 0.05, R2 > 0.05)"))
        dev.off()
        par(mar = mar.def)
        
        # Saving excel file of data used to create the module trait correlation heatmap
        moduleTraitCombined <- (function(df.corr, df.pval){
          rownames(df.corr) <- gsub("ME","", rownames(df.corr))
          rownames(df.pval) <- gsub("ME","", rownames(df.pval))
          colnames(df.corr) <- paste0(colnames(df.corr),"_CORR")
          colnames(df.pval) <- paste0(colnames(df.pval),"_PVAL")
          finaldf <- c()
          for (idx in 1:ncol(df.corr)) {
            finaldf <- cbind(finaldf, df.corr[,idx], df.pval[,idx])
            colnames(finaldf)[c(2*idx-1, 2*idx)] <- c(colnames(df.corr[,idx, drop=F]), colnames(df.pval[,idx, drop=F]))
          }
          return(finaldf)
        })(moduleTraitCor_f, moduleTraitPvalue_f)
        write.xlsx(list(Combined=moduleTraitCombined, CorrCoeff=moduleTraitCor_f, R2=moduleTraitCor_f^2, Pvalues=moduleTraitPvalue_f), paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,"-pval_0.05_R2_0.05.xlsx"), sep = "/"), row.names = T, overwrite = T, headerStyle = createStyle(textRotation = 45), colWidths = "auto") 
      }
      
      # Filtering data to only include module trait correlations with pval < 0.05 and R2 > 2.5%
      moduleTraitCor_f <- moduleTraitCor
      moduleTraitCor_f[moduleTraitCor^2 < 0.025] <- 0
      moduleTraitCor_f[moduleTraitPvalue > 0.05] <- 0
      if (all(moduleTraitCor_f == 0)) {
        cat(getTime(), "No significant module trait correlations for minimum module size of",minModuleSize,"at pval < 0.05 and R^2 > 2.5% \n")
      } else {
        flag_row <- rowSums(moduleTraitCor_f) != 0
        flag_col <- colSums(moduleTraitCor_f) != 0
        moduleTraitCor_f <- moduleTraitCor_f[flag_row, flag_col, drop=F]
        moduleTraitPvalue_f <- moduleTraitPvalue[flag_row, flag_col, drop=F]
        textMatrix <- paste(signif(moduleTraitCor_f, 2), "\n(", signif(moduleTraitPvalue_f, 1), ")", sep = "")
        dim(textMatrix) <- dim(moduleTraitCor_f)
        yLabels <- substring(rownames(moduleTraitCor_f),3)
        ySymbols <- table(mergedColors)[yLabels]
        ySymbols <- paste(yLabels, "(" , ySymbols, ")" , sep = "")
        
        pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,"-pval_0.05_R2_0.025.pdf"), sep = "/"), width = 14, height = 10)
        # Display the correlation values within a heatmap plot
        mar.def <- par()$mar
        par(mar = mar.def + c(0,2,0,0))
        labeledHeatmap(Matrix = moduleTraitCor_f,
                       xLabels = colnames(moduleTraitCor_f),
                       yLabels = rownames(moduleTraitCor_f),
                       ySymbols = ySymbols,
                       colorLabels = FALSE,
                       colors = blueWhiteRed(50),
                       textMatrix = textMatrix,
                       setStdMargins = FALSE,
                       cex.text = MTCHlabSize,
                       cex.lab.x = MTCHxSize,
                       cex.lab.y = 0.6,
                       zlim = c(-1,1),
                       main = paste("Module-trait correlations (pval < 0.05, R2 > 0.025)"))
        dev.off()
        par(mar = mar.def)
        
        # Saving excel file of data used to create the module trait correlation heatmap
        moduleTraitCombined <- (function(df.corr, df.pval){
          rownames(df.corr) <- gsub("ME","", rownames(df.corr))
          rownames(df.pval) <- gsub("ME","", rownames(df.pval))
          colnames(df.corr) <- paste0(colnames(df.corr),"_CORR")
          colnames(df.pval) <- paste0(colnames(df.pval),"_PVAL")
          finaldf <- c()
          for (idx in 1:ncol(df.corr)) {
            finaldf <- cbind(finaldf, df.corr[,idx], df.pval[,idx])
            colnames(finaldf)[c(2*idx-1, 2*idx)] <- c(colnames(df.corr[,idx, drop=F]), colnames(df.pval[,idx, drop=F]))
          }
          return(finaldf)
        })(moduleTraitCor_f, moduleTraitPvalue_f)
        write.xlsx(list(Combined=moduleTraitCombined, CorrCoeff=moduleTraitCor_f, R2=moduleTraitCor_f^2, Pvalues=moduleTraitPvalue_f), paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("module_trait_correlation_MMS",minModuleSize,"-pval_0.05_R2_0.025.xlsx"), sep = "/"), row.names = T, overwrite = T, headerStyle = createStyle(textRotation = 45), colWidths = "auto") 
      }
      
      # Step 4: Compute gene significance (GS) and module membership (MM) #---------------
      
      topn <- 10 # top n features to be labelled in the MM x GS scatterplots
      bar.topn <- 20 # number of top for bar plot
      
      # GS: Absolute correlation between individual expression variables with clinical variables
      # MM: Correlation between individual expression variables with each module eigen vector
      #GS x MM plots for significant module trait pairs
      pvalThreshold <- 0.05
      sigPairs <- apply(moduleTraitPvalue, c(1,2), FUN =  function(x) x < pvalThreshold)
      sigPairs <- data.frame(which(sigPairs, arr.ind = TRUE))
      colnames(sigPairs) <- c("module", "trait")
      modNames <- substring(names(MEs), 3)
      #MM
      geneModuleMembership <- as.data.frame(cor(d, MEs, method = corMethod , use = "p")) # correlation of cytokines with module eigen values
      MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
      names(geneModuleMembership) <- paste("MM", modNames, sep="")
      names(MMPvalue) <- paste("p.MM", modNames, sep="")
      
      if (nrow(sigPairs) != 0) {
        p.abs <- list()
        p.actual <- list()
        bar.plot <- list()
        bar.plot.data <- list()
        
        for (r in 1:nrow(sigPairs)){
          module <- modNames[sigPairs$module[r]]
          trait <- colnames(moduleTraitPvalue)[sigPairs$trait[r]]
          y <- as.data.frame(traitData[, trait]) # trait
          names(y) <- trait
          traitFlag <- rep(0, nClin)
          names(traitFlag) <- clinNames
          # GS
          geneTraitSignificance <- as.data.frame(cor(d, y, method = corMethod, use = "p")) # correlation of expr with trait, y
          GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
          names(geneTraitSignificance) <- paste("GS.", trait, sep="")
          names(GSPvalue) = paste("p.GS.", trait, sep="")
          
          column <- sigPairs$module[r] #match(module, modNames)
          moduleGenes <- moduleColors==module
          tempdata <- data.frame(geneModuleMembership[moduleGenes, ])
          
          # Bar plots of correlation of top n features within each module with each trait
          tempdata$geneTraitCor <- geneTraitSignificance[moduleGenes, 1]
          tempdata$GS <- abs(geneTraitSignificance[moduleGenes, 1])
          lastColInd <- ncol(tempdata)
          tempdata <- tempdata[, c(column, (lastColInd-1), lastColInd)]
          tempdata$MM <- abs(tempdata[,1])
          #colnames(tempdata)[1] <- "MM"
          tempdata$label <- rownames(tempdata)
          
          tempdata <- tempdata[, c(grep("MM", colnames(tempdata)),
                                   grep("geneTrait", colnames(tempdata)),
                                   grep("GS", colnames(tempdata)),
                                   grep("label", colnames(tempdata)))] #tempdata[, c(2, "geneTraitCor", "MM", "GS", "label")]
          
          # bar plots
          bardata <- tempdata[order(tempdata$GS, decreasing = TRUE),  ]
          bar.plot.data[[r]] <- bardata
          tempname <- paste0(trait,"_",module)
          if(nchar(tempname) > 27) {tempname <- substr(tempname, 1, 27)}
          names(bar.plot.data)[r] <- paste0(tempname,"_",r)
          
          if (nrow(bardata) > bar.topn) {
            bardata <- bardata[1:bar.topn,]
          }
          bar.plot[[r]] <- ggplot(bardata, aes(x=reorder(label, GS) , y=geneTraitCor)) +
            geom_bar(stat = "identity", fill = module, colour = "black", width = 0.8) + #, colour = ifelse(bardata$geneTraitCor > 0, "red", "blue")
            scale_fill_manual(values = module) +
            xlab(paste0("Top correlated features in ", module, " module")) +
            scale_y_continuous(name = paste0("correlation with trait ", trait), limits = c(-1,1)) +
            coord_flip() +
            theme_bw(20) +
            theme(axis.text.y = element_text(size = barplotYLabSize, angle = 45),
                  axis.text.x = element_text(size = 8),
                  axis.title.y = element_text(size = 11, face="bold"),
                  axis.title.x = element_text(size = 11, face="bold"))
          
          # scatter plots - MM x GS
          if (nrow(tempdata) > topn) {
            index1 <- sort(tempdata$MM, index.return = TRUE, decreasing = TRUE)
            index1 <- index1$ix[1:topn]
            index2 <- sort(tempdata$GS, index.return = TRUE, decreasing = TRUE)
            index2 <- index2$ix[1:topn]
            idxs <- union(index1, index2)
            idxs <- idxs[!is.na(idxs)]
            tempdata$label[-idxs] <- "" 
          }
          colnames(tempdata)[1] <- "MMactualCor"
          
          # MM x GS
          ## absolute correlations
          p.abs[[r]] <- ggplot(data = tempdata,
                               aes(x = MM, y = GS, label = label)) +
            geom_point(shape = 21, size = 4, colour = "black", fill = module) +
            labs(title = "MM vs GS") + #paste("module membership (MM) vs gene significance (GS)\n")#,
            #subtitle = paste("cor = ", round(cor(x,y, method = corMethod),3), "; p-val = ", round(moduleTraitPvalue[m,trt], 2))
            #) +
            geom_smooth(formula = y~x, method=lm, se=FALSE, linetype="dashed", color="black") +
            xlab(paste("Abs cor of features with", module, "module (MM)")) +
            ylab(paste("Abs cor of features with", trait, "(GS)")) +
            geom_text_repel(size = 3, max.overlaps = Inf) +
            theme_bw(25) +
            theme(plot.title = element_text(hjust = 0.5)#,
                  #plot.subtitle = element_text(hjust = 0.5)
            )
          
          ## actual correlations
          p.actual[[r]] <- ggplot(data = tempdata, aes(x = MMactualCor, y = geneTraitCor, label = label)) +
            geom_point(shape = 21, size = 4, colour = "black", fill = module) +
            labs(title = "MM vs GS") +
            geom_smooth(formula = y~x, method=lm, se=FALSE, linetype="dashed", color="black") +
            xlab(paste("Cor of features with", module, "module")) +
            ylab(paste("Cor of features with", trait)) +
            geom_text_repel(size = 3, max.overlaps = Inf) +
            theme_bw(25) +
            theme(plot.title = element_text(hjust = 0.5))
          
          # csv file with GS and MM info for each trait
          if (!traitFlag[trait]){
            # Create the starting data frame
            geneInfo0 <- cbind(exprId = colnames(d), geneTraitSignificance, GSPvalue, moduleColors)
            # Order modules by their significance for y
            modOrder <- order(-abs(cor(MEs, y, method = corMethod, use = "p")))
            # Add module membership information in the chosen order
            for (mod in 1:ncol(geneModuleMembership))
            {
              oldNames = names(geneInfo0)
              geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                                     MMPvalue[, modOrder[mod]])
              names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                                   paste("p.MM.", modNames[modOrder[mod]], sep=""))
            }
            # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
            geneOrder <- order(geneInfo0$moduleColors, -abs(geneInfo0[, 2]))
            geneInfo <- geneInfo0[geneOrder, ]
            
            write.csv(geneInfo, file = paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("GS_MM_Info_", trait, "_MMS", minModuleSize, ".csv"), sep = "/"), row.names = F)
          } # end of traitFlag
        } # end of for in r sig.pairs
        
        pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("barplots_topCorrFeatures_MMS",minModuleSize,".pdf"), sep = "/"))
        for (i in seq_along(bar.plot)) {
          print(bar.plot[[i]]) 
        }
        dev.off()
        
        pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("scatterplots_MMvsGS_absCorr_MMS",minModuleSize,".pdf"), sep = "/")) # for 2-pnel, , width = 12
        for (i in seq_along(p.abs)) {
          print(p.abs[[i]])
        }
        dev.off()
        
        pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("scatterplots_MMvsGS_corr_MMS",minModuleSize,".pdf"), sep = "/")) # for 2-pnel, , width = 12
        for (i in seq_along(p.actual)) {
          print(p.actual[[i]])
        }
        dev.off()
        
        write.xlsx(bar.plot.data, file = paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("barplots_topCorrFeatures_MMS",minModuleSize,".xlsx"), sep = "/"), overwrite = TRUE, rowNames = TRUE, colWidths = "auto")
      }
      
      # Step 5: correlation heatmap #---------------
      
      heatSortOrder <- order(match(moduleColors, substring(colnames(MEs),3)))
      corMat <- as.matrix(cor(d, method = corMethod, use = "p"))
      col_list <- unique(moduleColors)
      names(col_list) <- col_list
      
      pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("correlation_heatmap_MMS",minModuleSize,".pdf"), sep = "/"))
      # correlation heatmap
      ht <- Heatmap(corMat,
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    show_row_names = FALSE, show_column_names = FALSE,
                    row_order = heatSortOrder, column_order = heatSortOrder,
                    column_title = "correlation heatmap",
                    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                    col = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red")),
                    top_annotation = HeatmapAnnotation(module = moduleColors, col = list(module = col_list)),
                    left_annotation = HeatmapAnnotation(which = "row", module = moduleColors, col = list(module = col_list), show_legend = FALSE),
                    show_heatmap_legend = TRUE,
                    name = "corr")
      draw(ht)
      dev.off()
      
      pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("correlation_heatmap_MMS",minModuleSize,"-abs.pdf"), sep = "/"))
      # abs(correlation) heatmap
      ht <- Heatmap(abs(corMat),
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    show_row_names = FALSE, show_column_names = FALSE,
                    row_order = heatSortOrder, column_order = heatSortOrder,
                    column_title = "absolute correlation heatmap",
                    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                    col = colorRamp2(c(0,1), c("white", "red")),
                    top_annotation = HeatmapAnnotation(module = moduleColors, col = list(module = col_list)),
                    left_annotation = HeatmapAnnotation(which = "row", module = moduleColors, col = list(module = col_list), show_legend = FALSE),
                    show_heatmap_legend = TRUE,
                    name = "corr")
      draw(ht)
      dev.off()
      
      pdf(paste(paste0(folder,"/minModuleSize",minModuleSize), paste0("correlation_heatmap_MMS",minModuleSize,"-sq.pdf"), sep = "/"))
      # correlation^2 heatmap
      ht <- Heatmap(corMat^2,
                    cluster_rows = FALSE, cluster_columns = FALSE,
                    show_row_names = FALSE, show_column_names = FALSE,
                    row_order = heatSortOrder, column_order = heatSortOrder,
                    column_title = "squared correlation heatmap",
                    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
                    col = colorRamp2(c(0,1), c("white", "red")),
                    top_annotation = HeatmapAnnotation(module = moduleColors, col = list(module = col_list)),
                    left_annotation = HeatmapAnnotation(which = "row", module = moduleColors, col = list(module = col_list), show_legend = FALSE),
                    show_heatmap_legend = TRUE,
                    name = "corr")
      draw(ht)
      dev.off()

      cat(getTime(), "WGCNA for minimum module size of",minModuleSize,"complete \n")
    }, error=function(e){cat(getTime(), "ERROR at MMS =", minModuleSize, ":", conditionMessage(e), "Only grey module found \n")})
  }
  sink()
}
