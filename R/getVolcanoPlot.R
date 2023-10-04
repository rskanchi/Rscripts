# get volcano plot; the input data is an output from the differential expression analysis
# xName = variable name for log FC
# yName = variable name for either pvalue or adjusted pvalue
# rowLabels = variable name for the col containing genes, or the features

getVolcanoPlot <- function(data, xName, yName, featureCol, folder = NULL,
                           pcutoffs = c(0.25, 0.1, 0.05), # 
                           titleText = NULL,
                           xLabel = "logFC", yLabel = "FDR", featureLabel = "Genes"){ # yLabel = "P.Value" or "adj P.Value"
  library(stringi)
  library(ggplot2)
  
  data <- data[, c(featureCol, xName, yName)]
  colnames(data) <- c("name", "x", "y")
  
  # volcano plot and output table
  deColors <- c("blue", "red", "black")
  names(deColors) <- c("DOWN", "UP", "NO")
  
  # DE columns
  DElist <- lapply(pcutoffs, FUN = function(p) {
    de <- rep("NO", nrow(data))
    de[data$x > 0 & data$y < p] <- "UP"
    de[data$x < 0 & data$y < p] <- "DOWN"
    return(de)
  })
  DElist <- as.data.frame(stri_list2matrix(DElist))
  colnames(DElist) <- paste0("DE_", pcutoffs)
  
  DEtable <- apply(DElist, 2, table)
  DEtable <- DEtable[c("UP", "DOWN", "NO"),]
  
  data <- data.frame(data, DElist)
  xlim <- c(-floor(max(abs(data$x), na.rm = TRUE)+1), floor(max(abs(data$x), na.rm = TRUE)+1))
  
  volcano_plots <- lapply(pcutoffs, FUN = function(p){
    tempData <- data[, c("name", "x", "y", paste0("DE_", p))]
    colnames(tempData)[4] <- "DE"
    tempPlot <- ggplot(data = tempData, aes(x = x, y = -log10(y), col = DE)) + #, label = delabel
      xlim(xlim = xlim) +
      geom_point(size = 1) +
      scale_color_manual(values = deColors) +
      #geom_text_repel(size = text_size, max.overlaps = Inf) +
      #geom_hline(yintercept = -log10(thresholdAdjPVAL), lty = linetype) +
      xlab("log2 fold change") +
      ylab(paste0("-log10(", yLabel, ")")) +
      labs(title = titleText,
           subtitle = paste0(paste("Total =", nrow(tempData),
                            "| UP =", sum(tempData$DE == "UP", na.rm = TRUE),
                            "| DOWN =", sum(tempData$DE == "DOWN", na.rm = TRUE),
                            sep = " "), "\n", paste0(yLabel, " threshold = ", p))) +
      theme_bw(15)
    
    return(tempPlot)
  })
  
  names(volcano_plots) <- names(DElist)

  if (!is.null(folder)){
    if (is.null(titleText)) titleText <- ""
    createDir(paste(folder, "volcano", sep = "/"))
    pdf(file = paste(folder, "volcano", paste0("volcanoplot_", yLabel, "_", featureLabel, "_", titleText, ".pdf"), sep = "/"))
    print(volcano_plots)
    dev.off()
  }
  
  write.csv(DEtable, file = paste(folder, "volcano", paste0("DEtable_", yLabel, "_", featureLabel, "_", titleText, ".csv"), sep = "/"))
} # end of getVolcaoPlot

