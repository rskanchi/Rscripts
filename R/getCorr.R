# corrFunctions.R
library(ggplot2)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(ggcorrplot)

# function to compute correlation, p values and adjusted pvalues 
# two sets of variables could be specified
# one set used as rows - xNames (must be specified)
# another as columns - yNames (could be null)
# padj.n = 
# "self" when self correlation (one matrix) | n*(n-1)/2 number of p values where n is the number of xNames
# "columnwise" when you would like to compute padj column-wise (y matrix column-wise) | number of rows of d for each column
# "overall" when computing padj based on p values from two different matrices | n(xNames) * n(yNames)
compute.corr <- function(d, xNames, yNames = NULL, corMethod = "spearman", padjMethod = "BH", folder = NULL, padj.n = "overall"){ 
  library(openxlsx)
  if (is.null(yNames)) {
    d <- d[, xNames]
    yNames <- xNames
    #if (padj.n == "self") pn <- ncol(d)*(ncol(n)-1)/2
  } else  {
    d <- d[, c(xNames, yNames)]
    #if (padj.n == "overall") pn <- nrow(d) * ncol(d)
    #else if (padj.n == "columnwise") pn <- nrow(d)
  } # end of else
  
  nX <- length(xNames)
  nY <- length(yNames)
  dimNames <- list(yNames, xNames)
  cMat <- matrix(ncol = nX, nrow = nY, dimnames = dimNames)
  pMat <- matrix(ncol = nX, nrow = nY, dimnames = dimNames)
  padjMat <- matrix(ncol = nX, nrow = nY, dimnames = dimNames)
  
  for (c in 1:nX){
    for (r in 1:nY){
      single <- cor.test(d[, xNames[c]], d[, yNames[r]], use = "pairwise",
                         method = corMethod, exact = FALSE, alternative = "two.sided")
      cMat[r,c] <- single$estimate
      pMat[r,c] <- single$p.value
    } #end row
  } # end col
  #padjMat <- apply(pMat, 2, FUN = function(pcol) p.adjust(pcol, method = padjMethod))
  
  # number of p values based on which to compute the padj
  pn <- switch(padj.n,
               "overall" = nrow(cMat)*ncol(cMat),
               "columnwise" = nrow(cMat),
               "self" = (nrow(cMat)-1)*nrow(cMat)/2
  ) # end of switch
  
  padjMat <- apply(pMat, c(1,2), FUN = function(px) p.adjust(px, method = padjMethod, n = pn))
  
  if(!is.null(folder)){ if (!dir.exists(folder)){ dir.create(folder, recursive = TRUE) } } 
  openxlsx::write.xlsx(list("corr" = cMat, "pvalues" = pMat, "padjmat" = padjMat, "padj.n" = c(padj.n, pn)), 
                       file = paste(folder, paste("OutputTables_corr_", corMethod, ".xlsx", sep = ""), sep = "/"), 
                       rowNames = TRUE, overwrite = TRUE,
                       colWidth = "auto",
                       firstRow = TRUE, firstCol = TRUE,
                       headerStyle = createStyle(textDecoration = "BOLD", textRotation = 45)
  )
  return(list("corr" = cMat , "p" = pMat, "padj" = padjMat, "p.adj_type" = padj.n, "p.adj_n" = pn))
} # end of function compute.corr

# get corr Heatmap
getCorrHeatmap <- function(res, padjThreshold = 1, fontsize = 15, xTitle = "", yTitle = "", folder,
                           top_bar = TRUE, right_bar = TRUE, pdf.width = 10, pdf.height = 10,
                           ht_row_fontsize = 10, ht_col_fontsize = 10,
                           corMethod = "spearman"){
  
  col_fun <- colorRamp2(c(-1,0,1), c("blue", "white", "red"))
  xNames <- colnames(res$corr)
  yNames <- rownames(res$corr)
  
  if (padjThreshold < 1){
    indices <- res$padj >= padjThreshold
    res$corr[indices] <- 0
    
    # remove rows and cols with all 0 corr
    rowcorSum <- apply(indices, 1, FUN = function(x) sum(!x))
    rowKeep <- which(rowcorSum > 0)
    colcorSum <- apply(indices, 2, FUN = function(x) sum(!x))
    colKeep <- which(colcorSum > 0)
    
    res$corr <- res$corr[rowKeep, colKeep]
    res$p <- res$p[rowKeep, colKeep]
    res$padj <- res$padj[rowKeep, colKeep]
    
    heat_col_title <- paste(xTitle, ", n = ", length(colKeep), " | ", length(xNames), sep = "")
    heat_row_title <- paste(yTitle, ", n = ", length(rowKeep), " | ", length(yNames), sep = "")
  } else {
    heat_col_title <- paste(xTitle, ", n = ", length(xNames), sep = "")
    heat_row_title <- paste(yTitle, ", n = ", length(yNames), sep = "")
  } # end of if in padjThreshold
  
  
  # top annotation barplot
  posCor <- apply(as.matrix(res$corr), 2, FUN = function(x) sum(x > 0))
  negCor <- apply(as.matrix(res$corr), 2, FUN = function(x) sum(x < 0))
  if (top_bar){     ha <- HeatmapAnnotation(freq = anno_barplot(cbind(negCor, posCor), 
                                                                gp = gpar(fill = c("blue", "red"), col = "white"), #
                                                                border = FALSE, height = unit(3, "cm")), 
                                            show_annotation_name = FALSE)}
  
  gapanno <- rowAnnotation(emp = anno_empty(border = FALSE), show_annotation_name = FALSE)
  gapanno2 <- HeatmapAnnotation(emp = anno_empty(border = FALSE), show_annotation_name = FALSE)
  
  # right annotation barplot
  posCor_rt <- apply(as.matrix(res$corr), 1, FUN = function(x) sum(x > 0))
  negCor_rt <- apply(as.matrix(res$corr), 1, FUN = function(x) sum(x < 0))
  if (right_bar){     ha_rt <- rowAnnotation(freq = anno_barplot(cbind(negCor_rt, posCor_rt), 
                                                                 gp = gpar(fill = c("blue", "red"), col = "white"), #
                                                                 border = FALSE, width = unit(3, "cm")), #height = unit(3, "cm")), 
                                             show_annotation_name = FALSE)}
  
  # heatmap
  ht <- draw(
    gapanno +  
      Heatmap(res$corr, name = corMethod,
              col = col_fun, # color
              show_row_names = TRUE, cluster_rows = TRUE, #row_names_side = "left",
              row_names_gp = gpar(fontsize = ht_row_fontsize),
              show_column_names = TRUE, cluster_columns = TRUE, 
              column_names_gp = gpar(fontsize = ht_col_fontsize),  
              top_annotation = if (top_bar) c(ha, gapanno2) else NULL, 
              column_title = heat_col_title, column_title_side = "bottom",
              column_title_gp = gpar(fontsize = fontsize), #column_names_max_height = unit(3, "cm"),
              row_title = heat_row_title, 
              row_title_gp = gpar(fontsize = fontsize),
              right_annotation = if (right_bar) ha_rt else NULL,
              border = TRUE
      ))
  
  xNamesOrder <- xNames[column_order(ht)]
  yNamesOrder <- yNames[row_order(ht)]
  
  pdf(paste(folder, paste("Heatmap_corr_", corMethod, "_padj_", padjThreshold, ".pdf", sep = ""), sep = "/"), 
      width = pdf.width, height = pdf.height) #
  print(ht)
  dev.off()
  
  openxlsx::write.xlsx(list("corr" = res$corr, "pvalue" = res$p, "padj" = res$padj, 
                            "xNamesOrder" = xNamesOrder, "yNamesOrder" = yNamesOrder,
                            "freq_topBar" = data.frame(cbind(negCor, posCor)), 
                            "freq_rt_bar" = data.frame(cbind(negCor_rt, posCor_rt))
  ),
  rowNames = TRUE,
  file = paste(folder, 
               paste("SupportingData_Heatmap_padj_", padjThreshold,".xlsx", sep = ""), 
               sep = "/")
  )
  
} # end of function getCorrHeatmap


# get corr Heatmap based on p-values
getCorrHeatmap_p <- function(res, pvalue = 0.05, fontsize = 15, xTitle = "", yTitle = "", folder,
                             top_bar = TRUE, right_bar = TRUE, pdf.width = 10, pdf.height = 10,
                             ht_row_fontsize = 10, ht_col_fontsize = 10,
                             corMethod = "spearman"){
  
  col_fun <- colorRamp2(c(-1,0,1), c("blue", "white", "red"))
  xNames <- colnames(res$corr)
  yNames <- rownames(res$corr)
  
  indices <- res$p >= pvalue
  res$corr[indices] <- 0
  
  # remove rows and cols with all 0 corr
  rowcorSum <- apply(indices, 1, FUN = function(x) sum(!x))
  rowKeep <- which(rowcorSum > 0)
  colcorSum <- apply(indices, 2, FUN = function(x) sum(!x))
  colKeep <- which(colcorSum > 0)
  
  res$corr <- res$corr[rowKeep, colKeep]
  res$p <- res$p[rowKeep, colKeep]
  res$padj <- res$padj[rowKeep, colKeep]
  
  heat_col_title <- paste(xTitle, ", n = ", length(colKeep), " | ", length(xNames), sep = "")
  heat_row_title <- paste(yTitle, ", n = ", length(rowKeep), " | ", length(yNames), sep = "")
  
  # top annotation barplot
  posCor <- apply(as.matrix(res$corr), 2, FUN = function(x) sum(x > 0))
  negCor <- apply(as.matrix(res$corr), 2, FUN = function(x) sum(x < 0))
  if (top_bar){     ha <- HeatmapAnnotation(freq = anno_barplot(cbind(negCor, posCor), 
                                                                gp = gpar(fill = c("blue", "red"), col = "white"), #
                                                                border = FALSE, height = unit(3, "cm")), 
                                            show_annotation_name = FALSE)}
  
  gapanno <- rowAnnotation(emp = anno_empty(border = FALSE), show_annotation_name = FALSE)
  gapanno2 <- HeatmapAnnotation(emp = anno_empty(border = FALSE), show_annotation_name = FALSE)
  
  # right annotation barplot
  posCor_rt <- apply(as.matrix(res$corr), 1, FUN = function(x) sum(x > 0))
  negCor_rt <- apply(as.matrix(res$corr), 1, FUN = function(x) sum(x < 0))
  if (right_bar){     ha_rt <- rowAnnotation(freq = anno_barplot(cbind(negCor_rt, posCor_rt), 
                                                                 gp = gpar(fill = c("blue", "red"), col = "white"), #
                                                                 border = FALSE, width = unit(3, "cm")), #height = unit(3, "cm")), 
                                             show_annotation_name = FALSE)}
  
  # heatmap
  ht <- draw(
    gapanno +  
      Heatmap(res$corr, name = corMethod,
              col = col_fun, # color
              show_row_names = TRUE, cluster_rows = TRUE, #row_names_side = "left",
              row_names_gp = gpar(fontsize = ht_row_fontsize),
              show_column_names = TRUE, cluster_columns = TRUE, 
              column_names_gp = gpar(fontsize = ht_col_fontsize),  
              top_annotation = if (top_bar) c(ha, gapanno2) else NULL, 
              column_title = heat_col_title, column_title_side = "bottom",
              column_title_gp = gpar(fontsize = fontsize), #column_names_max_height = unit(3, "cm"),
              row_title = heat_row_title, 
              row_title_gp = gpar(fontsize = fontsize),
              right_annotation = if (right_bar) ha_rt else NULL,
              border = TRUE
      ))
  
  xNamesOrder <- xNames[column_order(ht)]
  yNamesOrder <- yNames[row_order(ht)]
  
  pdf(paste(folder, paste("Heatmap_corr_", corMethod, "_p_", pvalue, ".pdf", sep = ""), sep = "/"), 
      width = pdf.width, height = pdf.height) #
  print(ht)
  dev.off()
  
  openxlsx::write.xlsx(list("corr" = res$corr, "pvalue" = res$p, "padj" = res$padj, 
                            "xNamesOrder" = xNamesOrder, "yNamesOrder" = yNamesOrder,
                            "freq_topBar" = data.frame(cbind(negCor, posCor)), 
                            "freq_rt_bar" = data.frame(cbind(negCor_rt, posCor_rt))
  ),
  rowNames = TRUE,
  file = paste(folder, 
               paste("SupportingData_Heatmap_p_", pvalue,".xlsx", sep = ""), 
               sep = "/")
  )
  
} # end of function getCorrHeatmap_p


