# get jaccard heatmaps
# list 1 and 2 are the lists of elements; data frame columns as list elements

getJaccardHeatmaps <- function(list1, list2, 
                               dataname1, dataname2, # labels for identifying list 1 and 2 respectively
                               fontsize = 14, ht_row_fontsize = 10, ht_col_fontsize = 10,
                               top_bar = TRUE, right_bar = TRUE,
                               heat_row_title = NULL, heat_col_title = NULL,
                               pdf.width = 12, pdf.height = 12,
                               folder){
  
  lenList1 <- length(list1)
  lenList2 <- length(list2)
  namesList1 <- names(list1)
  namesList2 <- names(list2)
  
  nIntersection <- matrix(nrow = lenList1, ncol = lenList2, dimnames = list(namesList1, namesList2))
  nUnion <- matrix(nrow = lenList1, ncol = lenList2, dimnames = list(namesList1, namesList2))
  jSimilarity <- matrix(nrow = lenList1, ncol = lenList2, dimnames = list(namesList1, namesList2))
  jDistance <- matrix(nrow = lenList1, ncol = lenList2, dimnames = list(namesList1, namesList2))
  
  for (lrow in 1:lenList1){
    for (lcol in 1:lenList2){
      jTemp <- getJaccard(list1[[namesList1[lrow]]], list2[[namesList2[lcol]]])
      nIntersection[lrow, lcol] <- jTemp$nIntersection
      nUnion[lrow, lcol] <- jTemp$nUnion
      jSimilarity[lrow, lcol] <- jTemp$similarity
      jDistance[lrow, lcol] <- jTemp$distance
    } # end of for in l2
  } # end of for in l1
  
  col_fun <- colorRamp2(c(0,1), c("white", "red"))
  # top annotation barplot
  nGenesTop <- lengths(list2)
  if (top_bar){     ha <- HeatmapAnnotation(freq = anno_barplot(nGenesTop, 
                                                                gp = gpar(fill = "black", col = "white"), #
                                                                border = FALSE, height = unit(3, "cm")), 
                                            show_annotation_name = FALSE)}
  
  #gapanno <- rowAnnotation(emp = anno_empty(border = FALSE), show_annotation_name = FALSE)
  #gapanno2 <- HeatmapAnnotation(emp = anno_empty(border = FALSE), show_annotation_name = FALSE)
  
  # right annotation barplot
  nGenesRight <- lengths(list1)
  if (right_bar){     ha_rt <- rowAnnotation(freq = anno_barplot(nGenesRight, 
                                                                 gp = gpar(fill = "black", col = "white"), #
                                                                 border = FALSE, width = unit(3, "cm")), #height = unit(3, "cm")), 
                                             show_annotation_name = FALSE)}
  
  
  library(dendsort)
  row_dend <- dendsort(hclust(dist(jSimilarity)))
  col_dend <- dendsort(hclust(dist(t(jSimilarity))))
  #  top_annotation <- HeatmapAnnotation(freq = anno_boxplot(jSimilarity, height = unit(3, "cm")))
  
  # the heatmap has the number of intersection elements printed in the cells when intersection > 0
  ht <- draw(
    #gapanno +  
    Heatmap(jSimilarity, name = "Jaccard Similarity",
            col = col_fun, # color
            show_row_names = TRUE, cluster_rows = row_dend, #TRUE,
            row_names_gp = gpar(fontsize = ht_row_fontsize),
            show_column_names = TRUE, cluster_columns = col_dend, #TRUE, 
            column_names_gp = gpar(fontsize = ht_col_fontsize),  
            top_annotation = if (top_bar) ha else NULL, #c(ha, gapanno2)
            column_title = heat_col_title, column_title_side = "bottom",
            column_title_gp = gpar(fontsize = fontsize), #column_names_max_height = unit(3, "cm"),
            row_title = heat_row_title, 
            row_title_gp = gpar(fontsize = fontsize),
            right_annotation = if (right_bar) ha_rt else NULL,
            cell_fun = function(j, i, x, y, width, height, fill) {
              if (nIntersection[i,j] > 0)
                grid.text(sprintf("%.0f",nIntersection[i,j]), x, y, gp = gpar(fontsize = 5))
            },
            border = TRUE
    ))
  
  xNamesOrder <- namesList1[column_order(ht)]
  yNamesOrder <- namesList2[row_order(ht)]
  
  # the heatmap has the similarity printed in the cells when similarity > 0
  ht2 <- draw(
    #gapanno +  
    Heatmap(jSimilarity, name = "Jaccard Similarity",
            col = col_fun, # color
            show_row_names = TRUE, cluster_rows = row_dend, #TRUE,
            row_names_gp = gpar(fontsize = ht_row_fontsize),
            show_column_names = TRUE, cluster_columns = col_dend, #TRUE, 
            column_names_gp = gpar(fontsize = ht_col_fontsize),  
            top_annotation = if (top_bar) ha else NULL, #c(ha, gapanno2)
            column_title = heat_col_title, column_title_side = "bottom",
            column_title_gp = gpar(fontsize = fontsize), #column_names_max_height = unit(3, "cm"),
            row_title = heat_row_title, 
            row_title_gp = gpar(fontsize = fontsize),
            right_annotation = if (right_bar) ha_rt else NULL,
            cell_fun = function(j, i, x, y, width, height, fill) {
              if (jSimilarity[i,j] > 0)
                grid.text(sprintf("%.2f",jSimilarity[i,j]), x, y, gp = gpar(fontsize = 5))
            },
            border = TRUE
    ))
  
  folder <- paste(folder, paste("Jaccard_Similarity", dataname1, dataname2, sep = "_"), sep = "/")
  createDir(folder)
  pdf(paste(folder, paste("Heatmap_Jaccard_Similarity_", dataname1, "_", dataname2, ".pdf", sep = ""), sep = "/"), 
      width = pdf.width, height = pdf.height) #
  print(ht)
  print(ht2)
  dev.off()
  
  openxlsx::write.xlsx(list("nIntersection" = nIntersection, "nUnion" = nUnion, "jSimilarity" = jSimilarity, 
                            "xNamesOrder" = xNamesOrder, "yNamesOrder" = yNamesOrder,
                            "freq_topBar" = data.frame(nGenesTop), 
                            "freq_rt_bar" = data.frame(nGenesRight)
  ),
  rowNames = TRUE,
  file = paste(folder, 
               paste("SupportingData_Jaccard_Similarity_", dataname1, "_", dataname2,".xlsx", sep = ""), 
               sep = "/")
  )
} # end of function getJaccardHeatmaps
