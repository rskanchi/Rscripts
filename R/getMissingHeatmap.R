# function to plot a heatmap of missing values for a data frame


getMissingHeatmap <- function(data, folder){
  
  na.data <- apply(data, MARGIN = c(1,2), FUN = function(x) if (is.na(x)) 1 else 0)
  library(ComplexHeatmap)
  library(circlize)
  # top annotation barplot
  miss.n <- apply(na.data, 2, FUN = function(x) sum(x))
  miss.pct <- apply(na.data, 2, FUN = function(x) round(sum(x)*100/length(x),1))
  ha.n <- HeatmapAnnotation(freq = anno_barplot(miss.n,gp = gpar(fill = "grey", col = "white"), border = FALSE, 
                                                height = unit(2, "cm")))#, 
  #                          show_annotation_name = FALSE)
  gapanno <- HeatmapAnnotation(emp = anno_empty(border = FALSE, height = unit(0.2, "cm")), show_annotation_name = FALSE)
  
  ha.pct <- HeatmapAnnotation(pct = anno_barplot(miss.pct, 
                                                 gp = gpar(fill = "grey", col = "white"), 
                                                 border = FALSE, height = unit(1, "cm")))#, 
  # show_annotation_name = FALSE)
  htMap <-   Heatmap(na.data, name = "missing",
                     col = colorRamp2(c(0,1), c("white", "red")), # color
                     show_row_names = F, cluster_rows = TRUE, 
                     row_names_gp = gpar(fontsize = 3),
                     show_column_names = TRUE, cluster_columns = TRUE, 
                     column_names_gp = gpar(fontsize = 10),  
                     column_title_gp = gpar(fontsize = 10), 
                     row_title_gp = gpar(fontsize = 10),
                     top_annotation = c(ha.n, gapanno, ha.pct, gapanno),
                     border = TRUE
  )
  
  pdf(paste(folder, "missing_heatmap.pdf", sep = "/"))
  draw(htMap)
  dev.off()
  
  write.csv(cbind(miss.n, miss.pct), paste(folder, "missing_summary.csv", sep = "/"))

} # end of function getMissingHeatmap


