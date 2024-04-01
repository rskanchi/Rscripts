library(edgeR)
library(tidyverse)
library(RColorBrewer)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)

performRNAseqUQnormalization <- function(work_dir, files){
  setwd(work_dir) #setwd("/Users/coarfa/Research/IvanRosas/Loor/Exosome_Redux/")
  # files <-c("EVLP_combined_combatseq_counts.06_20_2023.csv")
  
  dlist <- map(files, read_csv)
  names(dlist)<-files %>% basename %>% str_replace(".csv",".UQ.csv")
  
  ylist<-map(dlist,~DGEList(counts = .x[,-1:-3],genes = .x[,1:3]))
  ylist <- map(ylist,~UQfilter(x = .x,samples = 16))
  ylist <- map(ylist,~calcNormFactors(object = .x,method = "upperquartile"))
  
  map2(.x = ylist,.y = c("protein.coding"),~UQheatmap(results = .x,label = "z-score log2CPM",column_cluster = T,row_cluster = T,nam = .y))
  map2(.x = ylist,.y = str_c("CPM.",names(ylist)),~write_csv(x = cbind(.x$genes,cpm(y = .x)),path = .y))
} # end of function performRNAseqUQnormalization

UQheatmap <-function(results,label,column_cluster,row_cluster,nam){
  results <- cbind(results$genes,cpm(results,log = T))
  b <- results[,-1:-3]
  b <- t(scale(t(b)))
  df_scale<-b
  
  if(!near(dim(df_scale)[1],0)){
    
    row.names(df_scale)<-results$GeneSymbol
    
    legend <- label
    x<-abs(max(df_scale))
    y<-abs(min(df_scale))
    scale<-max(c(x,y))
    
    col_heatmap <- colorRamp2(c((-scale)*.75, 0,(scale)*.75), c("BLUE2", "white", "RED"))
    xlsx.tbl<-df_scale%>%as.data.frame
    print(head(xlsx.tbl))
    print(ncol(xlsx.tbl))
    print("table")
    
    wb <- createWorkbook()
    addWorksheet(wb = wb,sheetName = "heatmap")
    writeData(wb = wb,sheet = "heatmap",x = xlsx.tbl,rowNames = T)
    conditionalFormatting(wb = wb,sheet = "heatmap",cols = 1:(ncol(xlsx.tbl)+1),rows =1:(nrow(xlsx.tbl)+1),style = c("blue","white","red"),rule = c(-2,0,2),type = "colourScale" )
    saveWorkbook(wb = wb,file =str_c(nam,".xlsx"),overwrite = T )
  } # end of if !near
} # end of heatmap

UQfilter <- function(x, samples){
  keep <- rowSums(cpm(y = x, log = F) > 10) > samples
  print(summary(keep))
  return(x[keep, , keep.lib.sizes = FALSE])
} # end of UQfilter
