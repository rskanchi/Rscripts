# association between outcome and predictor data, both categorical
# chi sq test of association is performed 
# bar plots with counts and proportions generated

get.association.barplot <- function(xydata, x, y, pval, theme_size =15, output_folder = NULL){
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(EnvStats)
  library(gridExtra)
  library(cowplot)
  xydata <- na.omit(xydata)
  #nY <- table(xydata[, y])
  nXY <- table(xydata[,x], xydata[,y])
  rownames(nXY) <- paste(x, "=", rownames(nXY))
  colnames(nXY) <- paste(y, "=", colnames(nXY))
  nXY <- addmargins(nXY)
  
  p1 <- ggplot(xydata, aes(x = xydata[, y], fill = xydata[, x])) +
    geom_bar(width = 0.4) + #show.legend = FALSE, 
    xlab(y) + ylab(pvar) + labs(title = "count", fill = x) +
    theme_bw(theme_size)
  
  p2 <- ggplot(xydata, aes(x = xydata[, y], fill = xydata[, x])) +
    geom_bar(position = "fill", width = 0.4) + #
    xlab(y) + ylab(" ") + labs(title = "proportion", fill = x) + 
    theme_bw(theme_size)
  
  legend <- getPlotLegend(p1)
  ggsave(file = paste0(output_folder, "/chisq_association_barplots_", ovar, ".pdf"), width = 10, height = 6, units = "in",
         grid.arrange(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), 
                      tableGrob(nXY), legend,
                      top = paste("chisq association p-value =", scientific(pval)),
                      nrow = 2, heights = c(4,1))
  )
} # end of function get.association.barplot

getChiSq.association <- function(data, xNames, yNames, output_folder, ...){
  library(ggplot2)
  library(dplyr)
  createDir(output_folder)
  # make the variables factor if not already
  data <- data %>% mutate_at(c(xNames, yNames), as.factor)
  for(ovar in yNames){
    # chi sq test of association between the 
    dfChi <- t(sapply(xNames, FUN = function(x){
      noNA.data <- na.omit(data[, c(ovar, x)]) 
      noNA.data <- noNA.data %>% mutate_at(c(ovar, x), as.factor)
      chi <- chisq.test(noNA.data[, ovar], noNA.data[, x])
      c(chi$statistic, chi$parameter, chi$p.value)
    }))
    colnames(dfChi) <- c("chisq.statistic", "df", "pvalue")
    
    assoc.barplot <- lapply(xNames, FUN = function(pvar){get.association.barplot(xydata = data, x = pvar, y = ovar, pval = dfChi[pvar, "pvalue"], output_folder = output_folder)})
    write.csv(dfChi, file = paste0(output_folder, "/chisq_association_", ovar, ".csv"))
  } # end of for loop in ovar
} # end of function getChiSq.association

# function to extract legend from plot 
getPlotLegend <- function(plot) { 
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  legend <- plot_table$grobs[[legend_plot]] 
  return(legend) 
}
