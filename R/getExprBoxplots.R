# new after the getTransformedData function was written
# d is the expression data with sample rows and expression columns 
getExprBoxplots <- function(d, dataTransformation = NULL, nVarexpr = NULL,
                            outlier.size = 0.2, outlier.color = "red", box.outline.color = "grey",
                            all_expr_boxplot = TRUE,
                            folder = NULL){
  
  if (is.null(folder)) folder <- paste(getwd(), "output", "exprBoxplots", sep = "/") else
    folder <- paste(getwd(), "output", "exprBoxplots", sep = "/")
  createDir(folder)
  
  getFunctionMessage(text = paste("output path:", folder), func = "getExprBoxplots")
  
  # variability measure: sd
  stdev <- apply(d, 2, sd)
  d <- d[, names(stdev)[rev(order(stdev))]] # data sorted by stdev
  
  # full data boxplots
  if (all_expr_boxplot){
    # boxplot expr data
    boxplot_fulldata <- ggplot(stack(data.frame(d)), aes(x = ind, y = values)) + 
      geom_boxplot(color = box.outline.color, outlier.color = outlier.color, outlier.size = outlier.size) + 
      xlab(paste("#expr = ", ncol(d))) + ylab("expression data") +
      theme_classic() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) 
    
    jpeg(filename = paste(folder, "boxplots_exprData_full.jpg", sep = "/"))
    boxplot_fulldata
    dev.off()
    
    pdf(paste(folder, paste("boxplots_exprData_full.pdf", sep = ""), sep = "/"))
    print(boxplot_fulldata) # end of print
    dev.off()
  } # end of if all_expr_boxplot = TRUE
  
  
  if (!is.null(nVarexpr)){
    if (is.positive.integer(nVarexpr)){
      d.nVarexpr <- d[, 1:nVarexpr]
      
      boxplot_nVarexpr <- ggplot(stack(data.frame(d.nVarexpr)), aes(x = ind, y = values)) + 
        geom_boxplot(color = box.outline.color, outlier.color = outlier.color, outlier.size = outlier.size) + 
        xlab(paste("#expr = ", ncol(d.nVarexpr))) + ylab("expression data") +
        theme_classic() +
        theme(axis.ticks.x = element_blank(),
              axis.text.x = element_blank()) 
      
      jpeg(filename = paste(folder, paste0("boxplots_exprData_", nVarexpr, ".jpg"), sep = "/"))
      boxplot_nVarexpr
      dev.off()
      
      pdf(paste(folder, paste0("boxplots_exprData_", nVarexpr, ".pdf"), sep = "/"))
      print(boxplot_nVarexpr) # end of print
      dev.off()
    } # end of if, is.positive.integer()
  } # end of if in !is.null(nVarexpr)
  
  # expr data transformation
  if(!is.null(dataTransformation)){
    if (!is.null(nVarexpr)) d.trans <- getTransformedData(d = d[, 1:nVarexpr], dataTransformation = dataTransformation) else
      d.trans <- getTransformedData(d = d, dataTransformation = dataTransformation)
    
    yLabel <- switch (dataTransformation,
                      "log2" = "expression data (log2)",
                      "scale" = "expression data (scaled)",
                      "robustScale" = "expression data (median centered and IQR scaled)"
    ) # end of switch
    
    boxplot_trans <- ggplot(stack(data.frame(d.trans)), aes(x = ind, y = values)) + 
      geom_boxplot(color = box.outline.color, outlier.color = outlier.color, outlier.size = outlier.size) + 
      xlab(paste("#expr = ", ncol(d.trans))) + ylab(yLabel) +
      theme_classic() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) 
    
    jpeg(filename = paste(folder, paste0("boxplots_transformedExprData_", ncol(d.trans), ".jpg"), sep = "/"))
    boxplot_trans
    dev.off()
    
    pdf(paste(folder, paste0("boxplots_transformedExprData_", ncol(d.trans), ".pdf"), sep = "/"))
    print(boxplot_trans) # end of print
    dev.off()
  } # end of if in !is.null(dataTransformation)
  
  if (is.null(dataTransformation) & is.null(nVarexpr)) return(list("d" = d)) else
    if (!is.null(dataTransformation) & is.null(nVarexpr)) return(list("d" = d, "d.nVarexpr" = d.nVarexpr)) else
      if (is.null(dataTransformation) & !is.null(nVarexpr)) return(list("d" = d, "d.trans" = d.trans)) else
        return(list("d" = d, "d.nVarexpr" = d.nVarexpr, "d.trans" = d.trans))
  
  cat("Boxplots folder location: ", paste(folder))
} # end of function getExprBoxplots






