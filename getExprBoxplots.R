# R libraries
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(ggraph)){install.packages("ggraph")}


getExprBoxplots <- function(d, dataTransformation = NULL, nVarexpr = NULL,
                            outlier.size = 0.2, outlier.color = "red", box.outline.color = "grey",
                            folder = NULL){
  library(ggplot2)
  if (is.null(folder)) folder <- paste(getwd(), "Output_consensus_clustering", sep = "/")
  pdf(paste(folder, paste("boxplots_exprData.pdf", sep = ""), sep = "/"))
  print(
    ggplot(stack(d), aes(x = ind, y = values)) + 
      geom_boxplot(color = box.outline.color, outlier.color = outlier.color, outlier.size = outlier.size) + 
      xlab(paste("#expr = ", ncol(d))) + ylab("expression data") +
      theme_classic() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) 
  )
  # expr data transformation
  if (dataTransformation == "log2") {
    d <- log2(d+1)
    print(ggplot(stack(d), aes(x = ind, y = values)) + 
            geom_boxplot(color = box.outline.color, outlier.color = outlier.color, outlier.size = outlier.size) + 
            xlab(paste("#expr = ", ncol(d))) + ylab("expression data (log2)") +
            theme_classic() +
            theme(axis.ticks.x = element_blank(),
                  axis.text.x = element_blank()) 
    )
    
    d <- sweep(d, 2, apply(d, 2, median))
    mads <- apply(d, 2, mad)
    d <- d[, names(mads)[rev(order(mads))]] # data sorted by mad
  } else if (dataTransformation == "scale") {# end of if in dataTransformation == "log2"
    d <- data.frame(scale(d))
  } else { # no log2 or scale, then median centering the expr data
    d <- sweep(d, 2, apply(d, 2, median))
    mads <- apply(d, 2, mad)
    d <- d[, names(mads)[rev(order(mads))]] # data sorted by mad
  } 
  
  # plot the expr data
  if (dataTransformation == "log2") yLab <- "expression data (log2, median centered)" else
    if (dataTransformation == "scale") yLab <- "expression data (scaled)" else
      yLab <- "expression data (median centered)"
  print(
    ggplot(stack(d), aes(x = ind, y = values)) + 
      geom_boxplot(color = box.outline.color, outlier.color = outlier.color, outlier.size = outlier.size) + 
      xlab(paste("#expr = ", ncol(d))) + ylab(yLab) +
      theme_classic() +
      theme(axis.ticks.x = element_blank(),
            axis.text.x = element_blank()) 
  )
  # get highly variable genes
  if (!is.null(nVarexpr)) {
    
    lapply(nVarexpr, FUN =function(v){
      tempD <- d[, names(mads)[rev(order(mads))][1:v]]
      print(
        ggplot(stack(tempD), aes(x = ind, y = values)) + 
          geom_boxplot(color = box.outline.color, outlier.color = outlier.color, outlier.size = outlier.size) + 
          xlab(paste("#highly variable expr = ", ncol(tempD))) + ylab(yLab) +
          theme_classic() +
          theme(axis.ticks.x = element_blank(),
                axis.text.x = element_blank()) 
      )
    }) # end of lapply in nVarexpr
  } # end of !is.null(nVarexpr)
  dev.off()
  return(d)
  cat("Boxplots saved in file: ", paste(folder, paste("boxplots_exprData.pdf", sep = ""), sep = "/"))
} # end of function getExprBoxplots
