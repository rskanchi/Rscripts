# association between outcome and predictor data, both categorical
# chi sq test of association is performed 
# bar plots with counts and proportions generated

get.association.barplot <- function(xydata, x, y, pval, theme_size =15){
  library(ggplot2)
  library(patchwork)
  library(scales)
  
  
  xydata <- na.omit(xydata)
  nY <- table(xydata[, y])
  
  p1 <- ggplot(xydata, aes(x = xydata[, y], fill = xydata[, x])) +
    geom_bar(show.legend = FALSE) +
    xlab(y) + #ylab(pvar) +
    theme_bw(theme_size)
  
  p2 <- ggplot(xydata, aes(x = xydata[, y], fill = xydata[, x])) +
    geom_bar(position = "fill") +
    xlab(y) + ylab("proportion") +
    labs(fill = x) +
    stat_n_text() +
    #geom_text(aes(label = paste("n =", nY))) +
    theme_bw(theme_size)
  
  p1 + p2 + plot_annotation(title = paste("chisq p-value =", scientific(pval))) 
}

getChiSq.association <- function(data, xNames, yNames, output_folder, ...){
  library(ggplot2)

  createDir(output_folder)
  
  # make the variables factor if not already
  for (v in c(xNames, yNames)){
    if (!is.factor(data[, v])) data[,v] <- factor(data[, v])
  } # end of for in v
  
  for(ovar in yNames){
    # chi sq test of association between the 
    dfChi <- t(sapply(xNames, FUN = function(x){
      chi <- chisq.test(data[, ovar], data[, x])
      c(chi$statistic, chi$parameter, chi$p.value)
    }))
    colnames(dfChi) <- c("chisq.statistic", "df", "pvalue")
    
    assoc.barplot <- lapply(xNames, FUN = function(pvar){get.association.barplot(xydata = data, x = pvar, y = ovar, pval = dfChi[pvar, "pvalue"])})

    write.csv(dfChi, file = paste0(output_folder, "/chisq_association_", ovar, ".csv"))
    pdf(file = paste0(output_folder, "/chisq_association_barplots_", ovar, ".pdf"))
    print(assoc.barplot)
    dev.off()
  } # end of for loop in ovar
} # end of function getChiSq.association
