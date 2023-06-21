# association between outcome and predictor data, both categorical
# chi sq test of association is performed 
# bar plots with counts and proportions generated

get.association.barplot <- function(xydata, x, y, p){
  p1 <- ggplot(xydata, aes(x = xydata[, y], fill = xydata[, x])) +
    geom_bar(show.legend = FALSE) +
    xlab(y) + #ylab(pvar) +
    theme_bw(15)
  
  p2 <- ggplot(xydata, aes(x = xydata[, y], fill = xydata[, x])) +
    geom_bar(position = "fill") +
    xlab(y) + ylab("proportion") +
    labs(fill = x) +
    theme_bw(15)
  
  p1 + p2 + plot_annotation(title = paste("chisq p-value =", scientific(p))) 
}

getChiSq.association <- function(data, xNames, yNames, output_folder, ...){
  # later coding - check if the xNames and yNames are factors
  createDir(output_folder)
  for(ovar in yNames){
    # chi sq test of association between the 
    dfChi <- t(sapply(xNames, FUN = function(x){
      chi <- chisq.test(data[, ovar], data[, x])
      c(chi$statistic, chi$parameter, chi$p.value)
    }))
    colnames(dfChi) <- c("chisq.statistic", "df", "pvalue")
    
    assoc.barplot <- lapply(xNames, FUN = function(pvar){get.association.barplot(xydata = data, x = pvar, y = ovar, p = dfChi[pvar, "pvalue"])})

    write.csv(dfChi, file = paste0(output_folder, "/chisq_association_", ovar, ".csv"))
    pdf(file = paste0(output_folder, "/chisq_association_barplots_", ovar, ".pdf"))
    print(assoc.barplot)
    dev.off()
  } # end of for loop in ovar
} # end of function getChiSq.association
