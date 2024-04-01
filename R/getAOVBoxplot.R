# function to get box and violin plots for the specified continuous variables in a data set
# data with sample rows and columns including xNames and yNames
# specify the names of continuous traits in the argument yNames: yNames = colnames()[c()] or yNames = c("","","")

getAOVBoxplot <- function(data, xNames, yNames, folder = NULL,
                          theme_size = 20, text_size = 10,
                          pdf.height = 8, pdf.width = 8){
  library(ggplot2)
  library(gridExtra)
  library(scales) # scientific()
  if (is.null(folder)) folder <- "output/AOV_Boxplots"
  createDir(folder = folder)
  
  for (x in xNames){
    dfAOV <- lapply(yNames, FUN = function(y){
      aov(data[, y] ~ data[, x] , data = data)
    }) # end of sapply in y
    names(dfAOV) <- yNames
    
    box_plot <- lapply(yNames, FUN = function(y){
      aov.p <- summary(dfAOV[[y]])[[1]][1,5]
      tbl.n <- table(data[, x])
      ggplot(data, aes(x = factor(data[, x]), y=data[,y])) +
        geom_violin(trim = TRUE, width = 0.3, size = 1.5) +
        geom_boxplot(width=0.1, color="grey", alpha=0.1) +
        stat_summary(fun = "mean",
                     geom = "point",
                     color = "red") +
        xlab(x) + ylab(y) +
        labs(caption = paste(paste("AOV p-value: ", scientific(aov.p), sep = ""),
              paste(names(tbl.n), tbl.n, sep = " = ", collapse = ", "), sep = "\n")) +
        #labs(caption = paste("AOV p-value: ", scientific(aov.p), sep = "")) +
        theme_bw(theme_size)
    }) # end of lapply in y
    
    pdf(paste(folder, paste0("AOV_Boxplots_", x, ".pdf"), sep = "/"), height = pdf.height, width = pdf.width)
    print(box_plot)
    dev.off()
  } # end of for loop in x
  #return(box_plot)
} # end of function getAOVBoxplot
