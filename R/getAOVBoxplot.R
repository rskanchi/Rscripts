# function to get box and violin plots for the specified continuous variables in a data set
# data with sample rows and columns including xNames and yNames
# specify the names of continuous traits in the argument yNames: yNames = colnames()[c()] or yNames = c("","","")

getAOVBoxplot <- function(data, xNames, yNames, folder = NULL,
                         theme_size = 20, text_size = 10){
  library(ggplot2)
  library(gridExtra)
  
  if (is.null(folder)) folder <- "output/AOV_Boxplots"
  createDir(folder = folder)
  
  box_plot <- list()
  for (x in xNames){
    for (y in yNames){
      bname <- paste0("xy_", x, y) # boxplot name, combination of x and y
      
      aov <- aov(data[, y] ~ data[, x] , data = data)
      aov.p <- summary(aov)[[1]][1,5]
      
      box_plot[[bname]] <- ggplot(data, aes(x = factor(data[, x]), y=data[,y])) +
        geom_violin(trim = TRUE, width = 0.3, size = 1.5) +
        geom_boxplot(width=0.1, color="grey", alpha=0.1) +
        stat_summary(fun = "mean",
                     geom = "point",
                     color = "red") +
        xlab(x) + ylab(y) +
        labs(caption = paste("AOV p-value: ", round(aov.p, 4), sep = "")) +
        theme_bw(theme_size)
    } # end of for loop in t
  } # end of for loop in x
  
  pdf(paste(folder, "AOV_Boxplots.pdf", sep = "/"), height = 8, width = 8)
  print(box_plot)
  dev.off()
  
  return(list("boxplots" = box_plot))
} # end of function getAOVBoxplot