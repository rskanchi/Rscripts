# function to get histograms for the specified continuous columns of a data set
# data with sample rows and columns as traits
# specify the names of traits in the argument traits: traits = colnames()[c()] or traits = c("","","")

getHistogram <- function(data, traits, folder = NULL,
                   theme_size = 15, text_size = 7){
  library(ggplot2)
  library(gridExtra)

  if (is.null(folder)) folder <- "output/histograms"
  createDir(folder = folder)
  
  hist_plot <- list()
  vio_plot <- list()
  pdf(paste(folder, "histograms.pdf", sep = "/"), height = 8, width = 8)
  
  for (t in traits){
    vio_plot[[t]] <- ggplot(data, aes(x="", y=data[,t])) +
      geom_violin(trim = TRUE) +
      geom_boxplot(width=0.1, color="blue", alpha=0.2) +
      stat_summary(fun = "mean",
                   geom = "point",
                   color = "red") +
      xlab(t) + ylab("") +
      coord_flip() +
      theme_bw(theme_size)
    
    hist_plot[[t]] <- ggplot(data, aes(x = data[,t])) +
      geom_histogram(colour = "black", fill = "grey", binwidth = 5) + #bins = 15, 
      xlab(t) + ylab("Frequency") +
      theme_bw(theme_size)
    
    grid.arrange(vio_plot[[t]], hist_plot[[t]], ncol = 1, 
                 layout_matrix = rbind(c(1,1), c(2,2), c(2,2)))
    
  } # end of for loop in t
  
  #print(histograms)
  dev.off()
  
  return(list("violins" = vio_plot,  "histograms" = hist_plot))
} # end of function getHorizontalBarplot