# function to get a bar chart for the specified discrete columns of a data set
# data with sample rows and columns as traits
# specify the names of traits in the argument traits: traits = colnames()[c()] or traits = c("","","")

getHorizontalBarplot <- function(data, traits, folder = NULL,
                   theme_size = 15, text_size = 7){
  library(dplyr)
  library(ggplot2)
  library(scales)
  
  if (is.null(folder)) folder <- "output/barcharts"
  createDir(folder = folder)
  
  barcharts <- list()
  # get the frequencies and % for each trait and then plot
  for (t in traits){
    tempdf <- data %>% 
      group_by(data[,t]) %>% # Variable to be transformed
      count() %>% 
      ungroup() %>% 
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    colnames(tempdf)[1] <- "var"
    tempdf$var <- factor(tempdf$var, levels = tempdf$var, ordered = TRUE)
    barcharts[[t]] <- ggplot(data=tempdf, aes(x=var, y=perc*100)) +
      geom_segment( aes(x=var, xend=var, y=0, yend=perc*100), color = "grey", size = 1.5) + 
      geom_point(color = "black", size = 5) +
      xlab(t) + ylab("Size (%)") +
      coord_flip() +
      theme(
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      theme_bw(theme_size)

    colnames(tempdf)[1] <- t
    write.csv(tempdf, file = paste(folder, paste0("Supporting_data_Barcharts_", t, ".csv"), sep = "/"))
  } # end of for loop in t
  
  pdf(paste(folder, "barcharts.pdf", sep = "/"), height = 8, width = 8)
  print(barcharts)
  dev.off()
  
  return(barcharts)
} # end of function getHorizontalBarplot