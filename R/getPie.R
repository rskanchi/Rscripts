# function to get a pie chart for the specified discrete columns of a data set
# data with sample rows and columns as traits
# specify the names of traits in the argument traits: traits = colnames()[c()] or traits = c("","","")

getPie <- function(data, traits, folder = NULL,
                   theme_size = 15, text_size = 7){
  library(dplyr)
  library(ggplot2)
  library(scales)
  if (is.null(folder)) folder <- "output/piecharts"
  createDir(folder = folder)
  
  piecharts <- list()
  # get the frequencies and % for each trait
  for (t in traits){
    tempdf <- data %>% 
      group_by(data[,t]) %>% # Variable to be transformed
      count() %>% 
      ungroup() %>% 
      mutate(perc = `n` / sum(`n`)) %>% 
      arrange(perc) %>%
      mutate(labels = scales::percent(perc))
    colnames(tempdf)[1] <- "var"
    
    piecharts[[t]] <- ggplot(tempdf, aes(x = "", y = perc, fill = var)) +
      geom_col(color = "black") +
      geom_text(aes(label = labels),
                 position = position_stack(vjust = 0.5),
                 size = text_size,
                 show.legend = FALSE) +
      guides(fill = guide_legend(title = t)) +
      coord_polar(theta = "y") +
      #theme(text = element_text(size = 15)) +
      theme_void(theme_size)

    colnames(tempdf)[1] <- t
    write.csv(tempdf, file = paste(folder, paste0("Supporting_data_Piechart_", t, ".csv"), sep = "/"))
  } # end of for loop in t
  
  pdf(paste(folder, "Piecharts.pdf", sep = "/"), height = 8, width = 8)
  print(piecharts)
  dev.off()
  
  return(piecharts)
} # end of function getPie