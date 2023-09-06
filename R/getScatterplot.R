# to get scatter plot with correlation and R2 values, a best fit line

# add_fit = "none", "reg.line" for linear regression line, "loess" for local regression fitting

getScatterplot <- function(data, xName, yNames, byGroup = NULL, 
                           cor.method = "pearson", theme_size = 15, scatter_point_size = 1,
                           add_fit = "reg.line", # "none", "reg.line" "loess"
                           color_fit = "red",
                           group.by = NULL,
                           plot.title = NULL,
                           
                           
                           ...){
  library(ggpubr)
  
  ggscatter(data = data, x = data[, xName], y = data[, yNames], combine = TRUE,
            add.params = list(color = color_fit),
            #cor.coef = TRUE,
            size = scatter_point_size, title = plot.title, 
            add = add_fit, rug = TRUE) +
    #stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~`,`~"))) +
    stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., 
                                             ..rr.label.., sep = "~`,`~"))) +
    #stat_regline_equation(aes(label =  paste(..eq.label..))) +
    stat_cor(method = cor.method, digits = 3, label.x.npc = "center",
             aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) + # default top left   #, label.x = 3, label.y = 70
    theme_bw(theme_size)
  
  if (!is.null(facet.by)){
    ggscatter(data = data, x = data[, xName], y = data[, yNames], combine = TRUE,
              add.params = list(color = color_fit),
              #cor.coef = TRUE,
              size = scatter_point_size, title = plot.title, 
              facet.by = group.by,
              add = add_fit, rug = TRUE) +
      #stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~`,`~"))) +
      stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., 
                                               ..rr.label.., sep = "~`,`~"))) +
      #stat_regline_equation(aes(label =  paste(..eq.label..))) +
      stat_cor(method = cor.method, digits = 3, label.x.npc = "center",
               aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~"))) + # default top left   #, label.x = 3, label.y = 70
      theme_bw(theme_size)
  }
  
  
} # end of function getScatterplot
  


