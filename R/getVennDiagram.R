# function to get the proportional venn diagrams 
# data as a csv file with set variables 

# if signVariables are specified, venn diagrams are generated for all, UP and DOWN lists
# order of signVariables, if given, must be the same as setVariables
# order of setLabels, if given, must be the same as setVariables

getVennDiagram <- function(data, setVariables, signVariables = NULL, setLabels = NULL, setColors = NULL,
                           folder,
                           alpha = 0.7, lwd = 2, label.size = 15
                           , ...){
  library(eulerr)
  library(gplots)
  library(stringi)

  # all elelments in the sets used
  setList <- as.list(data[, setVariables])
  setList <- lapply(setList, FUN = function(x) x[!is.na(x)])
  if (!is.null(setLabels)) names(setList) <- setLabels
  saveVenn(setList = setList, folder = folder, label.size = label.size, lwd = lwd, setColors = setColors)
  
  # if the sign variables are given to get the UP and DOWN overlaps and venn diagrams
  if (!is.null(signVariables)){
    signList <- as.list(data[, signVariables])
    signList <- lapply(signList, FUN = function(x) x[!is.na(x)])

    # UP 
    setList_UP <- lapply(1:length(setList), FUN = function(x) {setList[[x]][signList[[x]] > 0]})
    names(setList_UP) <- names(setList)
    saveVenn(setList = setList_UP, folder = folder, vennLabel = "UP", label.size = label.size, lwd = lwd, setColors = setColors)
    
    # DOWN
    setList_DOWN <- lapply(1:length(setList), FUN = function(x) {setList[[x]][signList[[x]] < 0]})
    names(setList_DOWN) <- names(setList)
    saveVenn(setList = setList_DOWN, folder = folder, vennLabel = "DOWN", label.size = label.size, lwd = lwd, setColors = setColors)
  } # end of !is.null(signVariables)
} # end of function getVennDiagram


saveVenn <- function(setList, folder, vennLabel = "ALL", alpha = 0.7, lwd = 2, label.size = 15, setColors = NULL, ...){
  
  createDir(folder)
  nSets <- length(setList)
  if (is.null(setColors)){
    # colors: blue, lightorange, pink, orange, green, lightblue, yellow, black, grey
    clrs  <- c("#0072B2", "#E69F00", "#CC79A7", "#D55E00", "#009E73", "#56B4E9", "#F0E442", "#000000", "#999999")
    # pie(rep(1, 9), col = clrs)
    setColors <- clrs[1:nSets]
  }
  
  eulerData <- euler(setList)
  write.csv(x = data.frame(eulerData$original.values), 
            file = paste(folder, paste0("euler_data_", vennLabel, ".csv"), sep = "/"))
  pdf(file = paste(folder, paste0("venn_", vennLabel, ".pdf"), sep = "/")) #, width = 15, height = 15
  #par(mar = c(1,1,1,1))
  print(
      plot(eulerData, quantities = list(fontsize = label.size), 
           edges = list(col = setColors, alpha = alpha, lwd = lwd),
           fills = list(fill = setColors, alpha = alpha),
           legend = list(fontsize = label.size, cex = 1.2), labels = list(fontsize = label.size), 
           main = list(label = vennLabel, fontsize = label.size, cex = 1.2),
           adjust_labels = TRUE)
    )
  dev.off()
  
  vennData <- gplots::venn(setList, intersections = TRUE, showSetLogicLabel = TRUE, 
                           show.plot = FALSE, simplify = TRUE)
  intersectionList <- as.list(attr(vennData, "intersection"))
  #names(intersectionList)
  intersectionDF <- as.data.frame(stri_list2matrix(intersectionList), row.names = F)
  colnames(intersectionDF) <- names(intersectionList)
  write.table(intersectionDF, 
              file = paste(folder, paste0("venn_intersection_data_", vennLabel, ".csv"), sep = "/"), 
              quote = F, sep = ",", row.names = F)
} # end of saveVenn
