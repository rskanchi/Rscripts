# get Jaccard similarity and distance for two sets
# Ex: set1 and set2 can be numeric or char sets 

getJaccard <- function(set1, set2){
  intersection <- length(intersect(set1, set2))
  union <- length(set1) + length(set2) - intersection
  jSimilarity <- intersection/union
  jDistance <- 1 - jSimilarity
  return(list("nIntersection" = intersection, "nUnion" = union, "similarity" = jSimilarity, "distance" = jDistance))
} # end of function getJaccard
