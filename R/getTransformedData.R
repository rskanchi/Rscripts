# function to get transformed data
# dataTransformation = "log2", "scale", "robustScale" # default "log2"
  # robustScale: center using median and scale using either IQR (inter-quartile range) or MAD (median absolute dev), default = IQR

getTransformedData <- function(d, dataTransformation = "log2", 
                               scl.func = "IQR") { # scl = IQR or MAD, used only for dataTransformation = "robustScale"
  res <- switch (dataTransformation,
                 "log2" = log2(d+1),
                 "scale" = scale(d),
                 "robustScale" = {
                   ctr <- apply(d, MARGIN = 2, median, na.rm = TRUE) # for median centering
                   # for scaling using median absolute deviation, mad or "IQR"
                   if (scl.func == "IQR"){scl <- apply(d, MARGIN = 2, IQR, na.rm = TRUE)} else
                     if (scl.func == "MAD"){scl <- apply(d, MARGIN = 2, mad, na.rm = TRUE)}
                   sweep(sweep(d, MARGIN = 2, ctr, "-"), 
                         MARGIN = 2, scl, "/")
                   
                 } # end of case median_IQR
  ) # end of switch
  return(res)
} # end of function getTransformedData
