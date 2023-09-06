computeDiffExpr_limma <- function(d, group_var, group_names, # d = expr data + the group variable
                            thresholdAdjPVAL = 0.05, padj_method = "BH", topn = 10, plot = FALSE,
                            folder = NULL, file_id = NULL,   text_size = 3, rowsAsSamples = TRUE){
  library(limma)
  library(ggplot2)
  library(ggrepel)
  
  indices <- which(group_var %in% group_names)
  d <- d[indices,]
  group_var <- group_var[indices]
  group_var <- factor(group_var, labels = group_names, ordered = TRUE)

  # with covariate adjustment
  #design <- model.matrix(~ 0 + Group + data[, xcovariates[1]] + data[, xcovariates[2]] + data[, xcovariates[3]])
  design <- model.matrix(~ 0 + group_var)
  colnames(design) <- group_names
  contrast_name <- paste(group_names[2], "over", group_names[1], sep = "")
  contrast_expr <- paste(group_names[2], "-", group_names[1], sep = "")
  # automate makeContrasts
  cmd <- paste("cont.matrix <- makeContrasts(", contrast_name, "=", contrast_expr, ", levels = design)", sep = '"')
  eval(parse(text = cmd))
  
  # volcano plot and output table
  deColors <- c("blue", "red", "black")
  names(deColors) <- c("DOWN", "UP", "NO")
  #linetype <- "dashed"
  
  res <- matrix(ncol = 8)
  colnames(res) <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "diffexpressed", "delabel")
  if (rowsAsSamples) fit <- lmFit(t(d), design = design) else fit <- lmFit(d, design = design)
  
  fit.cont <- contrasts.fit(fit, cont.matrix)
  fit.cont <- eBayes(fit = fit.cont)
  
  res <- topTable(fit = fit.cont, number = ncol(d), coef = 1, adjust.method = padj_method)
  # column for differentially expressed, UP or DOWN
  res$diffexpressed <- "NO"
  res$diffexpressed[res$logFC > 0 & res$adj.P.Val < thresholdAdjPVAL] <- "UP"
  res$diffexpressed[res$logFC < 0 & res$adj.P.Val < thresholdAdjPVAL] <- "DOWN"
  
  diffexpr <- rownames(res)[which(res$diffexpressed != "NO")]
  
  # save files
  if (!is.null(folder)){
    if (!dir.exists(paste(folder))){dir.create(paste(folder), recursive = TRUE)}
    if (is.null(file_id)) file_id <- ""
    fileName <- paste(folder, paste(contrast_name, file_id, "padj", padj_method, thresholdAdjPVAL, sep = "_"), sep = "/")
    openxlsx::write.xlsx(res, file = paste(fileName, "DiffExp_summary.xlsx", sep = "_"), rowNames = TRUE, overwrite = TRUE)
  } # end of saving to file
  
  
  if (plot){
    res$delabel <- NA
    #res$delabel[res$diffexpressed != "NO"] <- rownames(res)[res$diffexpressed != "NO"]
    # many significantly UP and DOWN expr to label. So, just labelling the top n on each side
    top_up <- min(topn, sum(res$diffexpressed == "UP"))
    top_down <- min(topn, sum(res$diffexpressed == "DOWN"))
    
    res$delabel[which(res$diffexpressed == "UP")[1:top_up]] <- rownames(res)[which(res$diffexpressed == "UP")[1:top_up]]
    res$delabel[which(res$diffexpressed == "DOWN")[1:top_down]] <- rownames(res)[which(res$diffexpressed == "DOWN")[1:top_down]]
    
    xlim <- c(-floor(max(abs(res$logFC), na.rm = TRUE)+1), floor(max(abs(res$logFC), na.rm = TRUE)+1))
    
    tblGroup <- table(group_var)
    volcano_plot <- ggplot(data = res, aes(x = logFC, y = -log10(P.Value), col = diffexpressed, label = delabel)) +
      xlim(xlim = xlim) +
      geom_point(size = 1) +
      scale_color_manual(values = deColors) +
      geom_text_repel(size = text_size, max.overlaps = Inf) +
      #geom_hline(yintercept = -log10(thresholdAdjPVAL), lty = linetype) +
      xlab("log2 fold change") +
      labs(title = contrast_name,
           subtitle = paste("Total =", nrow(res),
                            "| DOWN =", sum(res$diffexpressed == "DOWN", na.rm = TRUE),
                            "| UP =", sum(res$diffexpressed == "UP", na.rm = TRUE),
                            sep = " "),
           caption = paste("sample sizes | ", paste(names(tblGroup), tblGroup, sep = " = ", collapse = ", "), "\n",
                           "p-value threshold: ", thresholdAdjPVAL, "\n",
                           "p-value adjustment method: ", padj_method, "\n",
                           "Top ", top_up, " UP and ", top_down, " DOWN genes labelled", sep = "")) +
      theme_bw(15)
    
    if (!is.null(folder)){
      if (is.null(file_id)) file_id <- ""
      pdf(file = paste(fileName, "volcanoplot.pdf", sep = "_"))
      print(volcano_plot)
      dev.off()
    }
    
    
    return(list(volcano = volcano_plot, result = res, diffexpr = diffexpr))
    
  } else return(list(result = res, diffexpr = diffexpr))
  
  
} # end of computeDiffExpr
