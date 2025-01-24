# ORA
performORA <- function(sig.d, # input gene list, log2FC or cor to identify up and down genes
                       pvaleCutoff = 0.05, # for identifying over represented pathways
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.25,
                       compendium, # name of the compendium ex "GOBP"
                       genesets.gmt, # gmt file with path from msigdb for the compendium of choice
                       minGSSize = NULL,
                       maxGSSize = NULL,
                       universe = NULL, # background genes: All genes identified in the study could be specified. If NULL, all genes in the database (TERM2GENE table) will be used.
                       output_folder = "output",
                       top_n = 10,
                       theme_size = 10,
                       ...
                       ){
  
  # functions
  library(clusterProfiler)
  library(dplyr)
  
  genesets <- read.gmt(gmtfile = genesets.gmt)
  genesets$term <- as.character(genesets$term)
  genesets$term <- gsub(" ", ".", genesets$term)
  if (compendium == "Epigenomics"){
    temp <- strsplit(genesets$gene, ",")
    genesets$gene <- sapply(temp, FUN = function(x) x[1])
  } # end of if in epigenomics
  
  genesets <- unique(genesets) # removing duplicated rows; can happen with user defined genesets
  pathways_n <- length(unique(genesets$term))
  
  # making sure the genesets are included  
  if(!is.null(universe)){
    uni_diff <- setdiff(universe, unique(genesets$gene))
    if (length(uni_diff) > 0){
      genesets2 <- data.frame(term = rep("manual", length(uni_diff)), gene = uni_diff)
      genesets <- rbind(genesets, genesets2)
      rm(genesets2)
    }}
  
  tbl_gs_name <- table(genesets$term)
  if(is.null(minGSSize)) minGSSize <- min(tbl_gs_name)
  if(is.null(maxGSSize)) maxGSSize <- max(tbl_gs_name)
  
  signature_name <- colnames(sig.d)[1]
  
  sig_ALL <- toupper(sig.d[[1]]) 
  sig_UP <- toupper(sig.d[[1]])[sig.d[[2]] > 0]
  sig_DOWN <- toupper(sig.d[[1]])[sig.d[[2]] < 0]
  
#  df_term2gene <- df_msigdbr_compendium[, c("gs_name", gene_key)]
#  head(df_term2gene)  
  
  if (length(sig_ALL) > 0){
    ora_ALL <- enricher(gene = sig_ALL, TERM2GENE = genesets, #df_term2gene, 
                        pvalueCutoff = pvaleCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
                        minGSSize = minGSSize, maxGSSize = maxGSSize
    )
  }

  if (length(sig_UP) > 0){
    ora_UP <- enricher(gene = sig_UP, TERM2GENE = genesets,# df_term2gene, 
                       pvalueCutoff = pvaleCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
                       minGSSize = minGSSize, maxGSSize = maxGSSize
    )
  }
    

  if (length(sig_DOWN) > 0){
    ora_DOWN <- enricher(gene = sig_DOWN, TERM2GENE = genesets, #df_term2gene, 
                         pvalueCutoff = pvaleCutoff, pAdjustMethod = pAdjustMethod, qvalueCutoff = qvalueCutoff,
                         minGSSize = minGSSize, maxGSSize = maxGSSize
    )
  }
  
  
  # results write to file
  res_ora <- list("ora_ALL" = ora_ALL@result, "ora_UP" = ora_UP@result, "ora_DOWN" = ora_DOWN@result)
  # remove manual
  res_ora <- lapply(res_ora, FUN = function(x){
    if (any(x$ID == "manual")){
      manual_index <- grep("manual", x$ID)
      x <- x[-manual_index, ]
    } else x <- x
  })
  
  tbl_gs_name <- data.frame(tbl_gs_name)
  colnames(tbl_gs_name) <- c("geneset", "nGenes")
  if (any(tbl_gs_name$geneset == "manual")){
    manual_index <- grep("manual", tbl_gs_name$geneset)
    tbl_gs_name <- tbl_gs_name[-manual_index, ]
  }

  
    
  output_folder <- paste(output_folder, signature_name, compendium, sep = "/")
  
  if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  write.csv(res_ora$ora_ALL, file = paste(output_folder, paste0(signature_name, "_ORA_ALL_", compendium, ".csv"), sep = "/"), row.names = F)
  write.csv(res_ora$ora_UP, file = paste(output_folder, paste0(signature_name, "_ORA_UP_", compendium, ".csv"), sep = "/"), row.names = F)
  write.csv(res_ora$ora_DOWN, file = paste(output_folder, paste0(signature_name, "_ORA_DOWN_", compendium, ".csv"), sep = "/"), row.names = F)
  write.csv(tbl_gs_name, file = paste(output_folder, paste0(compendium, "_geneset_size.csv"), sep = "/"), row.names = F)
  
  # visualization
  library(ggplot2)
  # filter genesets based on padj
  sig.res_ORA <- lapply(res_ora, FUN = function(x){
    x <- x %>% 
      select(pvalue, p.adjust, Count) %>%
      filter(p.adjust < pvaleCutoff)
    
    if (nrow(x) > 0){
      negLog10pvalue <- -log10(x[,"pvalue"])
      x <- cbind(x, negLog10pvalue)
    }
  })
  
  # remove null list elements (i.e. if ALL, UP or DOWN do not have sig OR genesets)
  sig.res_ORA <- sig.res_ORA[!sapply(sig.res_ORA, is.null)]
  
  
  lapply(sig.res_ORA, dim)
  lapply(sig.res_ORA, FUN = function(x) x[1:min(top_n, nrow(x)),])
  # create plots if not null
  if (length(sig.res_ORA) > 0){
  
  user_theme <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    
    # lollipop plots
    dir.create(paste(output_folder, "plotdata", sep = "/"))
    lolli_plots <- lapply(1:length(sig.res_ORA), FUN = function(i){
      x <- sig.res_ORA[[i]]
      x$gs <- rownames(x)
      sig_n <- nrow(x)
      top_n <- min(sig_n, top_n)
      x <- x[order(x$pvalue),]
      x <- x[1:top_n,]
      
      write.csv(x, paste(output_folder, "plotdata", paste0(names(sig.res_ORA)[i], ".csv"), sep = "/"))
      ggplot(data = x, aes(x=reorder(gs, negLog10pvalue), y=negLog10pvalue)) +
        geom_segment(aes(xend=gs, yend=0), linewidth = 0.8) +
        geom_point(mapping=aes(size=Count), color="black") +
        scale_size_area(name = "overlap") + #breaks = breaks, labels = c("min", "median", "max")
        coord_flip() +
        theme_bw(theme_size) + user_theme +
        xlab("") + ylab("-log10(pvalue)") + 
        labs(title = paste(compendium, signature_name, names(sig.res_ORA)[i], sep = " | "),
             subtitle = paste0("pathways = ", pathways_n, " | over-represented pathways = ", sig_n, "\n padj method = ", pAdjustMethod, " | padj threshold = ", pvaleCutoff))
    }) 
    
    names(lolli_plots) <- names(sig.res_ORA)
    
    pdf(paste(output_folder, paste0("lollipop_plots_top_pathways", signature_name, "_", compendium, ".pdf"), sep = "/"), height = 4, width = 12)
    print(lolli_plots)
    dev.off()


    # bar plots
    bar_plots <- lapply(1:length(sig.res_ORA), FUN = function(i){
      x <- sig.res_ORA[[i]]
      x$gs <- rownames(x)
      sig_n <- nrow(x)
      top_n <- min(sig_n, top_n)
      x <- x[order(x$pvalue),]
      x <- x[1:top_n,]
      

      ggplot(data = x, aes(x = reorder(gs, negLog10pvalue) , y = negLog10pvalue)) +
        geom_col(width = 0.5, aes(fill = Count)) +
        labs(x = "", y = "-log10(pval)") +
        scale_fill_gradient2(low = "grey", high = "black") +
        coord_flip() +
        theme_bw(theme_size) + user_theme + 
        xlab("") + ylab("-log10(pvalue)") +
        labs(title = paste(compendium, signature_name, names(sig.res_ORA)[i], sep = " | "),
             subtitle = paste0("pathways = ", pathways_n, " | over-represented pathways = ", sig_n, "\n padj method = ", pAdjustMethod, " | padj threshold = ", pvaleCutoff))
      
    }) 
    
    names(bar_plots) <- names(sig.res_ORA)
    
    pdf(paste(output_folder, paste0("bar_plots_top_pathways", signature_name, "_", compendium, ".pdf"), sep = "/"), height = 4, width = 12)
    print(bar_plots)
    dev.off()
  } # end of create plots if not null
 
  return(res_ora)
  
} # end of function performORA

