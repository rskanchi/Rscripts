# Function to perform time to event analysis such as survival analysis
# Two ways these data are different from linear regression outcome variables:
    # 1. Event times typically skewed right; violate normality assumption
    # 2. Censored data: may have unknown starting event (left censoring) and/or ending event (right censoring)

# Survival analysis models the time until an event occurs - compare time-to-event between different groups, or how time-to-event correlates with quantitative variables

# Models
# Parametric: treats time to event as continuous outcome with survival times following Weibull dist
# Discrete time : treats time to event as a series of person-periods with binary outcomes that follow a logit or cumulative log-log distribution
# Semi-parametric: Unspecified time to event dist estimated by non-parameteric methods coupled with covariate effects 
                  # following a parametric distribution. The Cox proportional hazards model is a semi-parametric model
                  # typically Kaplan-Meier plots to visualize survival curves and non-parametric log-rank tests to compare the curves among groups
                  # Cox proportional hazards model to to describe the effect of explanatory variables on survival
# Machine learning: trees, SVMs to optimize predictive power



# Resource: https://bookdown.org/mpfoley1973/survival/
# https://thriv.github.io/biodatasci2018/r-survival.html


# R: survival package to model, survminer to visualize, and gtsummary for summarize
# survfit: three tests for the overall significance of the model likelihood ratio, wald, score.
# Methods are asymptotically equivalent. The Likelihood ratio test has better behavior for small sample sizes, so it is generally preferred.

survContinuous <- function(expr, # predictors, to perform single predictor analysis
                           meta, time, status, covariates = NULL, # meta data including time, status and if available, the covariates
                           threshold_pval = 0.05, topn = 10, plot = FALSE,
                           text_size = 3, rowsAsSamples = TRUE, output_folder = NULL){
  
  library(tidyverse)
  library(scales)
  library(survival)
  library(survminer)
  library(gtsummary)
  library(ggrepel)
  library(forcats)
  
  res <- matrix(nrow = ncol(expr) , ncol = 15)
  
  for (i in 1:ncol(expr)){
    pred <- colnames(expr)[i]
    temp <- merge(expr[, pred, drop = F], meta[, c(time, status, covariates)], by = 0)
    if (is.null(covariates)) formula <- paste0("Surv(", time, ", ", status, ") ~ ", pred) else
      formula <- paste0("Surv(", time, ", ", status, ") ~ ", paste(pred, covariates, sep = " + "))
    cox <- coxph(as.formula(formula), data = temp)
    summ <- summary(cox)
    names(summary(cox))
    res[i,] <- c(summ[["logtest"]], summ[["coefficients"]][1,-2], summ[["conf.int"]][1,], summ[["n"]], summ[["nevent"]], summ[["concordance"]])
  } # end of for loop in pred
  
  colnames(res) <- c("LR_test", "LR_df", "LR_pval", "coeff", "se_coef", "z", "z_pval", "exp_coef", "exp_neg_coef", 
                     "CI_lower_95", "CI_upper_95", "n", "n_event", "C", "se_C")
  rownames(res) <- colnames(expr)
  res <- data.frame(res)
  res <- res[order(res$z_pval), ]
  #p.adjust(res$z_pval, method = "BH")
  # volcano plot
  deColors <- c("blue", "red", "black")
  names(deColors) <- c("DOWN", "UP", "NO")
  res$diffexpressed <- "NO"
  res$diffexpressed[res$coeff > 0 & res$z_pval < threshold_pval] <- "UP"
  res$diffexpressed[res$coeff < 0 & res$z_pval < threshold_pval] <- "DOWN"
  
  diffexpr <- rownames(res)[which(res$diffexpressed != "NO")]
  
  # save files
  if (is.null(output_folder)) output_folder <- paste("output/survival", time, sep = "/") else
    output_folder <- paste(output_folder, time, sep = "/")
  createDir(output_folder)
  
  fileName <- paste(output_folder, paste0("survivalAnalysis_summary_", time, "_pvalThreshold_", threshold_pval, ".csv"), sep = "/")
  write.csv(res, file = fileName, row.names = T)
  
  if (plot){
    res$delabel <- NA
    # many significantly UP and DOWN expr to label. So, just labelling the top n on each side
    top_up <- min(topn, sum(res$diffexpressed == "UP"))
    top_down <- min(topn, sum(res$diffexpressed == "DOWN"))
    
    res$delabel[which(res$diffexpressed == "UP")[1:top_up]] <- rownames(res)[which(res$diffexpressed == "UP")[1:top_up]]
    res$delabel[which(res$diffexpressed == "DOWN")[1:top_down]] <- rownames(res)[which(res$diffexpressed == "DOWN")[1:top_down]]
    
    xlim <- c(-floor(max(abs(res$coeff), na.rm = TRUE)+1), floor(max(abs(res$coeff), na.rm = TRUE)+1))
    
    tblGroup <- table(meta[, status])
    volcano_plot <- ggplot(data = res, aes(x = coeff, y = -log10(z_pval), col = diffexpressed, label = delabel)) +
      xlim(xlim = xlim) +
      geom_point(size = 1) +
      scale_color_manual(values = deColors) +
      geom_text_repel(size = text_size, max.overlaps = Inf) +
      #geom_hline(yintercept = -log10(thresholdAdjPVAL), lty = linetype) +
      xlab("coefficient estimate, cox PH model") + ylab("-log10(p-value)") +
      labs(title = paste0("time to event = ", time, "\n", "status = ", status),
           subtitle = paste("Total =", nrow(res),
                            "| DOWN =", sum(res$diffexpressed == "DOWN", na.rm = TRUE),
                            "| UP =", sum(res$diffexpressed == "UP", na.rm = TRUE)),
           caption = paste(status, " sample sizes | ", paste(names(tblGroup), tblGroup, sep = " = ", collapse = ", "), "\n",
                           "p-value threshold: ", threshold_pval, "\n",
                           #"p-value adjustment method: ", padj_method, "\n",
                           "Top ", top_up, " UP and ", top_down, " DOWN labelled", sep = "")) +
      theme_bw(15)
    
    # exp(coef) plot
    if (length(diffexpr) > 0) CIdata <- res[diffexpr, , drop = F] else CIdata <- res[1:topn,]
    
    CIdata$expr <- rownames(CIdata)
    CIdata <- CIdata[, c("exp_coef", "CI_lower_95", "CI_upper_95", "expr")]
    CIdata$expr <- fct_rev(factor(CIdata$expr, levels = CIdata$expr, ordered = T)) 
    
    CIplot <- ggplot(CIdata) +
      geom_segment(aes(x=expr, xend=expr, y=CI_lower_95, yend=CI_upper_95), color="#000000") +
      geom_point( aes(x=expr, y=exp_coef), color="#000000", size=3, shape = 8) +
      geom_point( aes(x=expr, y=CI_lower_95), color="#0072B2", size=3) +
      geom_point( aes(x=expr, y=CI_upper_95), color="#D55E00", size=3) +
      geom_hline(yintercept = 1, lty = "dashed", col = "black") +
      coord_flip()+
      theme_bw(15) +
      theme(legend.position = "none") +
      xlab("") + ylab("95% CI for exp(coeff)") +
      labs(title = paste0("time to event = ", time, "\n", "status = ", status),
           subtitle = paste("Total =", nrow(res),
                            "| DOWN =", sum(res$diffexpressed == "DOWN", na.rm = TRUE),
                            "| UP =", sum(res$diffexpressed == "UP", na.rm = TRUE)),
           caption = paste(status, " sample sizes | ", paste(names(tblGroup), tblGroup, sep = " = ", collapse = ", "), "\n",
                           "p-value threshold: ", threshold_pval))

      pdf(file = paste(output_folder, paste0("survivalAnalysis_", time, "_volcanoplot.pdf"), sep = "/"), width = 8)
      print(volcano_plot)
      dev.off()
      
      pdf(file = paste(output_folder, paste0("survivalAnalysis_", time, "_exp_coeff_CIplot.pdf"), sep = "/"), width = 8)
      print(CIplot)
      dev.off()
      
  } # end of if(plot)

} # end of function survContinuous




