# function to perform survival analysis and 
# generate KM curves for multiple (>1) groups; p-value based on logrank test
# input data: a data frame with time, event and sample grouping columns

generate_KM_curves <- function(data, time_column, event_column, grouping_cols, output_prefix = "results") {

  # Install necessary packages if not already installed
  if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
  if (!requireNamespace("survminer", quietly = TRUE)) install.packages("survminer")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("dplyr")
  
  library(survival)
  library(survminer)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  
  # ensure time and event columns are numeric
  data[[time_column]] <- as.numeric(data[[time_column]])
  data[[event_column]] <- as.numeric(data[[event_column]])
  
  for (group_col in grouping_cols) {
    cat(paste("Survival analysis started for '", group_col, "'\n", sep = ""))
    
    # check if the grouping column exists
    if (!(group_col %in% names(data))) {
      cat(paste("Warning: Grouping column '", group_col, "' not found in the data.\n", sep = ""))
      next # halts the current iteration and advance the for loop index
    }
    
    # if not a factor, make the group_col a factor
    if (!is.factor(data[[group_col]])) {data[[group_col]] <- factor(data[[group_col]])
    cat(paste("Notes: Grouping column '", group_col, "' is a factor now.\n", sep = ""))}
    
    # output 
    output_folder <- paste(output_prefix, group_col, sep = "/")
    if (!dir.exists(output_folder)) dir.create(output_folder, recursive = T)
    
    # define survival formula
    surv_formula <- as.formula(paste0("Surv(", time_column, ", ", event_column, ") ~ ", group_col))
    
    # Cox PH model
    coxph_model <- coxph(surv_formula, data = data)
    coxph_summary <- summary(coxph_model)
    logrank_test <- round(coxph_summary$sctest["test"], 2)
    logrank_p <- format(coxph_summary$sctest["pvalue"], scientific = TRUE, digits = 2)
    
    # console
    paste(str(data))
    print(surv_formula)
    print(coxph_model)
    print(coxph_summary)
    
    
    sink(file = paste(output_folder, paste0(group_col, "_coxph_summary.txt"), sep = "/"))
    paste(str(data))
    print(surv_formula)
    print(coxph_model)
    print(coxph_summary)
    sink()

    # Kaplan-Meier plot
    fit <- do.call("survfit", list(formula = surv_formula, data = quote(data)))
    plot_KM <- ggsurvplot(fit, data = data,
                          pval = paste0("Log-rank = ", logrank_test, "\n", "(p = ", logrank_p, ")"),
                          pval.size = 3,
                          risk.table = TRUE, risk.table.fontsize = 3,
                          xlab = "Time", ylab = "Survival Probability",
                          title = paste("Kaplan-Meier Plot - ", group_col),
                          ggtheme = theme_minimal(10))
    
    
    plot_KM_CI <- ggsurvplot(fit, data = data,
                             pval = paste0("Log-rank = ", logrank_test, "\n", "(p = ", logrank_p, ")"),
                             pval.size = 3,
                             risk.table = TRUE, risk.table.fontsize = 3,
                             xlab = "Time", ylab = "Survival Probability",
                             title = paste("Kaplan-Meier Plot - ", group_col),
                             ggtheme = theme_minimal(10),
                             conf.int = TRUE) # Add this line to include 95% CI

    write.csv(plot_KM$data.survplot, file = paste(output_folder, paste0(group_col, "_KM_survplot.csv"), sep = "/"))
    write.csv(plot_KM$data.survtable, file = paste(output_folder, paste0(group_col, "_KM_survtable.csv"), sep = "/"))
    
    pdf(file = paste(output_folder, paste0(group_col, "_KM.pdf"), sep = "/"), width = 5, height = 5, useDingbats = FALSE)
    grid.arrange(plot_KM$plot, plot_KM$table, nrow = 2, heights = c(0.7, 0.3))
    grid.arrange(plot_KM_CI$plot, plot_KM_CI$table, nrow = 2, heights = c(0.7, 0.3))
    dev.off()
    
    cat(paste("Survival analysis completed for '", group_col, "'\n", sep = ""))
  } # end of for in group_col
} # end of function survival_multipleGroups

