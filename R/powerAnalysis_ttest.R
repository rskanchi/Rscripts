# function to perform power analysis when sample means, SDs and the observed sample sizes (different) are available. 
# Powere analysis uses these inputs to compute the effect size.
# Power computations are performed for proposed sample sizes.

powerAnalysis_ttest <- function(mean1, mean2, sd1, sd2, n1, n2, proposed_n1, proposed_n2, sig_level = 0.05,
                                output_folder){
  
  dir.create(output_folder, recursive = TRUE)
  sink(paste(output_folder, "log.txt", sep = "/"))
  
  # Write a brief Methods section to the log
  cat("### Methods ###\n")
  cat("This analysis calculates the statistical power and minimum sample size required to detect a 
      difference between two groups.\n")
  cat("The effect size (Cohen's d) is computed using pooled standard deviation, 
      and the power is determined for the proposed sample sizes using a two-sample t-test when the sample sizes for each group are different.\n")
  cat("The actual sample sizes are used for the effect size and pooled standard deviation calculations and proposed sample sizes for the power calculation. \n")
  cat("An ROC analysis is performed using simulated data based on the input parameters (means, SDs, actual sample sizes) and the area under the curve (AUC) is computed.\n")
  cat("An ROC curve is plotted to visualize the classification performance.\n\n")
  cat("The following R packages were used in the analysis:\n")
  cat(" - pwr: for power and sample size calculations\n")
  cat(" - pROC: for calculating and visualizing ROC curves\n")
  cat(" - ggplot2: for generating publication-quality ROC plots\n\n")
  
  library(pwr) # power calculation
  library(pROC) # for ROC analysis
  library(ggplot2) # plotting
  
  cat("### Inputs ###\n")
  cat("mean1 =",mean1, "mean2 =", mean2, "SD1 =",sd1, "SD2 =", sd2, "n1 =", n1, "n2 =", n2, "proposed n1 =", proposed_n1, "proposed n2 =", proposed_n2, "\n")
  
  # Calculate effect size (Cohen's d)
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  effect_size <- (mean2 - mean1) / pooled_sd
  
  # Power calculation for proposed sample sizes
  power <- pwr.t2n.test(n1 = proposed_n1, n2 = proposed_n2, d = effect_size, sig.level = 0.05)$power
  cat("### Results ###\n")
  cat("Effect size (Cohen's d):", round(effect_size, 2), "\n")
  cat("Power for proposed sample sizes:", round(power, 2), "\n")
  cat("This indicates that there is a", round(power * 100, 0),"% chance of correctly detecting a true difference 
  in means between the two groups if it exists, assuming a significance level of 0.05.\n") 
  cat("Conversely, there is a", round((1 - power) * 100, 0), "% chance of failing to detect
      a true effect (Type II error) with these sample sizes.\n\n")
  
  # Calculate sample size for 80% power
  sample_size <- pwr.t.test(d = effect_size, sig.level = 0.05, power = 0.8, type = "two.sample")$n
  cat("Minimum sample size required per group for 80% power:", ceiling(sample_size), "\n")
  
  # Create a range of power values (e.g., from 0.1 to 0.99)
  power_values <- seq(0.2, 0.99, by = 0.1)
  
  # Calculate the minimum sample size for each power value
  sample_sizes <- sapply(power_values, function(p) {
    pwr.t.test(d = effect_size, sig.level = sig_level, power = p, type = "two.sample")$n
  })
  
  # Create a data frame for plotting
  power_data <- data.frame(
    Power = power_values,
    SampleSize = ceiling(sample_sizes)  # Round up to the nearest integer
  )
  write.csv(power_data, paste(output_folder, "power_min_sample_size.csv", sep = "/"), row.names = FALSE)
  # Create the plot using ggplot2
  power_samplesize_plot <- ggplot(power_data, aes(x = Power, y = SampleSize)) +
    geom_line(color = "black", linewidth = 1, linetype = "dashed") +    # Line plot
    geom_point(color = "red", size = 4) +    # Points on the plot
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +  # x-axis breaks at every 0.1
    labs(title = "Sample Size & Power for Two-Sample t-test",
         x = "Power",
         y = "Minimum Sample Size") +
    theme_minimal(base_size = 12) +           # Clean theme with larger text
    theme(plot.title = element_text(hjust = 0.5)) # Center the title
  
  # ROC analysis
  # Simulate data for groups
  set.seed(12345)
  group_1 <- rnorm(n1, mean = mean1, sd = sd1)
  group_2 <- rnorm(n2, mean = mean2, sd = sd2)
  
  # Combine into a data frame
  data <- data.frame(
    value = c(group_1, group_2),
    group = factor(c(rep("1", n1), rep("2", n2)))
  )
  #data$group <- factor(data$group, levels = c(1,2), )
  
  # Calculate ROC curve
  roc_obj <- roc(data$group, data$value)
  auc_value <- auc(roc_obj)
  perf <- get_AUC_performance(auc_value)
  cat("AUC for classification:", round(auc_value, 2), "\n\n")
  perf <- get_AUC_performance(as.numeric(auc_value))
  cat("The AUC value", round(auc_value, 2), "suggests that the ability of the model to distinguish 
      between the two severity groups based on the simulated data is:", perf, "\n")
  cat("The above AUC classification is based on the following:\n")
  cat("auc_value >= 0.9 ~ Excellent\n")
  cat("auc_value >= 0.8 and < 0.9 ~ Good\n")
  cat("auc_value >= 0.7 and < 0.8 ~ Fair\n")
  cat("auc_value >= 0.5 and < 0.7 ~ Poor\n")
  cat("auc_value < 0.5 ~ Fail\n")
  #  cat("auc_value >= 0.6 and < 0.7 ~ Poor\n")
  #  cat("auc_value >= 0.5 and < 0.6 ~ Fail\n")
  #  cat("auc_value <0.5 ~ Worse than random guessing\n")
  
  # Plot ROC curve
  roc_data <- data.frame(
    FPR = 1 - roc_obj$specificities,  # False Positive Rate
    TPR = roc_obj$sensitivities       # True Positive Rate
  )
  
  # Plot ROC curve
  roc_data <-roc_data[order(roc_data$FPR, roc_data$TPR, decreasing = c(TRUE, FALSE)),]
  roc_plot <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
    geom_line(color = "blue", linewidth = 1) +                # ROC curve
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # Diagonal line
    labs(title = "ROC curve",
         x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    xlim(0,1) + ylim(0,1) +
    theme_minimal(base_size = 12)  
  
  # Add AUC to the plot
  roc_plot <- roc_plot +
    annotate("text", x = 0.6, y = 0.1, label = paste("AUC =", round(auc_value, 2)), size = 5, color = "black")
  
  pdf(file = paste0(output_folder, "/power_analysis_plots.pdf"))
  print(power_samplesize_plot)
  print(roc_plot)
  dev.off()
  
  sink()
} # end of function

# Function to classify AUC performance
get_AUC_performance <- function(auc_value) {
  if (auc_value >= 0.9) performance <-"Excellent" else
    if (auc_value >= 0.8 && auc_value < 0.9) performance <-"Good" else
      if (auc_value >= 0.7 && auc_value < 0.8) performance <-"Fair" else
        if (auc_value >= 0.5 && auc_value < 0.7) performance <-"Poor" else
          if (auc_value < 0.5) performance <-"Fail"
          #  cat("The AUC value", round(auc_value, 2), "suggests that the ability of the model to distinguish 
          #      between the two severity groups based on the simulated data is:", performance, ".\n")
          return(performance)
          
}

