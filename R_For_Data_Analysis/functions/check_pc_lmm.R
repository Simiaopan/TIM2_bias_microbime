#-------------------------------------------------------------------------------
# check the lmm model on pc
#-------------------------------------------------------------------------------
check_pc_lmm <- function(data, 
                                  pc_names = c("PC1", "PC2", "PC3", "PC4", "PC5"),
                                  time_var = "Time",
                                  project_var = "Project",
                                  random_var = "Batch_ID") {
  library(lme4)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  
  plots_list <- list()
  
  for (pc in pc_names) {
    message("Checking model for ", pc, " ...")
    
    # Build linear mixed model
    fmla <- as.formula(paste0(pc, " ~ ", time_var, " * ", project_var, " + (1|", random_var, ")"))
    fit <- lmer(fmla, data = data, REML = TRUE)
    
    # Extract residuals and fitted values
    resid_vals <- resid(fit)
    fitted_vals <- fitted(fit)
    
    # (1) Residuals vs Fitted plot
    p1 <- ggplot(data.frame(Fitted = fitted_vals, Residuals = resid_vals),
                 aes(x = Fitted, y = Residuals)) +
      geom_point(alpha = 0.7, color = "#2A9D8F") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      theme_bw(base_size = 11) +
      labs(title = paste0(pc, " – Residuals vs Fitted"),
           x = "Fitted values", y = "Residuals") +
      theme(plot.title = element_text(face = "bold", size = 12),
            panel.grid = element_blank())
    
    # (2) Normal Q-Q plot
    qq_data <- data.frame(
      Theoretical = qqnorm(resid_vals, plot.it = FALSE)$x,
      Sample = qqnorm(resid_vals, plot.it = FALSE)$y
    )
    p2 <- ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
      geom_point(alpha = 0.7, color = "#E76F51") +
      geom_abline(slope = 1, intercept = 0, color = "gray40", linetype = "dashed") +
      theme_bw(base_size = 11) +
      labs(title = paste0(pc, " – Normal Q-Q"),
           x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme(plot.title = element_text(face = "bold", size = 12),
            panel.grid = element_blank())
    
    # Combine the two plots into one row
    row_plot <- p1 + p2 + plot_layout(ncol = 2)
    plots_list[[pc]] <- row_plot
  }
  
  # Combine all PCs into one full-page plot (each PC as one row)
  final_plot <- wrap_plots(plots_list, ncol = 1)
  return(final_plot)
}
