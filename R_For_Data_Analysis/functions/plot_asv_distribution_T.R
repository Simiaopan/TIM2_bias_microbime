plot_core_asv_T <- function(core_asv_table, subject_col = "Batch_ID", p_title) {
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(rlang)
  
  # Step 1: calculate summed relative abundance per sample
  core_sum_df <- core_asv_table %>%
    group_by(Sample, TimePoint, !!sym(subject_col)) %>%
    summarise(core_sum = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
  
  # Step 2: build model formula dynamically
  fml <- as.formula(paste("core_sum ~ TimePoint + (1|", subject_col, ")"))
  
  # Step 3: fit linear mixed model
  m <- lmer(fml, data = core_sum_df)
  
  # Step 4: extract overall p-value for TimePoint
  anova_res <- anova(m)
  p_val <- anova_res$`Pr(>F)`[1]
  
  # Step 5: estimated marginal means (for plotting)
  em_df <- as.data.frame(emmeans(m, ~ TimePoint))
  
  # 🔹 Step 6: calculate ASV richness per TimePoint
  asv_richness <- core_asv_table %>%
    group_by(TimePoint) %>%
    summarise(n_ASV = n_distinct(OTU), .groups = "drop")
  
  # Step 7: plotting
  p <- ggplot(core_sum_df, aes(x = TimePoint, y = core_sum, fill = TimePoint)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.1, size = 2) +  # scatter transparency
    geom_line(data = em_df, aes(x = TimePoint, y = emmean, group = 1),
              color = "red", size = 1) +
    geom_point(data = em_df, aes(x = TimePoint, y = emmean),
               color = "red", size = 3) +
    geom_text(
      data = asv_richness,
      aes(x = TimePoint, y = max(core_sum_df$core_sum) * 1.03, label = n_ASV),
      size = 4.2,
      fontface = "bold",
      color = "black"
    ) +
    labs(
      title = "Change of core-ASV summed relative abundance over time",
      x = "TimePoint",
      y = "Sum of core-ASV relative abundance per sample"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none") +
    annotate("text",
             x = 0.6,
             y = max(core_sum_df$core_sum) * 1.1,
             label = sprintf("p = %.3g", p_val),
             hjust = 0,
             size = 5,
             color = "red")
  
  # Step 8: print model summary
  cat("---------------------------------------------------\n")
  cat("Linear Mixed Model:\n")
  print(summary(m))
  cat(sprintf("\nOverall p-value for TimePoint = %.4g\n", p_val))
  cat("---------------------------------------------------\n")
  
  print(p)
  
  return(list(model = m, emmeans = em_df, richness = asv_richness, plot = p))
}
