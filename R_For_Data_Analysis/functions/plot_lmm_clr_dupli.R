#-------------------------------------------------------------------------------
# plot lmm trend of within/out duplication dists
#-------------------------------------------------------------------------------
plot_lmm_clr_dupli <- function(lmm_obj, time_range = c(0,72)) {
  library(ggplot2)
  
  fixed_table <- lmm_obj$fixed
  pred_df     <- lmm_obj$pred
  
  coefs <- fixed_table[fixed_table$Effect == "Time", ]
  
  p <- ggplot(pred_df, aes(x = Time, y = fit, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
    geom_line(size = 1.2) +
    scale_color_manual(values = c("Within" = "#1f78b4", "Without" = "#a6cee3")) +
    scale_fill_manual(values = c("Within" = "#1f78b4", "Without" = "#a6cee3")) +
    theme_bw(base_size = 14) +
    scale_x_continuous(breaks = c(0,24,48,72), limits = c(0,72)) +
    labs(x = "Time (h)", y = "Predicted Distance", color = "Group", fill = "Group") +
    theme(legend.position = "top")
  
  ann_text <- do.call(rbind, lapply(split(fixed_table, fixed_table$Group), function(tab) {
    intercept <- tab$Estimate[tab$Effect == "(Intercept)"]
    slope     <- tab$Estimate[tab$Effect == "Time"]
    pval      <- tab$p[tab$Effect == "Time"]
    y_pos     <- ifelse(tab$Group[1] == "Within", max(pred_df$fit) * 0.8, max(pred_df$fit) * 0.9)
    data.frame(
      Time = 5,
      Distance = y_pos,
      Label = sprintf("%s: y = %.2f + %.3f*Time, p = %.3g", 
                      tab$Group[1], intercept, slope, pval),
      Group = tab$Group[1]
    )
  }))
  
  p <- p + geom_text(data = ann_text, aes(label = Label, color = Group),
                     hjust = 0, vjust = -0.5, size = 4.2, show.legend = FALSE)
  
  return(p)
}

