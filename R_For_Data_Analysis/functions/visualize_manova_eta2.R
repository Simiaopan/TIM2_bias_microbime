#-------------------------------------------------------------------------------
# visualize MANOVA: need to arrange pillai results by yourself 
#-------------------------------------------------------------------------------
visualize_manova_eta2 <- function(pillai_df, color_vec = NULL) {
  library(ggplot2)
  library(dplyr)
  
  # ---- cal η² ----
  df_plot <- pillai_df %>%
    filter(Term != "(Intercept)") %>%
    mutate(
      eta2 = Stat / (Stat + 1),
      Term = factor(Term, levels = rev(unique(Term)))
    )
  
  # ---- color ----
  if (is.null(color_vec)) {
    color_vec <- c("#2A9D8F", "#E9C46A", "#E76F51") 
  }
  if (length(color_vec) < nrow(df_plot)) {
    color_vec <- rep(color_vec, length.out = nrow(df_plot))
  }
  
  # ---- plot ----
  p <- ggplot(df_plot, aes(x = Term, y = eta2, fill = Term)) +
    geom_col(show.legend = FALSE, width = 0.6) +
    geom_text(
      aes(label = sig, y = eta2 + max(eta2) * 0.03),
      size = 6,
      fontface = "bold"
    ) +
    coord_flip(clip = "off") +  
    scale_fill_manual(values = color_vec) +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, max(df_plot$eta2) * 1.15)
    ) +
    labs(
      title = expression("Multivariate explained variance ("*eta^2*")"),
      x = NULL,
      y = expression(eta^2)
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),  
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      axis.text.y = element_text(face = "bold", size = 13),
      axis.text.x = element_text(size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.margin = margin(10, 30, 10, 10) 
    )
  
  print(p)
  return(df_plot)
}
