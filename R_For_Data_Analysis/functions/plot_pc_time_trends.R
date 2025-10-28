plot_pc_time_trends <- function(res_time_consistency, 
                                plot_pc_title = "Project-specific time slopes across principal components") {
  library(ggplot2)
  library(dplyr)
  
  # Prepare data: add significance symbols
  slopes_plot <- res_time_consistency$slopes %>%
    mutate(
      Signif = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        p.value < 0.1   ~ ".",
        TRUE ~ "ns"
      )
    )
  
  # Fixed color mapping
  manual_cols <- c("#2A9D8F", "#8AB17D", "#E9C46A", "#E76F51", "#6D597A")
  pcs_all <- sort(unique(slopes_plot$PC))
  cols_use <- manual_cols[seq_along(pcs_all)]
  names(cols_use) <- pcs_all
  
  # Plot
  p <- ggplot(slopes_plot,
              aes(x = Project, y = Time.trend, color = PC)) +
    # Error bars (bottom layer)
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  position = position_dodge(width = 0.6),
                  width = 0.25,
                  linewidth = 0.5,
                  color = "gray60") +
    # Points (upper layer)
    geom_point(position = position_dodge(width = 0.6), 
               size = 2.6, alpha = 0.9) +
    # Significance symbols
    geom_text(aes(label = Signif),
              position = position_dodge(width = 0.6),
              hjust = -0.3, vjust = -0.6,
              size = 3.5, color = "black") +
    # Zero reference line
    geom_hline(yintercept = 0,
               linetype = "dashed",
               color = "gray40",
               linewidth = 0.4) +
    coord_flip() +
    facet_wrap(~ PC, scales = "free_y", ncol = 3) +
    scale_color_manual(values = cols_use) +
    theme_bw(base_size = 13) +
    theme(
      legend.position   = "none",
      strip.text        = element_text(size = 13, face = "bold", color = "gray20"),
      strip.background  = element_blank(),   # Remove background and border of facet titles
      panel.grid.major  = element_blank(),   # Remove major grid lines
      panel.grid.minor  = element_blank(),   # Remove minor grid lines
      panel.border      = element_rect(color = "gray70", fill = NA, linewidth = 0.5),
      panel.spacing.x   = unit(0.5, "cm"),
      panel.spacing.y   = unit(1.0, "cm"),   # Increase facet spacing
      aspect.ratio      = 0.9,               # Adjust vertical aspect ratio
      axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1)  
    ) +
    labs(
      x = "Project",
      y = "Estimated time slope (PC score / time)",
      title = plot_pc_title
    )
  
  return(p)
}
