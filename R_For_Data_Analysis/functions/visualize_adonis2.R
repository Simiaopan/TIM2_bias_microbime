#-------------------------------------------------------------------------------
# visualize adonis results
#-------------------------------------------------------------------------------
visualize_adonis2 <- function(adonis_res, color_vec = NULL) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
  })
  
  adonis_df <- adonis_res %>%
    as.data.frame() %>%
    dplyr::mutate(Variable = rownames(.)) %>%
    dplyr::filter(!Variable %in% c("Residual", "Total")) %>%
    dplyr::select(Variable, Df, R2, `Pr(>F)`) %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        `Pr(>F)` <= 0.001 ~ "***",
        `Pr(>F)` <= 0.01  ~ "**",
        `Pr(>F)` <= 0.05  ~ "*",
        `Pr(>F)` <= 0.1   ~ ".",
        TRUE              ~ "ns"
      )
    )
  
  var_order <- c("TimePoint", "Project", "TimePoint:Project")
  var_order <- var_order[var_order %in% adonis_df$Variable]
  
  adonis_df <- adonis_df %>%
    dplyr::mutate(Variable = factor(Variable, levels = rev(var_order)))
  
  if (is.null(color_vec)) {
    color_vec <- c("#00BFC4", "#F8766D", "#7CAE00")  # TimePoint, Project, Interaction
  }
  color_vec <- rep(color_vec, length.out = length(levels(adonis_df$Variable)))
  
  p <- ggplot2::ggplot(adonis_df, ggplot2::aes(x = Variable, y = R2, fill = Variable)) +
    ggplot2::geom_col(show.legend = FALSE, width = 0.6) +
    ggplot2::geom_text(
      ggplot2::aes(label = sig, y = R2 + max(R2) * 0.03),
      size = 6,
      fontface = "bold"
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_fill_manual(values = color_vec) +
    ggplot2::scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, max(adonis_df$R2) * 1.15)
    ) +
    ggplot2::labs(
      title = expression("Explained variance ("*R^2*") from PERMANOVA"),
      x = NULL,
      y = expression(R^2)
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.text.y = ggplot2::element_text(face = "bold", size = 13),
      axis.text.x = ggplot2::element_text(size = 12),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 16),
      plot.margin = ggplot2::margin(10, 30, 10, 10)
    )
  
  print(p)
  return(adonis_df)
}

