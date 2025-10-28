#-------------------------------------------------------------------------------
# make RPCA biplot (sample points + species arrow)
#-------------------------------------------------------------------------------
plot_rpca_biplot <- function(rpca_results,
                             dataset_name,
                             color_col = "Project",
                             shape_col = "TimePoint",
                             path_group = NULL,
                             add_species = TRUE,
                             top_n_species = 20,
                             add_ellipse = TRUE,
                             color_vec = NULL,
                             shape_vec = NULL) {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  
  res <- rpca_results[[dataset_name]]
  df <- res$biplot
  eigvals <- res$eigvals
  prop <- res$proportion
  species_df <- res$species
  
  # ---- choose top species ----
  if (add_species) {
    species_df <- as.data.frame(species_df) %>%
      mutate(score = abs(PC1) + abs(PC2)) %>%
      arrange(desc(score)) %>%
      slice_head(n = top_n_species)
  }
  
  # ---- axis label ----
  xlab_txt <- if (!is.null(prop)) sprintf("PC1 [%.1f%%]", prop$PC1 * 100) else "PC1"
  ylab_txt <- if (!is.null(prop)) sprintf("PC2 [%.1f%%]", prop$PC2 * 100) else "PC2"
  
  p <- ggplot(df, aes(x = PC1, y = PC2))
  
  # path
  if (!is.null(path_group)) {
    p <- p + geom_path(aes(group = .data[[path_group]], color = .data[[color_col]]),
                       linewidth = 0.5, alpha = 0.6)
  }
  
  # ellipse
  if (add_ellipse) {
    p <- p + stat_ellipse(aes(group = .data[[color_col]], fill = .data[[color_col]]),
                          geom = "polygon", alpha = 0.15,
                          level = 0.90, type = "t", linewidth = 0.4)
  }
  
  # points
  if (!is.null(shape_col)) {
    p <- p + geom_point(aes(color = .data[[color_col]], shape = .data[[shape_col]]),
                        size = 2.2, alpha = 0.9)
  } else {
    p <- p + geom_point(aes(color = .data[[color_col]]),
                        size = 2.2, alpha = 0.9)
  }
  
  # ---- species arrow and name ----
  if (add_species) {
    p <- p +
      geom_segment(
        data = species_df,
        aes(x = 0, y = 0, xend = PC1, yend = PC2),
        arrow = arrow(length = unit(0.15, "cm")),
        inherit.aes = FALSE, color = "grey40", alpha = 0.7
      ) +
      geom_text_repel(
        data = species_df,
        aes(x = PC1, y = PC2, label = FeatureID),
        inherit.aes = FALSE,
        size = 3.5, color = "grey20",
        box.padding = 0.6,
        point.padding = 0.4,
        max.overlaps = 100,
        segment.color = "grey60"
      )
  }
  
  # ---- extent teh plot----
  p <- p +
    expand_limits(
      x = c(min(df$PC1, species_df$PC1, na.rm = TRUE) * 1.3,
            max(df$PC1, species_df$PC1, na.rm = TRUE) * 1.3),
      y = c(min(df$PC2, species_df$PC2, na.rm = TRUE) * 1.3,
            max(df$PC2, species_df$PC2, na.rm = TRUE) * 1.3)
    )
  
  p <- p +
    coord_equal() +
    theme_bw() +
    labs(title = paste0("RPCA - ", dataset_name),
         x = xlab_txt, y = ylab_txt) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12, face = "bold"))
  
  if (!is.null(color_vec)) {
    p <- p + scale_color_manual(values = color_vec) +
      scale_fill_manual(values = color_vec)
  }
  if (!is.null(shape_vec) && !is.null(shape_col)) {
    p <- p + scale_shape_manual(values = shape_vec)
  }
  
  return(p)
}
