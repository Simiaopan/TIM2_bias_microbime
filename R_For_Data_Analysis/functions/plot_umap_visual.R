#-------------------------------------------------------------------------------
# make pure UMAP biplot (sample points only, no HDBSCAN)
#-------------------------------------------------------------------------------
plot_umap_visual <- function(rpca_results,
                             dataset_name,
                             color_col = NULL,    # Column for color (e.g., Project/Group)
                             shape_col = NULL,    # Column for shape (e.g., TimePoint)
                             path_group = NULL,   # Optional line connection (e.g., time series)
                             add_ellipse = FALSE, # Whether to draw group-wise ellipses
                             color_vec = NULL,    # Custom color palette
                             shape_vec = NULL,    # Custom shape palette
                             title_prefix = "UMAP Visualization - ") {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(scales)
  
  # ---- Input checks ----
  if (!dataset_name %in% names(rpca_results)) {
    stop(paste0("Dataset '", dataset_name, "' not found in rpca_results"))
  }
  res <- rpca_results[[dataset_name]]
  
  if (is.null(res$umap) && is.null(res$umap_visual)) {
    stop(paste0("No UMAP results found in rpca_results[['", dataset_name, "']]."))
  }
  
  # Prefer umap_visual if available
  umap_res <- if (!is.null(res$umap_visual)) res$umap_visual else res$umap
  
  df <- res$biplot
  umap_coords <- umap_res$layout
  if (ncol(umap_coords) < 2) {
    stop("UMAP result has less than 2 dimensions — cannot plot 2D projection.")
  }
  colnames(umap_coords) <- paste0("UMAP", seq_len(ncol(umap_coords)))
  plot_df <- df
  plot_df[, colnames(umap_coords)] <- umap_coords
  
  # ---- Base layer ----
  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2))
  
  # Optional path connection (e.g., time-series or trajectory)
  if (!is.null(path_group)) {
    p <- p + geom_path(aes(group = .data[[path_group]], 
                           color = .data[[color_col]]),
                       linewidth = 0.5, alpha = 0.6)
  }
  
  # Ellipse (only for groups with >= 5 samples)
  if (add_ellipse && !is.null(color_col)) {
    ellipse_groups <- plot_df %>%
      group_by(.data[[color_col]]) %>%
      filter(n() >= 5)
    
    if (nrow(ellipse_groups) > 0) {
      p <- p + stat_ellipse(data = ellipse_groups,
                            aes(group = .data[[color_col]], 
                                fill = .data[[color_col]]),
                            geom = "polygon", alpha = 0.15,
                            level = 0.68, type = "t", linewidth = 0.4)
    }
  }
  
  # ---- Sample points ----
  if (!is.null(color_col) && !is.null(shape_col)) {
    p <- p + geom_point(aes(color = .data[[color_col]], 
                            shape = .data[[shape_col]]),
                        size = 2.2, alpha = 0.9)
  } else if (!is.null(color_col)) {
    p <- p + geom_point(aes(color = .data[[color_col]]),
                        size = 2.2, alpha = 0.9)
  } else {
    p <- p + geom_point(size = 2.2, alpha = 0.9, color = "steelblue")
  }
  
  # ---- Custom color and shape scales ----
  if (!is.null(color_vec) && !is.null(color_col)) {
    p <- p + scale_color_manual(values = color_vec) +
      scale_fill_manual(values = color_vec)
  }
  if (!is.null(shape_vec) && !is.null(shape_col)) {
    p <- p + scale_shape_manual(values = shape_vec)
  }
  
  # ---- Final style ----
  p <- p +
    coord_equal() +
    theme_bw() +
    labs(title = paste0(title_prefix, dataset_name),
         x = "UMAP1", y = "UMAP2") +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "right")
  
  return(p)
}

