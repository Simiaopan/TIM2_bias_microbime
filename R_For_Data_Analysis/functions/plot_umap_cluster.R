#-------------------------------------------------------------------------------
# make UMAP-HDBSCAN biplot (sample points + optional cluster boundaries)
#-------------------------------------------------------------------------------
plot_umap_cluster <- function(rpca_results,
                              dataset_name,
                              color_col = "Cluster",
                              shape_col = NULL,
                              path_group = NULL,
                              add_ellipse = FALSE,
                              highlight_noise = FALSE,
                              color_vec = NULL,
                              shape_vec = NULL,
                              title_prefix = "UMAP + HDBSCAN") {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(scales)
  
  # ---- 取出数据 ----
  if (!dataset_name %in% names(rpca_results)) {
    stop(paste0("Dataset '", dataset_name, "' not found in rpca_results"))
  }
  res <- rpca_results[[dataset_name]]
  
  if (is.null(res$umap) || is.null(res$hdbscan_umap)) {
    stop(paste0("UMAP or HDBSCAN results not found in rpca_results[['", dataset_name, "']]."))
  }
  
  df <- res$biplot
  umap_coords <- res$umap$layout
  if (ncol(umap_coords) < 2) {
    stop("UMAP result has less than 2 dimensions — cannot plot 2D projection.")
  }
  colnames(umap_coords) <- paste0("UMAP", seq_len(ncol(umap_coords)))
  plot_df <- df
  plot_df[, colnames(umap_coords)] <- umap_coords
  
  # ---- 基础图层 ----
  p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2))
  
  # 路径连接
  if (!is.null(path_group)) {
    p <- p + geom_path(aes(group = .data[[path_group]], color = .data[[color_col]]),
                       linewidth = 0.5, alpha = 0.6)
  }
  
  # 椭圆（仅对样本数≥10且非噪声）
  if (add_ellipse) {
    ellipse_groups <- plot_df %>%
      group_by(.data[[color_col]]) %>%
      filter(n() >= 10 & .data[[color_col]] != "Noise")
    
    if (nrow(ellipse_groups) > 0) {
      p <- p + stat_ellipse(data = ellipse_groups,
                            aes(group = .data[[color_col]], fill = .data[[color_col]]),
                            geom = "polygon", alpha = 0.15,
                            level = 0.68, type = "t", linewidth = 0.4)
    }
  }
  
  # 噪声与非噪声
  if ("Cluster" %in% names(plot_df)) {
    noise_df <- plot_df %>% filter(Cluster == "Noise")
    non_noise_df <- plot_df %>% filter(Cluster != "Noise")
  } else {
    noise_df <- plot_df[0, ]
    non_noise_df <- plot_df
  }
  
  if (highlight_noise) {
    # 非噪声点
    if (!is.null(shape_col)) {
      p <- p + geom_point(data = non_noise_df,
                          aes(color = .data[[color_col]], shape = .data[[shape_col]]),
                          size = 2.2, alpha = 0.9)
    } else {
      p <- p + geom_point(data = non_noise_df,
                          aes(color = .data[[color_col]]),
                          size = 2.2, alpha = 0.9)
    }
    # 噪声点
    if (nrow(noise_df) > 0) {
      p <- p + geom_point(data = noise_df,
                          color = "grey70", size = 1.8, alpha = 0.6)
    }
  } else {
    # 只绘制非噪声
    if (!is.null(shape_col)) {
      p <- p + geom_point(data = non_noise_df,
                          aes(color = .data[[color_col]], shape = .data[[shape_col]]),
                          size = 2.2, alpha = 0.9)
    } else {
      p <- p + geom_point(data = non_noise_df,
                          aes(color = .data[[color_col]]),
                          size = 2.2, alpha = 0.9)
    }
  }
  
  # ---- 色板和形状 ----
  if (!is.null(color_vec)) {
    p <- p + scale_color_manual(values = color_vec) +
      scale_fill_manual(values = color_vec)
  } else if ("Cluster" %in% names(plot_df)) {
    cluster_levels <- unique(plot_df$Cluster)
    if ("Noise" %in% cluster_levels) {
      cluster_levels <- c("Noise", setdiff(cluster_levels, "Noise"))
    }
    base_colors <- c("grey70", hue_pal()(length(cluster_levels) - 1))
    color_map <- setNames(base_colors, cluster_levels)
    p <- p + scale_color_manual(values = color_map) +
      scale_fill_manual(values = color_map)
  }
  
  if (!is.null(shape_vec) && !is.null(shape_col)) {
    p <- p + scale_shape_manual(values = shape_vec)
  }
  
  # ---- 统一风格 ----
  p <- p +
    coord_equal() +
    theme_bw() +
    labs(title = paste0(title_prefix, " - ", dataset_name),
         x = "UMAP1", y = "UMAP2") +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "right")
  
  return(p)
}