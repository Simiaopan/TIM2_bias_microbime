#-------------------------------------------------------------------------------
# RPCA biplot（样本点 + 可选物种箭头 + 可选聚类着色/圈边界）
#-------------------------------------------------------------------------------
plot_rpca_biplot_cluster <- function(rpca_results,
                                     dataset_name,
                                     color_col = "Project",
                                     shape_col = "TimePoint",
                                     path_group = NULL,
                                     add_species = TRUE,
                                     top_n_species = 20,
                                     add_ellipse = TRUE,               # 元数据分组的 ellipse
                                     add_cluster = TRUE,               # 是否用 cluster 上色
                                     show_cluster_ellipse = TRUE,      # 是否圈出 cluster 边界
                                     cluster_col = NULL,               # 聚类列名（自动或手动）
                                     color_vec = NULL,                 # 元数据颜色
                                     cluster_color_vec = NULL,         # 聚类颜色
                                     shape_vec = NULL,
                                     ellipse_level = 0.90,
                                     ellipse_exclude = c("Noise")) {
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  
  res <- rpca_results[[dataset_name]]
  df  <- res$biplot
  prop <- res$proportion
  species_df <- res$species
  
  # ---- 找聚类列：自动或手动 ----
  if (is.null(cluster_col)) {
    candidates <- c("Cluster","cluster","HDBSCAN","hdbscan","kmeans","KM","cl","cluster_id")
    cluster_col <- intersect(candidates, colnames(df))
    cluster_col <- if (length(cluster_col)) cluster_col[1] else NULL
  }
  has_cluster <- !is.null(cluster_col) && cluster_col %in% colnames(df)
  
  # ---- 物种箭头（可选） ----
  if (add_species && !is.null(species_df)) {
    species_df <- as.data.frame(species_df) %>%
      mutate(score = abs(PC1) + abs(PC2)) %>%
      arrange(desc(score)) %>%
      slice_head(n = top_n_species)
  } else {
    add_species <- FALSE
  }
  
  # ---- 轴标签 ----
  xlab_txt <- if (!is.null(prop) && "PC1" %in% names(prop)) sprintf("PC1 [%.1f%%]", prop$PC1 * 100) else "PC1"
  ylab_txt <- if (!is.null(prop) && "PC2" %in% names(prop)) sprintf("PC2 [%.1f%%]", prop$PC2 * 100) else "PC2"
  
  # ---- 确定颜色使用逻辑 ----
  color_used <- color_col
  if (add_cluster && has_cluster) {
    color_used <- cluster_col
    message(sprintf("→ Coloring by cluster column: %s", cluster_col))
  } else {
    message(sprintf("→ Coloring by metadata column: %s", color_col))
  }
  
  p <- ggplot(df, aes(x = PC1, y = PC2))
  
  # ---- 时间轨迹 ----
  if (!is.null(path_group)) {
    p <- p + geom_path(aes(group = .data[[path_group]], color = .data[[color_used]]),
                       linewidth = 0.5, alpha = 0.6)
  }
  
  # ---- 元数据分组的置信椭圆 ----
  if (add_ellipse && color_used == color_col) {
    p <- p + stat_ellipse(aes(group = .data[[color_col]], fill = .data[[color_col]]),
                          geom = "polygon", alpha = 0.15,
                          level = ellipse_level, type = "t", linewidth = 0.4)
  }
  
  # ---- cluster 边界圈 ----
  if (show_cluster_ellipse && has_cluster) {
    ell_df <- df
    if (!is.null(ellipse_exclude)) {
      ell_df <- ell_df[!(ell_df[[cluster_col]] %in% ellipse_exclude), , drop = FALSE]
    }
    if (nrow(ell_df) > 0 && length(unique(ell_df[[cluster_col]])) > 0) {
      if (!is.null(cluster_color_vec)) {
        p <- p + stat_ellipse(
          data = ell_df,
          aes(x = PC1, y = PC2, group = .data[[cluster_col]], color = .data[[cluster_col]]),
          geom = "path", linewidth = 0.8, alpha = 0.8,
          type = "t", level = ellipse_level, inherit.aes = FALSE
        ) +
          scale_color_manual(values = cluster_color_vec, na.value = "grey60")
      } else {
        p <- p + stat_ellipse(
          data = ell_df,
          aes(x = PC1, y = PC2, group = .data[[cluster_col]]),
          color = "grey40", geom = "path", linewidth = 0.8, alpha = 0.7,
          type = "t", level = ellipse_level, inherit.aes = FALSE
        )
      }
    }
  }
  
  # ---- 样本点 ----
  if (!is.null(shape_col)) {
    p <- p + geom_point(aes(color = .data[[color_used]], shape = .data[[shape_col]]),
                        size = 2.0, alpha = 0.7)
  } else {
    p <- p + geom_point(aes(color = .data[[color_used]]),
                        size = 2.0, alpha = 0.7)
  }
  
  # ---- 物种箭头与标签 ----
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
        box.padding = 0.6, point.padding = 0.4,
        max.overlaps = 100, segment.color = "grey60"
      )
  }
  
  # ---- 扩边 ----
  x_min <- min(df$PC1, if (add_species) species_df$PC1 else df$PC1, na.rm = TRUE)
  x_max <- max(df$PC1, if (add_species) species_df$PC1 else df$PC1, na.rm = TRUE)
  y_min <- min(df$PC2, if (add_species) species_df$PC2 else df$PC2, na.rm = TRUE)
  y_max <- max(df$PC2, if (add_species) species_df$PC2 else df$PC2, na.rm = TRUE)
  
  p <- p +
    expand_limits(
      x = c(x_min * 1.4, x_max * 1.4),
      y = c(y_min * 1.4, y_max * 1.4)
    ) +
    coord_equal() +
    theme_bw() +
    labs(title = paste0("RPCA - ", dataset_name),
         x = xlab_txt, y = ylab_txt,
         color = color_used) +
    theme(panel.grid = element_blank(),
          plot.title = element_text(size = 12, face = "bold"))
  
  # ---- 手动调色 ----
  if (add_cluster && has_cluster && color_used == cluster_col) {
    if (!is.null(cluster_color_vec)) {
      p <- p + scale_color_manual(values = cluster_color_vec, na.value = "grey70")
    } else {
      lev <- unique(df[[cluster_col]])
      if ("Noise" %in% lev) {
        other_lev <- setdiff(lev, "Noise")
        pal <- scales::hue_pal()(length(other_lev))
        vals <- c("Noise" = "grey70", setNames(pal, other_lev))
        p <- p + scale_color_manual(values = vals, na.value = "grey70")
      } else {
        p <- p + scale_color_hue()
      }
    }
  } else if (!add_cluster && color_used == color_col && !is.null(color_vec)) {
    p <- p + scale_color_manual(values = color_vec, na.value = "grey70") +
      scale_fill_manual(values = color_vec, na.value = "grey70")
  }
  
  if (!is.null(shape_vec) && !is.null(shape_col)) {
    p <- p + scale_shape_manual(values = shape_vec)
  }
  
  return(p)
}
