#-------------------------------------------------------------------------------
# PLOT CLUSTER WITH HIERARCHY
#-------------------------------------------------------------------------------
plot_cluster_hierarchy <- function(
    biplot_tbl,
    cluster_long,
    fill_col = "Project",
    fill_colors = NULL,
    title = "Cluster composition and sample size per Cluster",
    link = c("average","complete","ward.D2"),
    hc = NULL   # ← support hclust object or character vector
) {
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  link <- match.arg(link)
  
  #------------------- Data preparation -------------------
  cluster_layers <- biplot_tbl %>%
    dplyr::mutate(Cluster = as.character(Cluster)) %>%
    dplyr::distinct(Cluster, SuperCluster, SubCluster) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  
  counts <- cluster_long %>%
    as.data.frame() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    dplyr::count(Cluster, .data[[fill_col]], name = "n") %>%
    dplyr::rename(FillGroup = .data[[fill_col]]) %>%
    dplyr::inner_join(cluster_layers, by = "Cluster")
  
  counts <- counts %>%
    dplyr::mutate(
      SuperCluster = dplyr::if_else(is.na(SuperCluster), "Noise", SuperCluster),
      SubCluster   = dplyr::if_else(is.na(SubCluster),   "Noise", SubCluster)
    )
  
  #------------------- Cluster order based on centroid similarity -------------------
  existing_clusters <- unique(counts$Cluster)
  
  get_leaf_order <- function() {
    if (inherits(hc, "hclust")) {
      if (is.null(hc$labels)) {
        pcs <- grep("^PC\\d+$", names(biplot_tbl), value = TRUE)
        centroids <- biplot_tbl %>%
          dplyr::mutate(Cluster = as.character(Cluster)) %>%
          dplyr::filter(Cluster != "Noise") %>%
          dplyr::group_by(Cluster) %>%
          dplyr::summarise(dplyr::across(dplyr::all_of(pcs), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
        labs <- centroids$Cluster
      } else {
        labs <- hc$labels
      }
      ord <- hc$order
      return(as.character(labs[ord]))
    }
    if (is.character(hc)) return(hc)
    
    pcs <- grep("^PC\\d+$", names(biplot_tbl), value = TRUE)
    if (length(pcs) == 0) return(character(0))
    
    centroids <- biplot_tbl %>%
      dplyr::mutate(Cluster = as.character(Cluster)) %>%
      dplyr::filter(Cluster != "Noise") %>%
      dplyr::group_by(Cluster) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(pcs), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
    
    if (nrow(centroids) < 2) return(centroids$Cluster)
    
    mat <- as.matrix(centroids[, pcs, drop = FALSE])
    mat[!is.finite(mat)] <- 0
    rownames(mat) <- centroids$Cluster
    hc2 <- stats::hclust(stats::dist(mat), method = link)
    as.character(hc2$labels[hc2$order])
  }
  
  ord_raw <- get_leaf_order()
  ord <- intersect(ord_raw, existing_clusters)
  cluster_levels <- ord
  if ("Noise" %in% existing_clusters) {
    cluster_levels <- c(cluster_levels, "GAP", "Noise")
  }
  
  #------------------- Insert GAP for visual separation -------------------
  counts_gap <- counts %>%
    tibble::add_row(Cluster = "GAP", FillGroup = NA, n = 0,
                    SuperCluster = NA, SubCluster = NA) %>%
    dplyr::mutate(Cluster = factor(Cluster, levels = cluster_levels))
  
  #------------------- Main bar plot -------------------
  p_bar <- ggplot2::ggplot(counts_gap, ggplot2::aes(x = Cluster, y = n, fill = FillGroup)) +
    ggplot2::geom_col(color = "grey30") +
    ggplot2::scale_fill_manual(values = fill_colors, na.value = "white") +
    ggplot2::labs(
      title = title,
      x = "Cluster (grouped by hierarchical structure)",
      y = "Number of samples", fill = fill_col
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, face = "bold", size = 13),
      legend.position = "right",
      legend.key.size = grid::unit(0.6, "lines"),
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11)
    )
  
  #------------------- Bottom bands for SubCluster -------------------
  band_sub <- counts_gap %>%
    dplyr::distinct(Cluster, SubCluster) %>%
    dplyr::mutate(SubCluster = ifelse(Cluster == "Noise", NA, SubCluster))  # Noise 填白
  
  p_band_sub <- ggplot2::ggplot(band_sub, ggplot2::aes(x = Cluster, y = 1, fill = SubCluster)) +
    ggplot2::geom_tile(color = "white", size = 0.3) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_brewer(palette = "Pastel1", na.value = "white", na.translate = FALSE) +
    ggplot2::labs(fill = "SubCluster") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.5, "lines"),
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 10),
      legend.margin = ggplot2::margin(t = -3),
      plot.margin = ggplot2::margin(-12, 5, -15, 5)
    )
  
  #------------------- Bottom bands for SuperCluster -------------------
  band_super <- counts_gap %>%
    dplyr::distinct(Cluster, SuperCluster) %>%
    dplyr::mutate(SuperCluster = ifelse(Cluster == "Noise", NA, SuperCluster))  # Noise 填白
  
  p_band_super <- ggplot2::ggplot(band_super, ggplot2::aes(x = Cluster, y = 1, fill = SuperCluster)) +
    ggplot2::geom_tile(color = "white", size = 0.3) +
    ggplot2::theme_void() +
    ggplot2::scale_fill_brewer(palette = "Set2", na.value = "white", na.translate = FALSE) +
    ggplot2::labs(fill = "SuperCluster") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.5, "lines"),
      legend.text = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_text(size = 10),
      legend.margin = ggplot2::margin(t = -3),
      plot.margin = ggplot2::margin(-18, 5, 0, 5)
    )
  
  #------------------- Combine all panels -------------------
  p_final <- p_bar / p_band_sub / p_band_super + patchwork::plot_layout(heights = c(5, 0.4, 0.4))
  return(p_final)
}
