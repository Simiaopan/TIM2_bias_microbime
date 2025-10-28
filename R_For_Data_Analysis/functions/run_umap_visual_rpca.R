#-------------------------------------------------------------------------------
# RUN UMAP for visualization only (no HDBSCAN)
#-------------------------------------------------------------------------------
run_umap_visual_rpca <- function(rpca_results,
                                 dataset_name,
                                 n_pcs = 15,          # RPCA输入维度
                                 n_neighbors = 20,
                                 min_dist = 0.15,
                                 n_umap = 2,          # 输出维度（2或3用于可视化）
                                 color_by = "Cluster", # 着色变量，可选
                                 point_size = 2,
                                 alpha = 0.8,
                                 plot = FALSE) {
  # ---- Load packages ----
  library(umap)
  library(ggplot2)
  library(scales)
  
  # ---- Extract RPCA results ----
  if (!dataset_name %in% names(rpca_results)) {
    stop(paste0("Dataset '", dataset_name, "' not found in rpca_results"))
  }
  
  biplot_df <- rpca_results[[dataset_name]]$biplot
  
  # ---- Select PCA coordinates ----
  pc_cols <- grep("^PC[0-9]+$", colnames(biplot_df), value = TRUE)
  if (length(pc_cols) < n_pcs) {
    stop(paste0("Not enough PCs available. Found only ", length(pc_cols)))
  }
  coords <- as.matrix(biplot_df[, pc_cols[1:n_pcs]])
  coords <- scale(coords)
  
  # ---- Step 1: UMAP ----
  message("Running UMAP for visualization ...")
  umap_res <- umap(coords,
                   n_neighbors = n_neighbors,
                   min_dist = min_dist,
                   n_components = n_umap)
  umap_coords <- umap_res$layout
  colnames(umap_coords) <- paste0("UMAP", seq_len(ncol(umap_coords)))
  
  # ---- Save results ----
  rpca_results[[dataset_name]]$umap_visual <- umap_res
  rpca_results[[dataset_name]]$biplot <- cbind(biplot_df, umap_coords)
  
  # ---- Step 2: Visualization ----
  if (plot && n_umap >= 2) {
    if (!color_by %in% colnames(biplot_df)) {
      warning(paste0("Column '", color_by, "' not found in biplot_df. Defaulting to 'grey'."))
      plot_df <- data.frame(umap_coords[, 1:2, drop = FALSE])
      plot_df$ColorVar <- "All"
    } else {
      plot_df <- data.frame(umap_coords[, 1:2, drop = FALSE],
                            ColorVar = biplot_df[[color_by]])
    }
    
    p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = ColorVar)) +
      geom_point(size = point_size, alpha = alpha) +
      scale_color_manual(values = c("Noise" = "grey70",
                                    hue_pal()(length(unique(plot_df$ColorVar))))) +
      theme_bw() +
      labs(title = paste0("UMAP Visualization - ", dataset_name),
           subtitle = paste0("PCA=", n_pcs,
                             ", UMAP_dim=", n_umap,
                             ", color_by=", color_by),
           x = "UMAP1", y = "UMAP2") +
      theme(legend.position = "right")
    
    print(p)
  }
  
  return(rpca_results[[dataset_name]])
}
