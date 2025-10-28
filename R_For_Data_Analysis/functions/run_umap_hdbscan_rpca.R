#-------------------------------------------------------------------------------
# RUN UMAP and HDBSCAN
#-------------------------------------------------------------------------------
run_umap_hdbscan_rpca <- function(rpca_results, 
                                  dataset_name,
                                  n_pcs = 15,             # input dimension
                                  n_neighbors = 20,
                                  min_dist = 0.15,
                                  n_umap = 6,           # UMAP output dim
                                  n_hdbscan_dim = 6,     # dim for HDBSCAN (≤ n_umap)
                                  minPts = 8,
                                  plot = FALSE) {
  # ---- Load packages ----
  library(umap)
  library(dbscan)
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
  message("Running UMAP ...")
  umap_res <- umap(coords,
                   n_neighbors = n_neighbors,
                   min_dist = min_dist,
                   n_components = n_umap)
  umap_coords <- umap_res$layout
  colnames(umap_coords) <- paste0("UMAP", seq_len(ncol(umap_coords)))
  
  # ---- Step 2: HDBSCAN ----
  if (n_hdbscan_dim > n_umap) {
    stop("n_hdbscan_dim cannot be greater than n_umap")
  }
  message(paste0("Running HDBSCAN on first ", n_hdbscan_dim, " UMAP dimensions ..."))
  hdb <- hdbscan(umap_coords[, 1:n_hdbscan_dim, drop = FALSE], minPts = minPts)
  cluster_labels <- ifelse(hdb$cluster == 0, "Noise", as.character(hdb$cluster))
  biplot_df$Cluster <- factor(cluster_labels)
  
  # ---- save results ----
  rpca_results[[dataset_name]]$umap <- umap_res
  rpca_results[[dataset_name]]$hdbscan_umap <- hdb
  rpca_results[[dataset_name]]$biplot <- cbind(biplot_df, umap_coords)
  
  # ---- Step 3: calculate info: noise propotion ----
  noise_ratio <- sum(hdb$cluster == 0) / length(hdb$cluster)
  message(paste0("Noise ratio: ", round(noise_ratio * 100, 2), "% (minPts = ", minPts, ")"))
  message(paste0("Detected clusters: ", length(unique(hdb$cluster[hdb$cluster != 0]))))
  
  # ---- optional: plot ----
  if (plot && n_umap >= 2) {
    plot_df <- data.frame(umap_coords[, 1:2, drop = FALSE],
                          Cluster = biplot_df$Cluster)
    p <- ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(values = c("Noise" = "grey70",
                                    hue_pal()(length(unique(hdb$cluster))))) +
      theme_bw() +
      labs(title = paste0("UMAP + HDBSCAN - ", dataset_name),
           subtitle = paste0("PCA=", n_pcs, 
                             ", UMAP=", n_umap,
                             ", HDBSCAN_dim=", n_hdbscan_dim,
                             ", minPts=", minPts),
           x = "UMAP1", y = "UMAP2") +
      theme(legend.position = "right")
    print(p)
  }
  
  return(rpca_results[[dataset_name]])
}
