#-------------------------------------------------------------------------------
# Run HDBSCAN on RPCA results
#-------------------------------------------------------------------------------
run_hdbscan_rpca <- function(rpca_results, 
                             dataset_name,
                             n_pcs = 3,
                             minPts = 5,
                             plot = FALSE) {
  library(dbscan)
  library(ggplot2)
  
  # extract data
  if (!dataset_name %in% names(rpca_results)) {
    stop(paste0("Dataset '", dataset_name, "' not found in rpca_results"))
  }
  
  biplot_df <- rpca_results[[dataset_name]]$biplot
  
  # check
  pc_cols <- grep("^PC[0-9]+$", colnames(biplot_df), value = TRUE)
  if (length(pc_cols) < n_pcs) {
    stop(paste0("Not enough PCs available. Found only ", length(pc_cols)))
  }
  
  coords <- as.matrix(biplot_df[, pc_cols[1:n_pcs]])
  
  # Run HDBSCAN
  hdb <- hdbscan(coords, minPts = minPts)
  biplot_df$Cluster <- factor(ifelse(hdb$cluster == 0, "Noise", hdb$cluster))
  
 
  rpca_results[[dataset_name]]$biplot <- biplot_df
  rpca_results[[dataset_name]]$hdbscan <- hdb
  
  
  if (plot && n_pcs >= 2) {
    p <- ggplot(biplot_df, aes_string(x = pc_cols[1], y = pc_cols[2], color = "Cluster")) +
      geom_point(size = 2, alpha = 0.8) +
      scale_color_manual(values = c("Noise" = "grey70",
                                    scales::hue_pal()(length(unique(hdb$cluster))))) +
      theme_bw() +
      labs(title = paste0("HDBSCAN clustering on RPCA (", dataset_name, ")"),
           subtitle = paste0("minPts = ", minPts, ", PCs = ", n_pcs),
           x = pc_cols[1], y = pc_cols[2]) +
      theme(legend.position = "bottom")
    
    print(p)
  }
  
  return(rpca_results)
}
