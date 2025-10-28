#-------------------------------------------------------------------------------
# Add hierarchical clustering layers to RPCA biplot tables
# - Exclude "Noise" from hclust (ordering & cutree on non-noise centroids)
# - Keep Noise samples in output but set SuperCluster/SubCluster to NA
# - All other behavior mirrors the old function
#-------------------------------------------------------------------------------
add_hclust_to_HDBSCAN <- function(
    RPCA, 
    pcs = paste0("PC", 1:15),
    min_super = 2, 
    min_sub = 3,
    link = c("average", "complete", "ward.D2")
) {
  suppressPackageStartupMessages({ library(dplyr) })
  link <- match.arg(link)
  
  message(">>> Start hierarchical clustering ON HDBSCAN clusters (for each RPCA subset) ...")
  out_list <- RPCA
  
  for (dataset in names(RPCA)) {
    message("\n------------------------------------------------------")
    message("Processing dataset: ", dataset)
    
    # --- basic checks (same as old) ---
    if (!is.list(RPCA[[dataset]]) || is.null(RPCA[[dataset]]$biplot)) {
      warning(dataset, ": no $biplot found, skipped.")
      next
    }
    bp <- RPCA[[dataset]]$biplot
    if (!"Cluster" %in% names(bp)) {
      warning(dataset, ": no 'Cluster' column in biplot, skipped.")
      next
    }
    pcs_exist <- pcs[pcs %in% names(bp)]
    if (length(pcs_exist) == 0) {
      warning(dataset, ": no requested PC columns present, skipped.")
      next
    } else if (length(pcs_exist) < length(pcs)) {
      message("  (info) Some PCs missing; using PCs found: ",
              paste(pcs_exist, collapse = ", "))
    }
    
    # --- Step 1: prepare data, exclude Noise for clustering ONLY ---
    bp <- bp %>% mutate(Cluster = as.character(Cluster))
    bp_nonnoise <- bp %>% filter(Cluster != "Noise")
    
    if (n_distinct(bp_nonnoise$Cluster) < 2) {
      warning(dataset, ": <2 non-noise clusters, skipped.")
      next
    }
    
    # --- Step 2: compute centroids on non-noise clusters ---
    centroids <- bp_nonnoise %>%
      group_by(Cluster) %>%
      summarise(across(all_of(pcs_exist), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
    
    # drop invalid rows
    valid_rows <- rowSums(is.finite(as.matrix(centroids[, pcs_exist, drop = FALSE]))) > 0
    centroids <- centroids[valid_rows, , drop = FALSE]
    
    n_clusters <- nrow(centroids)
    if (n_clusters < 2) {
      warning(dataset, ": <2 valid centroids, skipped.")
      next
    }
    
    # --- Step 3: prepare data for clustering ---
    centroids_df <- as.data.frame(centroids, stringsAsFactors = FALSE)
    rownames(centroids_df) <- centroids_df$Cluster
    cent_mat <- as.matrix(centroids_df[, pcs_exist, drop = FALSE])
    
    if (any(!is.finite(cent_mat))) {
      for (j in seq_len(ncol(cent_mat))) {
        col <- cent_mat[, j]
        col[!is.finite(col)] <- 0
        cent_mat[, j] <- col
      }
    }
    
    # --- Step 4: perform hierarchical clustering ---
    message("  Running hclust (method = ", link, ") ...")
    d  <- dist(cent_mat)
    hc <- hclust(d, method = link)
    
    # --- Step 5: define cut levels ---
    k_super <- max(min_super, round(sqrt(n_clusters)))
    k_sub   <- max(min_sub,   round(log2(n_clusters)))
    k_super <- min(max(1, k_super), n_clusters - 1)
    k_sub   <- min(max(1, k_sub),   n_clusters - 1)
    
    super_groups <- cutree(hc, k = k_super)
    sub_groups   <- cutree(hc, k = k_sub)
    
    if (is.null(names(super_groups)) || length(names(super_groups)) == 0)
      names(super_groups) <- rownames(cent_mat)
    if (is.null(names(sub_groups)) || length(names(sub_groups)) == 0)
      names(sub_groups) <- rownames(cent_mat)
    
    cluster_ids <- rownames(cent_mat)
    
    cluster_hierarchy <- data.frame(
      Cluster      = cluster_ids,
      SuperCluster = paste0("G", as.integer(super_groups[cluster_ids])),
      SubCluster   = paste0("S", as.integer(sub_groups[cluster_ids])),
      stringsAsFactors = FALSE
    )
    
    # --- Step 6: merge back to samples ---
    message("  Merging cluster hierarchy back to biplot ...")
    cols_to_drop <- intersect(names(bp),
                              c("SuperCluster", "SubCluster",
                                "SuperCluster.x", "SubCluster.x",
                                "SuperCluster.y", "SubCluster.y"))
    
    message("    Columns to drop: ", paste(cols_to_drop, collapse = ", "))
    message("    select() version: ", as.character(utils::packageVersion("dplyr")))
    
    # make sure select() from dplyr
    message("    select() source: ", attr(get("select", envir = asNamespace("dplyr")), "srcref") %||% "dplyr::select")
    
    bp_out <- tryCatch({
      bp %>%
        dplyr::select(-dplyr::all_of(cols_to_drop)) %>%
        dplyr::left_join(cluster_hierarchy, by = "Cluster") %>%
        dplyr::mutate(
          SuperCluster = if_else(Cluster == "Noise", NA_character_, SuperCluster),
          SubCluster   = if_else(Cluster == "Noise", NA_character_, SubCluster)
        )
    }, error = function(e) {
      message("  ⚠️ Error in select/join/mutate step: ", e$message)
      stop("Failed in dataset ", dataset, ": ", e$message)
    })
    
    # --- Step 7: write back ---
    out_list[[dataset]]$biplot            <- bp_out
    out_list[[dataset]]$hclust_cluster    <- hc
    out_list[[dataset]]$cluster_centroids <- centroids_df
    out_list[[dataset]]$cluster_distance  <- d
    
    message("  ✅ Finished dataset: ", dataset)
  }
  
  message("\n>>> Done processing all RPCA subsets.")
  return(out_list)
}
