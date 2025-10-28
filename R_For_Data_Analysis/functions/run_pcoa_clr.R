#-------------------------------------------------------------------------------
# run CLR based pcoa, get coordinates (2D, 3D, all)
#-------------------------------------------------------------------------------
run_pcoa_clr <- function(ps.obj,
                         color_col = "Treat_ID",
                         shape_col = "TimePoint",
                         path_group = "Batch_ID",
                         time_levels = c("T0","T2","T4","T8","T24","T48","T72")) {
  
  require(phyloseq)
  require(vegan)
  require(dplyr)
  
  meta <- as(sample_data(ps.obj), "data.frame")
  need <- c(color_col, shape_col, path_group)
  miss <- setdiff(need, colnames(meta))
  if (length(miss) > 0) {
    stop("Missing columns in metadata: ", paste(miss, collapse = ", "))
  }
  
  otu <- as(otu_table(ps.obj), "matrix")
  if (taxa_are_rows(ps.obj)) {
    otu <- t(otu)
  }
  
  # Euclidean dists = Aitchison dists
  dist_clr <- dist(otu, method = "euclidean")
  
  # PCoA
  ord <- cmdscale(dist_clr, k = nrow(otu) - 1, eig = TRUE)
  
  coords_all <- ord$points
  colnames(coords_all) <- paste0("Axis", seq_len(ncol(coords_all)))
  
  coords_2 <- coords_all[, 1:2, drop = FALSE]
  coords_3 <- if (ncol(coords_all) >= 3) coords_all[, 1:3, drop = FALSE] else NULL
  
  samp <- rownames(coords_all)
  meta.sub <- meta[samp, , drop = FALSE]
  
  if (shape_col %in% colnames(meta.sub)) {
    meta.sub[[shape_col]] <- factor(meta.sub[[shape_col]],
                                    levels = time_levels,
                                    ordered = TRUE)
  }
  
  # variance explination
  eig <- ord$eig / sum(ord$eig) * 100
  expl_2 <- eig[1:2]
  
  # combine coords with metadata
  df_2 <- dplyr::bind_cols(as.data.frame(coords_2), meta.sub)
  df_2 <- df_2[order(df_2[[path_group]], df_2[[shape_col]]), , drop = FALSE]
  
  if (!is.null(coords_3)) {
    df_3 <- dplyr::bind_cols(as.data.frame(coords_3), meta.sub)
    df_3 <- df_3[order(df_3[[path_group]], df_3[[shape_col]]), , drop = FALSE]
  } else {
    df_3 <- NULL
  }
  
  # results
  results <- list(
    coords_2   = df_2,          # first2 + metadata
    coords_3   = df_3,          # first3 + metadata
    coords_all = coords_all,    # all coordinates
    eig        = eig,           # all explained variance
    expl_2     = expl_2,        
    dist_name  = "clr_euclidean"
  )
  
  return(results)
}
