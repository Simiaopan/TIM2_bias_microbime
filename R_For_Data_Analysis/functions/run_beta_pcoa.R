#-------------------------------------------------------------------------------
# run beta-diversity based pcoa, get coordinates (2D, 3D, all)
#-------------------------------------------------------------------------------
run_beta_pcoa <- function(ps.obj,
                     dists,
                     color_col = "Treat_ID",
                     shape_col = "TimePoint",
                     path_group = "Batch_ID",
                     time_levels = c("T0", "T2","T4","T8","T24","T48","T72")) {
  
  meta <- as(sample_data(ps.obj), "data.frame")
  need <- c(color_col, shape_col, path_group)
  miss <- setdiff(need, colnames(meta))
  if (length(miss) > 0) {
    stop("Missing columns in metadata: ", paste(miss, collapse = ", "))
  }
  
  results <- list()
  
  for (dist_name in dists) {
    # check tree
    if (dist_name %in% c("unifrac", "wunifrac")) {
      tr <- phy_tree(ps.obj, errorIfNULL = FALSE)
      if (is.null(tr)) {
        warning("Skipping distance ", dist_name, " because phylogenetic tree is missing.")
        next
      }
    }
    
    ord <- tryCatch(
      ordinate(ps.obj, method = "PCoA", distance = dist_name),
      error = function(e) { 
        warning("PCoA failed for ", dist_name, ": ", conditionMessage(e)) 
        return(NULL) 
      }
    )
    if (is.null(ord)) next
    
    # extract all coordinates
    if (!is.null(ord$points)) {
      coords_all <- ord$points
    } else if (!is.null(ord$vectors)) {
      coords_all <- ord$vectors
    } else {
      warning("Skipping distance ", dist_name, " (no usable coordinates found).")
      next
    }
    
    # ensure column names
    colnames(coords_all) <- paste0("Axis", seq_len(ncol(coords_all)))
    
    # subset coords (2D & 3D)
    coords_2 <- coords_all[, 1:2, drop = FALSE]
    coords_3 <- if (ncol(coords_all) >= 3) coords_all[, 1:3, drop = FALSE] else NULL
    
    samp <- rownames(coords_all)
    meta.sub <- meta[samp, , drop = FALSE]
    
    if (shape_col %in% colnames(meta.sub)) {
      meta.sub[[shape_col]] <- factor(meta.sub[[shape_col]],
                                      levels = time_levels,
                                      ordered = TRUE)
    }
    
    # explained variance
    eig <- suppressWarnings({
      if (!is.null(ord$values$Relative_eig)) {
        ord$values$Relative_eig * 100
      } else if (!is.null(ord$values$Rel_corr_eig)) {
        ord$values$Rel_corr_eig * 100
      } else {
        rep(NA_real_, ncol(coords_all))
      }
    })
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
    
    results[[dist_name]] <- list(
      coords_2 = df_2,          # first2 + metadata
      coords_3 = df_3,          # first3 + metadata
      coords_all = coords_all,  # all coordinates
      eig = eig,                # all explained variance
      expl_2 = expl_2,          # variance explained by first2
      dist_name = dist_name
    )
  }
  
  return(results)
}
