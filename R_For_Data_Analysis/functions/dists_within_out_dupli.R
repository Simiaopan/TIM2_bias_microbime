#-------------------------------------------------------------------------------
# Compute within- vs between- duplication distances (+ Project label)
# Also perform per-timepoint significance test (Wilcoxon rank-sum)
#-------------------------------------------------------------------------------
dists_within_out_dupli <- function(coords_all,
                                   meta,
                                   treat_col   = "Treat_ID",
                                   batch_col   = "Batch_ID",
                                   tp_col      = "TimePoint",
                                   time_col    = "Time",
                                   proj_col    = "Project",      
                                   dist_method = "euclidean") {
  library(dplyr)
  
  #------------------------
  # Align sample metadata
  #------------------------
  stopifnot(all(rownames(coords_all) %in% rownames(meta)))
  meta <- meta[rownames(coords_all), , drop = FALSE]
  
  #------------------------
  # Distance matrix
  #------------------------
  Dmat <- as.matrix(dist(coords_all, method = dist_method))
  
  treat <- meta[[treat_col]]
  batch <- as.character(meta[[batch_col]])   
  tp    <- meta[[tp_col]]
  time  <- meta[[time_col]]
  proj  <- meta[[proj_col]]
  
  res_list <- list()
  
  #------------------------
  # Compute pairwise distances per timepoint
  #------------------------
  for (lev in unique(tp)) {
    idx <- which(tp == lev)
    if (length(idx) < 2) next
    
    subD     <- Dmat[idx, idx, drop = FALSE]
    subTreat <- treat[idx]
    subBatch <- batch[idx]
    subTime  <- time[idx]
    subProj  <- proj[idx]
    
    pairs <- combn(seq_along(idx), 2)
    
    # Within duplication
    mask_within <- subTreat[pairs[1,]] == subTreat[pairs[2,]] &
      subBatch[pairs[1,]] != subBatch[pairs[2,]]
    d_within <- subD[cbind(pairs[1, mask_within], pairs[2, mask_within])]
    
    df_within <- data.frame(
      Distance = d_within,
      Group    = "Within",
      TimePoint= lev,
      Time     = subTime[pairs[1, mask_within]],
      PairID   = paste0(subBatch[pairs[1, mask_within]], "_vs_", subBatch[pairs[2, mask_within]]),
      Project  = subProj[pairs[1, mask_within]]   
    )
    
    # Between duplication
    mask_between <- subTreat[pairs[1,]] != subTreat[pairs[2,]]
    d_between <- subD[cbind(pairs[1, mask_between], pairs[2, mask_between])]
    
    df_between <- data.frame(
      Distance = d_between,
      Group    = "Between",
      TimePoint= lev,
      Time     = subTime[pairs[1, mask_between]],
      PairID   = paste0(subBatch[pairs[1, mask_between]], "_vs_", subBatch[pairs[2, mask_between]]),
      Project  = NA  
    )
    
    res_list[[lev]] <- rbind(df_within, df_between)
  }
  
  # Merge all timepoints
  res_df <- do.call(rbind, res_list)
  rownames(res_df) <- NULL
  
  #------------------------
  # Statistical comparison: Within vs Between at each timepoint
  #------------------------
  res_stats <- res_df %>%
    group_by(TimePoint) %>%
    summarise(
      n_within  = sum(Group == "Within"),
      n_between = sum(Group == "Between"),
      p_value   = tryCatch(
        wilcox.test(Distance[Group == "Within"],
                    Distance[Group == "Between"],
                    exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      signif = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        p_adj < 0.1   ~ ".",
        TRUE          ~ "ns"
      )
    )
  
  #------------------------
  # Return results as a list
  #------------------------
  return(list(
    distances = res_df,   # full pairwise distance table
    stats     = res_stats  # per-timepoint test summary with significance marks
  ))
}
