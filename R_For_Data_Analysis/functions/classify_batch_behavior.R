#-------------------------------------------------------------------------------
# Classify each Batch_ID into one of four temporal behavior types
#-------------------------------------------------------------------------------
classify_batch_behavior <- function(cluster_wide,
                                    project_col = "Project",
                                    batch_col   = "Batch_ID",
                                    noise_label = "Noise") {
  library(dplyr)
  library(tidyr)
  
  # ---- Clean input ----
  df <- cluster_wide %>%
    as.data.frame() %>%
    mutate(across(everything(), as.character))
  
  # ---- Identify and order time columns (e.g., T0 < T24 < T48 < T72) ----
  time_cols <- grep("^T", names(df), value = TRUE)
  time_nums <- suppressWarnings(as.integer(gsub("\\D+", "", time_cols)))
  if (any(!is.na(time_nums))) time_cols <- time_cols[order(time_nums)]
  
  # ---- Rowwise evaluation ----
  out <- df %>%
    rowwise() %>%
    mutate(
      # Extract the ordered cluster sequence
      clusters_vec = list(c_across(all_of(time_cols))),
      path = paste(clusters_vec, collapse = "→"),
      
      # Extract specific time points (safe for missing ones)
      t0  = clusters_vec[1],
      t24 = clusters_vec[min(2, length(clusters_vec))],
      t48 = clusters_vec[min(3, length(clusters_vec))],
      t72 = clusters_vec[length(clusters_vec)],
      
      # Apply explicit step-by-step classification
      behavior = {
        v <- clusters_vec
        
        # ---- 1) Converging ----
        if (t0 == noise_label && t72 != noise_label) {
          "Converging"
          
          # ---- 2) Diverging ----
        } else if (t0 != noise_label && t72 == noise_label) {
          "Diverging"
          
          # ---- 3) Other (Noise → Noise) ----
        } else if (t0 == noise_label && t72 == noise_label) {
          "Other"
          
          # ---- 4) Stable: T0=T24 且 (T48=T24 或 T72=T24)，且非 Noise ----
        } else if (!is.na(t0) && !is.na(t24) &&
                   t0 == t24 &&
                   ((length(v) >= 3 && t48 == t24) ||
                    (length(v) >= 4 && t72 == t24)) &&
                   t24 != noise_label) {
          "Stable"
          
          # ---- 5) Dynamic transition (remaining cases) ----
        } else {
          "Dynamic transition"
        }
      }
    ) %>%
    ungroup()
  
  # ---- Select relevant columns ----
  keep_cols <- c(batch_col, project_col)
  keep_cols <- keep_cols[keep_cols %in% names(out)]
  
  out %>%
    as.data.frame() %>%
    dplyr::select(all_of(keep_cols), path, behavior)
}
