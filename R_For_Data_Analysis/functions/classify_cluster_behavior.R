#-------------------------------------------------------------------------------
# Classify cluster-level behavior across time (T0, T24, T48, T72)
# Rules:
# - If n_T0==0 & n_T72>0  -> Emergent
# - If n_T0>0  & n_T72==0 -> Dissolved
# - If n_T0==0 & n_T72==0 -> Other (never-present)
# - Else (T0 & T72 both present):
#     * if all timepoint sizes ~ equal (within epsilon) -> Stable
#     * else compare only T0 vs T72:
#         n_T72 > n_T0 + eps -> Converging
#         n_T72 < n_T0 - eps -> Dissolved
#        #|n_T72 - n_T0| <= eps -> Oscillatory
#-------------------------------------------------------------------------------
classify_cluster_behavior <- function(cluster_long,
                                      time_col   = "TimePoint",
                                      cluster_col= "Cluster",
                                      batch_col  = "Batch_ID",
                                      epsilon    = 0) {
  library(dplyr); library(tidyr)
  
  dat <- cluster_long %>%
    as.data.frame() %>%
    mutate(across(everything(), as.character))
  
  # keep only required columns
  dat <- dat %>%
    dplyr::select(any_of(c(cluster_col, time_col, batch_col))) %>%
    distinct()
  
  # Ensure time ordering
  time_levels <- c("T0","T24","T48","T72")
  dat[[time_col]] <- factor(dat[[time_col]], levels = time_levels)
  
  # Count samples per Cluster × TimePoint
  cnt <- dat %>%
    group_by(.data[[cluster_col]], .data[[time_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    complete(
      !!rlang::sym(cluster_col),
      !!rlang::sym(time_col),
      fill = list(n = 0)
    )
  
  # Wide table
  wide <- cnt %>%
    pivot_wider(names_from = !!rlang::sym(time_col), values_from = n) %>%
    mutate(across(all_of(time_levels), ~ replace_na(., 0)))
  
  # Helper to test equality
  all_equal_eps <- function(x, eps) {
    x <- x[!is.na(x)]
    if (length(x) <= 1) return(TRUE)
    (max(x) - min(x)) <= eps
  }
  
  out <- wide %>%
    rowwise() %>%
    mutate(
      n_T0  = as.numeric(coalesce(`T0`,  0)),
      n_T24 = as.numeric(coalesce(`T24`, 0)),
      n_T48 = as.numeric(coalesce(`T48`, 0)),
      n_T72 = as.numeric(coalesce(`T72`, 0)),
      behavior = {
        if (n_T0 == 0 && n_T72 > 0) {
          "Emergent"
        } else if (n_T0 > 0 && n_T72 == 0) {
          "Dissolved"
        } else if (n_T0 == 0 && n_T72 == 0) {
          "Intermediate"
        } else {
          sizes <- c(n_T0, n_T24, n_T48, n_T72)
          if (all_equal_eps(sizes, epsilon)) {
            "Stable"
          } else if (n_T72 > n_T0 + epsilon) {
            "Converging"
          } else if (n_T72 < n_T0 - epsilon) {
            "Dissolved"
          } else {
            "Oscillatory"
          }
        }
      }
    ) %>%
    ungroup() %>%
    dplyr::select(all_of(c(cluster_col, "n_T0", "n_T24", "n_T48", "n_T72", "behavior"))) %>%
    rename_with(~ "Cluster", all_of(cluster_col))
  
  out
}
