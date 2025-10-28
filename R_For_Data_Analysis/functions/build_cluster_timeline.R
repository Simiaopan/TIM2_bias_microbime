#-------------------------------------------------------------------------------
# Build_cluster_timeline (safe): keep metadata for viz, but DO NOT let them affect analysis
# - Behavior/analysis uses ONLY (Batch_ID, TimePoint, Cluster)
# - Project is attached per Batch_ID (most frequent if conflicting)
# - RUN is attached per (Batch_ID, TimePoint) using a deterministic rule:
#     pick the most frequent RUN; if tie, pick lexicographically smallest
#-------------------------------------------------------------------------------
build_cluster_timeline <- function(biplot_df, 
                                   batch_col   = "Batch_ID", 
                                   time_col    = "TimePoint", 
                                   cluster_col = "Cluster",
                                   meta_cols   = c("Project", "RUN")) {
  
  # 0) Coerce to data.frame and character to avoid factor/list quirks
  df <- biplot_df %>%
    as.data.frame() %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  
  # Debug: 检查核心三列的类型
  cat("\n--- Debug check: column structure before counting ---\n")
  print(sapply(df[c(batch_col, time_col, cluster_col)], class))
  print(sapply(df[c(batch_col, time_col, cluster_col)], is.list))
  cat("--- End of check ---\n\n")
  
  # 1) Clean & order time labels (T0, T24, T48, T72)
  #    Order by extracted numbers if present; otherwise keep unique order
  time_vals  <- unique(df[[time_col]])
  time_nums  <- suppressWarnings(as.integer(gsub("\\D+", "", time_vals)))
  if (any(!is.na(time_nums))) {
    ord <- order(time_nums, na.last = NA)
    time_levels <- time_vals[ord]
  } else {
    time_levels <- sort(time_vals)
  }
  df[[time_col]] <- factor(df[[time_col]], levels = time_levels)
  
  # 2) Build the CORE long table: (Batch_ID, TimePoint, Cluster) ONLY
  core_long <- df %>%
    dplyr::select(dplyr::all_of(c(batch_col, time_col, cluster_col))) %>%
    dplyr::distinct()
  
  # Sanity check: ensure uniqueness per (Batch, Time)
  dup_check <- core_long %>%
    dplyr::count(.data[[batch_col]], .data[[time_col]]) %>%
    dplyr::filter(n > 1)
  if (nrow(dup_check) > 0) {
    warning("Detected multiple Cluster assignments for some (Batch_ID, TimePoint). ",
            "Keeping the first by default. Check your input if this is unexpected.")
    core_long <- core_long %>%
      dplyr::group_by(.data[[batch_col]], .data[[time_col]]) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }
  
  # 3) Attach PROJECT (per Batch_ID) if present in meta_cols
  attach_project <- any(meta_cols %in% "Project") && ("Project" %in% names(df))
  if (attach_project) {
    proj_map <- df %>%
      dplyr::select(dplyr::all_of(c(batch_col, "Project"))) %>%
      dplyr::filter(!is.na(.data[["Project"]])) %>%
      dplyr::count(.data[[batch_col]], .data[["Project"]], name = "n") %>%
      dplyr::group_by(.data[[batch_col]]) %>%
      dplyr::arrange(dplyr::desc(n), .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!batch_col, Project)
    
    # Warn if a Batch_ID maps to multiple Projects
    proj_conflict <- df %>%
      dplyr::distinct(.data[[batch_col]], .data[["Project"]]) %>%
      dplyr::count(.data[[batch_col]]) %>%
      dplyr::filter(n > 1)
    if (nrow(proj_conflict) > 0) {
      warning("Some Batch_IDs map to multiple Projects. Using the most frequent per Batch_ID.")
    }
    
    core_long <- core_long %>%
      dplyr::left_join(proj_map, by = setNames(batch_col, batch_col))
  }
  
  # 4) Attach RUN (per Batch_ID × TimePoint) if present in meta_cols
  attach_run <- any(meta_cols %in% "RUN") && ("RUN" %in% names(df))
  if (attach_run) {
    run_map <- df %>%
      dplyr::select(dplyr::all_of(c(batch_col, time_col, "RUN"))) %>%
      dplyr::filter(!is.na(.data[["RUN"]])) %>%
      dplyr::count(.data[[batch_col]], .data[[time_col]], .data[["RUN"]], name = "n") %>%
      dplyr::group_by(.data[[batch_col]], .data[[time_col]]) %>%
      dplyr::arrange(dplyr::desc(n), .data[["RUN"]], .by_group = TRUE) %>%  # mode then lexicographically smallest
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(!!batch_col, !!time_col, RUN)
    
    core_long <- core_long %>%
      dplyr::left_join(run_map, by = setNames(c(batch_col, time_col), c(batch_col, time_col)))
  }
  
  # 5) The WIDE table for analysis MUST be built from CORE ONLY (no meta)
  cluster_wide <- core_long %>%
    dplyr::select(dplyr::all_of(c(batch_col, time_col, cluster_col))) %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(time_col), values_from = !!rlang::sym(cluster_col))
  
  # 6) Return
  out <- list(
    long = core_long,   # includes optional Project/RUN for viz; unique per (Batch, Time)
    wide = cluster_wide # analysis-ready (only batch × time → cluster)
  )
  message("Cluster timeline created safely: $wide unaffected by metadata; $long carries Project/RUN for viz.")
  return(out)
}
