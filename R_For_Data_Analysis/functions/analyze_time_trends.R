#===============================================================================
# analyze_time_trends
#-------------------------------------------------------------------------------
# Inputs:
#   res_df     : MaAsLin2 result table (contains feature, metadata, coef, stderr, qval, etc.)
#   sub_df     : Subset data used for “real FC + project-level statistics” (e.g., asv_raw_sub or genus_raw_sub)
#   feature_col: Column name of the feature in sub_df ("OTU" or "Genus")
#   timepoint_col: Column name for time points in sub_df (default "TimePoint", includes T0/T72)
#   project_col  : Column name for projects in sub_df (default "Project")
#   relab_col    : Column name for relative abundance in sub_df (default "RelAbundance")
#   time_effect_label: Main effect name in res_df (default "Time_c")
#   hours       : Time span for conversion (default 72)
#   alpha       : Confidence level (default 0.05 → 95% CI)
#   eps         : Pseudocount used for real FC calculation (to avoid division by zero)
#
# Returns:
#   list(
#     merged : Summary table with one row per feature (model FC, 95% CI, real FC, project consistency ratio)
#     per_project: Optional table with feature × Project mean values and directions at T0/T72
#   )
#===============================================================================
analyze_time_trends <- function(res_df, sub_df,
                                feature_col = "Genus",
                                timepoint_col = "TimePoint",
                                project_col = "Project",
                                relab_col = "RelAbundance",
                                time_effect_label = "Time_c",
                                hours = 72,
                                alpha = 0.05,
                                eps = 1e-8,
                                fc72_cutoff = 1,
                                fcreal_cutoff = 1,
                                nz_ratio_cutoff = 0.5,
                                consistency_cutoff = 0.7) {
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
  })
  
  # ---------- Safety check ----------
  need_cols <- c(feature_col, timepoint_col, project_col, relab_col)
  miss <- setdiff(need_cols, names(sub_df))
  if (length(miss)) stop(sprintf("sub_df is missing columns: %s", paste(miss, collapse=", ")))
  
  req_res <- c("feature","metadata","coef","stderr","qval","N","N.not.0")
  miss_res <- setdiff(req_res, names(res_df))
  if (length(miss_res)) stop(sprintf("res_df is missing columns: %s", paste(miss_res, collapse=", ")))
  
  z <- qnorm(1 - alpha/2)
  
  # ---------- 1) Model-level (main effect) ----------
  res_main <- res_df %>%
    filter(.data$metadata == time_effect_label) %>%
    transmute(
      feature, coef, stderr, pval, qval, N, N.not.0,
      nz_ratio = ifelse(N > 0, N.not.0 / N, NA_real_),
      coef_low  = coef - z * stderr,
      coef_high = coef + z * stderr,
      FC72        = exp(coef      * hours),
      FC72_low95  = exp(coef_low  * hours),
      FC72_high95 = exp(coef_high * hours),
      log2FC72        = (coef      * hours) / log(2),
      log2FC72_low95  = (coef_low  * hours) / log(2),
      log2FC72_high95 = (coef_high * hours) / log(2),
      ci95_cross_0 = (coef_low <= 0 & coef_high >= 0),
      ci95_cross_1 = (FC72_low95 <= 1 & FC72_high95 >= 1),
      FDR_lmm = qval,
      dir_model = case_when(
        !ci95_cross_1 & FC72 > 1 ~ "up",
        !ci95_cross_1 & FC72 < 1 ~ "down",
        TRUE ~ "nonsignificant"
      )
    )
  
  # ---------- 2) sub_df: overall real log2FC (T0 vs T72) ----------
  sub_core <- sub_df %>%
    select(feature = all_of(feature_col),
           TimePoint = all_of(timepoint_col),
           RelAbundance = all_of(relab_col)) %>%
    mutate(TimePoint = as.character(TimePoint)) %>%
    filter(TimePoint %in% c("T0", "T72"))
  
  real_fc_overall <- sub_core %>%
    group_by(feature, TimePoint) %>%
    summarise(mean_rel = mean(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = TimePoint, values_from = mean_rel) %>%
    mutate(T0 = ifelse(is.na(T0), 0, T0),
           T72 = ifelse(is.na(T72), 0, T72),
           log2FC_real = log2((T72 + eps) / (T0 + eps))) %>%
    select(feature, log2FC_real)
  
  # ---------- 3) sub_df: per-project direction proportion ----------
  per_project <- sub_df %>%
    select(feature = all_of(feature_col),
           Project = all_of(project_col),
           TimePoint = all_of(timepoint_col),
           RelAbundance = all_of(relab_col)) %>%
    mutate(TimePoint = as.character(TimePoint)) %>%
    filter(TimePoint %in% c("T0", "T72")) %>%
    group_by(feature, Project, TimePoint) %>%
    summarise(mean_rel = mean(RelAbundance, na.rm = TRUE), .groups = "drop_last") %>%
    pivot_wider(names_from = TimePoint, values_from = mean_rel) %>%
    mutate(T0 = ifelse(is.na(T0), 0, T0),
           T72 = ifelse(is.na(T72), 0, T72),
           dir_project = case_when(
             T72 > T0 ~ "up",
             T72 < T0 ~ "down",
             TRUE ~ "tie"
           )) %>%
    ungroup()
  
  per_feature_consistency <- per_project %>%
    group_by(feature) %>%
    summarise(
      n_projects_used = n(),
      n_up   = sum(dir_project == "up", na.rm = TRUE),
      n_down = sum(dir_project == "down", na.rm = TRUE),
      prop_up = ifelse(n_projects_used > 0, n_up / n_projects_used, NA_real_),
      consistency_prop = case_when(
        is.na(prop_up) ~ NA_real_,
        prop_up >= 0.5 ~ prop_up,
        TRUE ~ 1 - prop_up
      ),
      majority_dir = case_when(
        is.na(prop_up) ~ NA_character_,
        prop_up > 0.5  ~ "up",
        prop_up < 0.5  ~ "down",
        TRUE ~ "tie"
      ),
      .groups = "drop"
    )
  
  # ---------- 4) add mean relabundance ----------
  avg_rel <- sub_df %>%
    group_by(feature = !!sym(feature_col)) %>%
    summarise(RelAbundance_mean = mean(!!sym(relab_col), na.rm = TRUE), .groups = "drop")
  
  # ---------- 5) Merge ----------
  merged <- res_main %>%
    left_join(real_fc_overall, by = "feature") %>%
    left_join(per_feature_consistency, by = "feature") %>%
    left_join(avg_rel, by = "feature") %>%
    arrange(FDR_lmm, feature)
  
  # ---------- 6) Select Remarkable ----------
  merged <- merged %>%
    mutate(
      Remarkable = case_when(
        (abs(log2FC72)   >= fc72_cutoff) &
          (abs(log2FC_real) >= fcreal_cutoff) &
          (!ci95_cross_0) &
          (consistency_prop > consistency_cutoff) &
          (nz_ratio >= nz_ratio_cutoff) &
          (FDR_lmm < 0.05) ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  list(merged = merged, per_project = per_project)
}
