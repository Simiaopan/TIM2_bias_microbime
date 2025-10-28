#-------------------------------------------------------------------------------
# Comprehensive ASV filtering for longitudinal multi-project microbiome data
# (expects RelAbundance already computed)
#-------------------------------------------------------------------------------
filter_asv_time_project <- function(df,
                                       otu_col = "OTU",
                                       abundance_col = "Abundance",
                                       relab_col = "RelAbundance",
                                       time_col = "Time",
                                       timepoint_col = "TimePoint",
                                       project_col = "Project",
                                       prev_cutoff_time = 0.25,
                                       mean_cutoff = 1e-4,
                                       min_projects = 3,
                                       min_timepoints_per_project = 2) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(glue)
  })
  
  if (!relab_col %in% colnames(df)) {
    stop(glue("❌ The input data must contain a `{relab_col}` column. Please run `add_relative_abundance()` first."))
  }
  
  # average relative abundance
  message(glue("Step 1: Filter MeanRelAbundance ≥ {mean_cutoff}"))
  global_summary <- df %>%
    group_by(.data[[otu_col]]) %>%
    summarise(MeanRelAbundance = mean(.data[[relab_col]], na.rm = TRUE),
              .groups = "drop")
  
  keep_mean <- global_summary %>%
    filter(MeanRelAbundance >= mean_cutoff) %>%
    pull(.data[[otu_col]])
  df <- df %>% filter(.data[[otu_col]] %in% keep_mean)
  
  # filter by timepoint prevalence
  message(glue("Step 2: Filter by prevalence ≥ {prev_cutoff_time} in at least one time point ..."))
  time_prev <- df %>%
    group_by(.data[[otu_col]], .data[[timepoint_col]]) %>%
    summarise(Prev = mean(.data[[abundance_col]] > 0, na.rm = TRUE), .groups = "drop")
  
  keep_time <- time_prev %>%
    group_by(.data[[otu_col]]) %>%
    summarise(PassTime = any(Prev >= prev_cutoff_time), .groups = "drop") %>%
    filter(PassTime) %>%
    pull(.data[[otu_col]])
  
  df <- df %>% filter(.data[[otu_col]] %in% keep_time)
  
  # filter by project presenting
  message(glue("Step 3: Must appear in ≥ {min_projects} projects ..."))
  proj_cov <- df %>%
    group_by(.data[[otu_col]], .data[[project_col]]) %>%
    summarise(nonzero = any(.data[[abundance_col]] > 0), .groups = "drop") %>%
    group_by(.data[[otu_col]]) %>%
    summarise(n_projects = sum(nonzero), .groups = "drop")
  
  keep_proj <- proj_cov %>%
    filter(n_projects >= min_projects) %>%
    pull(.data[[otu_col]])
  
  df <- df %>% filter(.data[[otu_col]] %in% keep_proj)
  
  # filter by timepoint in project
  message(glue("Step 4: Each project must have ≥ {min_timepoints_per_project} nonzero time points ..."))
  proj_time_cov <- df %>%
    mutate(nonzero = .data[[abundance_col]] > 0) %>%
    group_by(.data[[otu_col]], .data[[project_col]], .data[[time_col]]) %>%
    summarise(nz = any(nonzero), .groups = "drop") %>%
    group_by(.data[[otu_col]], .data[[project_col]]) %>%
    summarise(n_time_nonzero = sum(nz), .groups = "drop") %>%
    group_by(.data[[otu_col]]) %>%
    summarise(n_projects_pass = sum(n_time_nonzero >= min_timepoints_per_project),
              .groups = "drop")
  
  keep_proj_time <- proj_time_cov %>%
    filter(n_projects_pass >= min_projects) %>%
    pull(.data[[otu_col]])
  
  df_filtered <- df %>% filter(.data[[otu_col]] %in% keep_proj_time)
  
  # Summary
  summary_df <- global_summary %>%
    left_join(proj_cov, by = otu_col) %>%
    left_join(proj_time_cov, by = otu_col)
  
  message("✅ Filtering completed.")
  message(glue("Retained {length(unique(df_filtered[[otu_col]]))} ASVs."))
  
  return(list(
    filtered_df = df_filtered,
    summary_df = summary_df,
    kept_OTUs = unique(df_filtered[[otu_col]])
  ))
}
