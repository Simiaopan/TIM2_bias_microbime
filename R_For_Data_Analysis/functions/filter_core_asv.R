#-------------------------------------------------------------------------------
# Filter core ASVs based on prevalence, relative abundance, and project occurrence at T0
#-------------------------------------------------------------------------------
filter_core_asv <- function(asv_raw,
                            prevalence_cutoff = 0.8,
                            relabundance_cutoff = 0,
                            min_project_count = 3) {
  
  library(dplyr)
  
  # Step 1: Subset only samples at TimePoint == "T0"
  asv_T0 <- asv_raw %>%
    filter(TimePoint == "T0")
  
  # Step 2: Calculate prevalence, mean relative abundance, and project occurrence per OTU
  prevalence_df <- asv_T0 %>%
    group_by(OTU) %>%
    summarise(
      prevalence = sum(Abundance > 0) / n_distinct(Sample),
      mean_rel_abundance = mean(RelAbundance, na.rm = TRUE),
      n_projects = n_distinct(Project),       # number of unique projects containing this OTU
      .groups = "drop"
    )
  
  # Step 3: Identify core ASVs based on thresholds
  core_otus <- prevalence_df %>%
    filter(prevalence >= prevalence_cutoff,
           mean_rel_abundance >= relabundance_cutoff,
           n_projects >= min_project_count) %>%   # must appear in ≥ X projects
    pull(OTU)
  
  # Step 4: Filter the original dataset to retain only the core ASVs
  asv_filtered <- asv_raw %>%
    filter(OTU %in% core_otus)
  
  # Step 5: Print summary information
  total_asv <- length(unique(asv_raw$OTU))
  kept_asv <- length(core_otus)
  cat(
    sprintf(
      "Retained %d core ASVs (%.2f%% of total).\n",
      kept_asv,
      kept_asv / total_asv * 100
    )
  )
  cat(sprintf("   Each retained ASV appeared in ≥ %d projects at T0.\n", min_project_count))
  
  # Step 6: Return results
  return(list(
    filtered_table = asv_filtered,
    keep_asv = core_otus
  ))
}
