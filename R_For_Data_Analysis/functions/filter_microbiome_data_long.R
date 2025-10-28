#-------------------------------------------------------------------------------
# filtering for long-format absolute abundance data
#-------------------------------------------------------------------------------
filter_microbiome_data_long <- function(df,
                                        level = "OTU",
                                        abundance_col = "Abundance",
                                        sample_col = "Sample",
                                        prev_cutoff = 0.2,
                                        mean_cutoff = 1e-4) {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(glue))
  
  #--- cal relative abundance ---
  df <- df %>%
    group_by(.data[[sample_col]]) %>%
    mutate(RelAbundance = .data[[abundance_col]] / sum(.data[[abundance_col]], na.rm = TRUE)) %>%
    ungroup()
  
  #--- cal QC parameters ---
  summary_df <- df %>%
    group_by(.data[[level]]) %>%
    summarise(
      Prevalence = mean(.data[[abundance_col]] > 0, na.rm = TRUE),
      MeanRelAbundance = mean(RelAbundance, na.rm = TRUE),
      ZeroProp = mean(.data[[abundance_col]] == 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  #--- filter ---
  keep_ids <- summary_df %>%
    filter(Prevalence >= prev_cutoff, MeanRelAbundance >= mean_cutoff) %>%
    pull(.data[[level]])
  
  filtered_df <- df %>%
    filter(.data[[level]] %in% keep_ids)
  
  #--- results ---
  n_total <- nrow(summary_df)
  n_keep <- length(keep_ids)
  
  message(glue("Retained {n_keep}/{n_total} taxa ({round(n_keep/n_total*100,1)}%) ",
               "with Prevalence ≥ {prev_cutoff}, MeanRelAbundance ≥ {mean_cutoff}."))
  
  return(list(
    filtered_df = filtered_df,
    summary_df = summary_df,
    kept_ids = keep_ids
  ))
}
