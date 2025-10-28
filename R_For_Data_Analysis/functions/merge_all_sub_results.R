#-------------------------------------------------------------------------------
# merge sub and all
#-------------------------------------------------------------------------------
merge_all_sub_results <- function(res_all, res_sub, feature_col = "feature") {
  suppressPackageStartupMessages(library(dplyr))
  

  all_core <- res_all %>%
    select(feature = !!sym(feature_col),
           Remarkable_all = Remarkable,
           log2FC72_all = log2FC72)
  
  sub_core <- res_sub %>%
    select(feature = !!sym(feature_col),
           Remarkable_sub = Remarkable,
           log2FC72_sub = log2FC72)
  

  merged <- full_join(all_core, sub_core, by = "feature") %>%
    mutate(
      Both_Remarkable = case_when(
        Remarkable_all & Remarkable_sub &
          sign(log2FC72_all) == sign(log2FC72_sub) ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  res_all_with_flag <- res_all %>%
    left_join(merged %>% select(feature, Both_Remarkable),
              by = setNames("feature", feature_col))
  
  res_sub_with_flag <- res_sub %>%
    left_join(merged %>% select(feature, Both_Remarkable),
              by = setNames("feature", feature_col))
  

  n_both <- sum(merged$Both_Remarkable, na.rm = TRUE)
  message(n_both, " features are Both-Remarkable across all and sub datasets.")
  

  if (n_both > 0) {
    dir_summary <- merged %>%
      filter(Both_Remarkable) %>%
      summarize(
        up = sum(log2FC72_all > 0, na.rm = TRUE),
        down = sum(log2FC72_all < 0, na.rm = TRUE)
      )
    message("Among them: ", dir_summary$up, " upregulated, ",
            dir_summary$down, " downregulated.")
  }
  
  list(all_with_both = res_all_with_flag,
       sub_with_both = res_sub_with_flag)
}
