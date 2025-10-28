#-------------------------------------------------------------------------------
# Prepare MaAsLin2 input
#-------------------------------------------------------------------------------
prep_maaslin_inputs <- function(df, feature_col) {
  df_sum <- df %>%
    group_by(Sample, !!sym(feature_col)) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  feat <- df_sum %>%
    tidyr::pivot_wider(id_cols = Sample,
                       names_from = !!sym(feature_col),
                       values_from = Abundance,
                       values_fill = 0) %>%
    tibble::column_to_rownames("Sample")
  
  meta <- df %>%
    distinct(Sample, Time, Project, Batch_ID, Treat_ID) %>%
    mutate(Time_c = as.numeric(scale(Time, center = TRUE, scale = FALSE))) %>%
    tibble::column_to_rownames("Sample")
  
  list(feat = feat, meta = meta)
}
