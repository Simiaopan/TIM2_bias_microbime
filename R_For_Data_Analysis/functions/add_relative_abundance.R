#-------------------------------------------------------------------------------
# Add relative abundance column to long-format microbiome data
#-------------------------------------------------------------------------------
add_relative_abundance <- function(df,
                                   abundance_col = "Abundance",
                                   sample_col = "Sample") {
  suppressPackageStartupMessages(library(dplyr))
  
  df <- df %>%
    group_by(.data[[sample_col]]) %>%
    mutate(RelAbundance = .data[[abundance_col]] / sum(.data[[abundance_col]], na.rm = TRUE)) %>%
    ungroup()
  
  return(df)
}
