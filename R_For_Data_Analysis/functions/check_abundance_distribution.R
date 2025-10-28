#-------------------------------------------------------------------------------
# check abundance and distribution
#-------------------------------------------------------------------------------
check_abundance_distribution <- function(df, level = "Genus", value_col = "Abundance") {
  library(ggplot2)
  library(dplyr)
  
  qc_tbl <- df %>%
    group_by(.data[[level]]) %>%
    summarise(
      Prevalence = mean(.data[[value_col]] > 0),
      MeanAbundance = mean(.data[[value_col]], na.rm = TRUE),
      ZeroProp = mean(.data[[value_col]] == 0),
      .groups = "drop"
    )
  
  p1 <- ggplot(qc_tbl, aes(x = MeanAbundance)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "white") +
    scale_x_log10() + theme_minimal() +
    labs(title = paste("Distribution of Mean Abundance (", level, ")", sep=""))
  
  p2 <- ggplot(qc_tbl, aes(x = Prevalence)) +
    geom_histogram(bins = 50, fill = "salmon", color = "white") +
    theme_minimal() +
    labs(title = paste("Prevalence Distribution (", level, ")", sep=""))
  
  p3 <- ggplot(qc_tbl, aes(x = ZeroProp)) +
    geom_histogram(bins = 50, fill = "lightgreen", color = "white") +
    theme_minimal() +
    labs(title = paste("Zero Proportion Distribution (", level, ")", sep=""))
  
  return(list(summary = qc_tbl, plots = list(p1, p2, p3)))
}
