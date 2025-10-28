#-------------------------------------------------------------------------------
# View the denoising report to select para for filtering
#-------------------------------------------------------------------------------

view_report <- function(file, plot = TRUE, summary = TRUE, log_scale = TRUE) {
  # Step 1: Read the denoising report
  stats <- read.delim(file, comment.char = "#", check.names = FALSE)
  rownames(stats) <- stats$`sample-id`
  stats$Sample <- stats$`sample-id`
  
  # Step 2: Define the steps of interest
  step_cols <- c("input", "filtered", "denoised", "merged", "non-chimeric")
  
  # Step 3: Check if all required columns exist
  missing_cols <- setdiff(step_cols, colnames(stats))
  if (length(missing_cols) > 0) {
    stop(paste("Missing expected columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Step 4: Convert to long format using pivot_longer
  df_long <- stats %>%
    select(Sample, all_of(step_cols)) %>%
    pivot_longer(cols = -Sample, names_to = "Step", values_to = "ReadCount")
  
  # Step 5: Plot boxplot (if requested)
  if (plot) {
    p <- ggplot(df_long, aes(x = Step, y = ReadCount)) +
      geom_boxplot(fill = "lightblue") +
      theme_minimal() +
      labs(title = "Read count per processing step", y = "Read count", x = "Step")
    
    if (log_scale) {
      p <- p + scale_y_log10() + labs(y = "Read count (log10)")
    }
    
    print(p)
  }
  
  # Step 6: Calculate summary statistics (if requested)
  if (summary) {
    step_counts <- stats[, step_cols]
    summary_stats <- data.frame(
      Step = step_cols,
      Mean = apply(step_counts, 2, mean),
      Median = apply(step_counts, 2, median),
      Min = apply(step_counts, 2, min),
      Max = apply(step_counts, 2, max)
    )
    # Round only numeric columns
    num_cols <- sapply(summary_stats, is.numeric)
    summary_stats[ , num_cols] <- round(summary_stats[ , num_cols], digits = 2)
    
    return(summary_stats)
  } else {
    return(invisible(NULL))
  }
}
