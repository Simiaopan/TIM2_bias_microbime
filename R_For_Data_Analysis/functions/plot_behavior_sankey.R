#-------------------------------------------------------------------------------
# Plot Sankey diagram (T0 → T72) for selected Batch behaviors
#-------------------------------------------------------------------------------
plot_behavior_sankey_T0_T72 <- function(behavior_df,
                                        cluster_long,
                                        behaviors_to_include = c("Converging", "Dynamic transition"),
                                        time_col = "TimePoint",
                                        cluster_col = "Cluster",
                                        batch_col = "Batch_ID",
                                        fill_by = "Cluster",
                                        title = "Cluster transitions (T0 → T72) of selected behaviors") {
  library(dplyr)
  library(ggplot2)
  library(ggalluvial)
  
  # 1. Clean data
  behavior_df <- as.data.frame(behavior_df)
  cluster_long <- as.data.frame(cluster_long)
  
  # 2. Select Batch_IDs based on behavior
  selected_batches <- behavior_df %>%
    filter(behavior %in% behaviors_to_include) %>%
    pull(all_of(batch_col)) %>%
    unique()
  
  # 3. Subset long table to selected Batch_IDs and only T0/T72
  sub_long <- cluster_long %>%
    filter(.data[[batch_col]] %in% selected_batches) %>%
    filter(.data[[time_col]] %in% c("T0", "T72")) %>%
    mutate(
      !!time_col := factor(.data[[time_col]], levels = c("T0", "T72")),
      !!cluster_col := as.character(.data[[cluster_col]])
    )
  
  # 4. Check that both T0 and T72 exist for each Batch_ID
  sub_long <- sub_long %>%
    group_by(.data[[batch_col]]) %>%
    filter(length(unique(.data[[time_col]])) == 2) %>%
    ungroup()
  
  # 5. Plot Sankey (T0 → T72)
  p <- ggplot(sub_long,
              aes(x = .data[[time_col]],
                  stratum = .data[[cluster_col]],
                  alluvium = .data[[batch_col]],
                  fill = .data[[fill_by]])) +
    geom_flow(stat = "alluvium", alpha = 0.6) +
    geom_stratum(width = 0.25, color = "grey30") +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Time Point",
      y = "Number of Batch_IDs",
      subtitle = paste("Included behaviors:", paste(behaviors_to_include, collapse = ", ")),
      fill = "Cluster"
    ) +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  message("Sankey plot created (T0 → T72) for behaviors: ", paste(behaviors_to_include, collapse = ", "))
  return(p)
}
