plot_cluster_sankey <- function(cluster_long,
                                time_col = "TimePoint",
                                cluster_col = "Cluster",
                                batch_col = "Batch_ID",
                                fill_by = "Cluster",
                                title = "Cluster transitions across time") {
  library(ggplot2)
  library(ggalluvial)
  

  stopifnot(all(c(time_col, cluster_col, batch_col) %in% colnames(cluster_long)))
  
  raw_levels <- unique(as.character(cluster_long[[cluster_col]]))
  is_num <- grepl("^[0-9]+$", raw_levels)
  num_lvls <- sort(as.numeric(raw_levels[is_num]))
  non_num_lvls <- raw_levels[!is_num]
  desired_levels <- c(as.character(num_lvls), non_num_lvls)
  

  cluster_long[[cluster_col]] <- factor(cluster_long[[cluster_col]], levels = desired_levels)

  if (fill_by != cluster_col) {
    cluster_long[[fill_by]] <- factor(cluster_long[[fill_by]],
                                      levels = desired_levels)
  }
  

  ggplot(cluster_long,
         aes(x = .data[[time_col]],
             stratum = .data[[cluster_col]],
             alluvium = .data[[batch_col]],
             fill = .data[[fill_by]])) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.6) +  
    geom_stratum(width = 0.25, color = "grey30") +                          
    scale_fill_manual(
      values = cluster_colors,
      breaks = desired_levels,
      drop = FALSE
    ) +
    theme_minimal(base_size = 14) +
    labs(
      x = "Time Point",
      y = "Number of Batch_IDs",
      title = title,
      subtitle = "Flow width ∝ number of Batch_IDs",
      fill = fill_by
    ) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}
