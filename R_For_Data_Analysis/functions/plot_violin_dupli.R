#-------------------------------------------------------------------------------
# Violin plot (Within vs Between) with optional trajectories and significance marks
# Significance marks are shown only above the Between group at each timepoint
#-------------------------------------------------------------------------------
plot_violin_dupli <- function(dupli_res, aes.ls,
                              mode = c("point", "line"),
                              proj_col_map = "col_alt",
                              time_levels = NULL,
                              show_box = TRUE) {
  library(ggplot2)
  library(dplyr)
  
  mode <- match.arg(mode)
  
  #----------------------------
  # Extract data
  #----------------------------
  df <- dupli_res$distances
  stats <- dupli_res$stats
  
  # Ensure time order
  if (!is.null(time_levels)) {
    df$TimePoint <- factor(df$TimePoint, levels = time_levels, ordered = TRUE)
    stats$TimePoint <- factor(stats$TimePoint, levels = time_levels, ordered = TRUE)
  }
  
  # Color palette
  proj_cols <- aes.ls[[proj_col_map]][["Project"]]
  
  #----------------------------
  # Base violin + boxplot
  #----------------------------
  p <- ggplot(df, aes(x = interaction(TimePoint, Group), y = Distance, fill = Group)) +
    geom_violin(alpha = 0.3, scale = "width", trim = FALSE)
  
  if (show_box) {
    p <- p + geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.4,
                          position = position_dodge(width = 0.9))
  }
  
  #----------------------------
  # Add points or lines
  #----------------------------
  if (mode == "point") {
    df$Color <- ifelse(df$Group == "Within", as.character(df$Project), "Between")
    color_map <- c(proj_cols, "Between" = "grey70")
    
    df_within  <- df[df$Group == "Within", ]
    df_between <- df[df$Group == "Between", ]
    
    p <- p +
      geom_jitter(data = df_between,
                  aes(color = Color),
                  width = 0.15, alpha = 0.3, size = 0.6) +
      geom_jitter(data = df_within,
                  aes(color = Color),
                  width = 0.15, alpha = 0.9, size = 1.8) +
      scale_color_manual(values = color_map, na.value = "grey70")
    
  } else { # mode == "line"
    df_within  <- df[df$Group == "Within", ]
    df_between <- df[df$Group == "Between", ]
    
    p <- p +
      geom_line(data = df_within,
                aes(x = interaction(TimePoint, Group),
                    y = Distance,
                    group = PairID,
                    color = Project),
                linewidth = 0.6, alpha = 0.7) +
      geom_line(data = df_between,
                aes(x = interaction(TimePoint, Group),
                    y = Distance,
                    group = PairID),
                color = "grey70", linewidth = 0.4, alpha = 0.5) +
      scale_color_manual(values = proj_cols, na.value = "grey70")
  }
  
  #----------------------------
  # Add significance marks
  #----------------------------
  # Compute position of y (slightly above Between max)
  y_pos_df <- df %>%
    filter(Group == "Between") %>%
    group_by(TimePoint) %>%
    summarise(y_pos = max(Distance, na.rm = TRUE) * 1.05)
  
  stats <- left_join(stats, y_pos_df, by = "TimePoint") %>%
    mutate(
      x_pos = paste0(TimePoint, ".Between")  # show above Between group only
    )
  
  p <- p +
    geom_text(
      data = stats,
      inherit.aes = FALSE,
      aes(x = x_pos, y = y_pos, label = signif),
      vjust = 0,
      size = 4.5,
      fontface = "bold",
      color = "black"
    )
  
  #----------------------------
  # Final theme and labels
  #----------------------------
  p <- p +
    scale_fill_manual(values = c("Within" = "#1f78b4", "Between" = "#a6cee3")) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = "TimePoint × Group", y = "Distance", color = "Project", fill = "Group")
  
  return(p)
}
