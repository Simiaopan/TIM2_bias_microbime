plot_netcomi_node_metric_boxplot <- function(node_table_list,
                                             stats_table_list,
                                             metrics = c("degree", "betweenness", "closeness", "eigenvector"),
                                             group_col = "candidate",
                                             group_labels = c("Background", "Candidate"),
                                             fill_cols = c("Background" = "#a6cee3", "Candidate" = "#1f78b4"), 
                                             x_pos = c(1, 1.45),        
                                             box_width = 0.20,         
                                             mean_size = 2.5,
                                             mean_alpha = 0,
                                             sign_size = 3.5,          
                                             add_points = TRUE,        
                                             point_size = 0.8,         
                                             point_alpha = 0.4,        
                                             jitter_width = 0.12,      
                                             x_expand = 0.35,          
                                             aspect_ratio = 1.2,  
                                             bracket_gap = 0.08,       
                                             bracket_height = 0.04,    
                                             text_gap = 0.04) {        
  
  library(ggplot2)
  library(dplyr)
  
  plot_list <- list()
  
  for (subset in names(node_table_list)) {
    
    df0 <- node_table_list[[subset]]
    stats0 <- stats_table_list[[subset]]
    
    plot_list[[subset]] <- list()
    
    for (met in metrics) {
      
      if (!met %in% colnames(df0)) {
        warning("Metric not found in subset ", subset, ": ", met)
        next
      }
      
      df <- df0 %>%
        dplyr::select(
          group_raw = dplyr::all_of(group_col),
          value = dplyr::all_of(met)
        ) %>%
        dplyr::filter(!is.na(group_raw), !is.na(value)) %>%
        dplyr::mutate(
          group = ifelse(group_raw, group_labels[2], group_labels[1]),
          group = factor(group, levels = group_labels),
          xpos = ifelse(group == group_labels[1], x_pos[1], x_pos[2])
        )
      
      stats <- stats0 %>%
        dplyr::filter(.data$metric == met)
      
      signif_label <- if (nrow(stats) > 0) stats$signif[1] else "ns"
      
      y_min_real <- min(df$value, na.rm = TRUE)
      y_max_real <- max(df$value, na.rm = TRUE)
      y_range_real <- y_max_real - y_min_real
      if (!is.finite(y_range_real) || y_range_real == 0) y_range_real <- 1
      
      line_y_low <- y_max_real + bracket_gap * y_range_real
      line_y_top <- line_y_low + bracket_height * y_range_real
      text_y     <- line_y_top + text_gap * y_range_real
      y_upper    <- text_y + 0.10 * y_range_real 
      
      anno_bracket <- data.frame(
        x = c(x_pos[1], x_pos[1], x_pos[2], x_pos[2]),
        y = c(line_y_low, line_y_top, line_y_top, line_y_low)
      )
      
      y_lower_limit <- y_min_real - 0.05 * y_range_real
      if(y_min_real == 0) {
        y_lower_limit <- -0.05 * y_range_real
      }
      
      p <- ggplot(df, aes(x = xpos, y = value, fill = group, group = group))
      
      if (add_points) {
        p <- p + geom_jitter(
          aes(color = group),
          position = position_jitterdodge(jitter.width = jitter_width, dodge.width = 0),
          size = point_size,
          alpha = point_alpha,
          show.legend = FALSE
        )
      }
      
      p <- p + geom_boxplot(
        width = box_width,
        alpha = 0.6,                 
        color = "black",             
        linewidth = 0.5,
        outlier.shape = if(add_points) NA else 16,
        outlier.color = "grey50",
        outlier.size = 1.0,
        fatten = 1.5                 
      ) +
        stat_summary(
          fun = mean,
          geom = "point",
          color = "black",
          size = mean_size,
          alpha = mean_alpha
        ) +
        geom_path(
          data = anno_bracket,
          aes(x = x, y = y),
          inherit.aes = FALSE,
          linewidth = 0.45,
          color = "black"
        ) +
        geom_text(
          data = data.frame(
            x = mean(x_pos),
            y = text_y,
            label = signif_label
          ),
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE,
          vjust = 0,
          size = sign_size,
          fontface = "bold"
        ) +
        scale_fill_manual(values = fill_cols) +
        scale_color_manual(values = fill_cols) + 
        scale_x_continuous(
          breaks = x_pos,
          labels = group_labels,
          limits = c(x_pos[1] - x_expand, x_pos[2] + x_expand),
          expand = expansion(mult = c(0.02, 0.02))
        ) +
        coord_cartesian(
          ylim = c(y_lower_limit, y_upper)
        ) +
        theme_classic(base_size = 12) +
        theme(
          axis.line = element_line(color = "black"),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"),
          aspect.ratio = aspect_ratio  
        ) +
        labs(
          title = paste0(subset, " - ", met),
          x = NULL,
          y = met
        )
      
      plot_list[[subset]][[met]] <- p
    }
  }
  
  return(plot_list)
}