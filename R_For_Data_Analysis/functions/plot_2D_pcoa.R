#-------------------------------------------------------------------------------
# make 2D plot for pcoa: will generate all dists plot
#-------------------------------------------------------------------------------
plot_2D_pcoa <- function(pcoa_results, # must be a list
                         color_col = "Treat_ID",
                         shape_col = "TimePoint",
                         path_group = "Batch_ID",
                         color_vec = NULL,
                         shape_vec = NULL,
                         add_path = TRUE,
                         add_ellipse = TRUE) {
  
  # make sure use the correct object
  color_col <- as.character(color_col)[1]
  shape_col <- as.character(shape_col)[1]
  path_group <- as.character(path_group)[1]
  
  plots <- list()
  
  for (dist_name in names(pcoa_results)) {
    res <- pcoa_results[[dist_name]]
    df   <- res$coords_2
    expl <- res$expl_2
    
    xlab_txt <- if (is.na(expl[1])) "Axis 1" else sprintf("Axis 1  [%.1f%%]", expl[1])
    ylab_txt <- if (is.na(expl[2])) "Axis 2" else sprintf("Axis 2  [%.1f%%]", expl[2])
    
    p <- ggplot(df, aes(x = Axis1, y = Axis2))
    
    # options: connection lines
    if (add_path) {
      p <- p + geom_path(aes(group = .data[[path_group]],
                             color = .data[[color_col]]),
                         linewidth = 0.5, alpha = 0.6)
    }
    
    # options: ellipse
    if (add_ellipse) {
      p <- p + stat_ellipse(aes(group = .data[[color_col]],
                                fill  = .data[[color_col]]),
                            geom = "polygon", alpha = 0.15,
                            level = 0.90, type = "t", linewidth = 0.4)
    }
    
    p <- p + geom_point(aes(color = .data[[color_col]],
                            shape = .data[[shape_col]]),
                        size = 2.2, alpha = 0.95) +
      coord_equal() +
      theme_bw() +
      labs(title = paste0("PCoA - ", dist_name),
           x = xlab_txt, y = ylab_txt) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(size = 11, face = "bold"))
    
    if (!is.null(color_vec)) {
      p <- p + scale_color_manual(values = color_vec) +
        scale_fill_manual(values  = color_vec)
    }
    if (!is.null(shape_vec)) {
      p <- p + scale_shape_manual(values = shape_vec)
    }
    
    plots[[dist_name]] <- p
  }
  
  return(plots)
}
