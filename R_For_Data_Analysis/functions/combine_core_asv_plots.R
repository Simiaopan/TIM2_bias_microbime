combine_core_asv_plots <- function(p90, p80, p50, p25,
                                   xlab = "TimePoint",
                                   ylab = "Sum of core-ASV relative abundance per sample",
                                   vspace = 0.12) {
  library(ggplot2)
  library(patchwork)
  
  
  strip_axes_titles <- function(p) {
    p +
      labs(title = NULL, x = NULL, y = NULL) +
      theme(
        legend.position = "none",
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.3)
      )
  }
  
  p1 <- strip_axes_titles(p90) + ggtitle("Core_ASV : 90%")
  p2 <- strip_axes_titles(p80) + ggtitle("Core_ASV : 80%")
  p3 <- strip_axes_titles(p50) + ggtitle("Core_ASV : 50%")
  p4 <- strip_axes_titles(p25) + ggtitle("Core_ASV : 25%")
  
  title_theme <- theme(
    plot.title = element_text(hjust = 0.5, vjust = -2, face = "bold", size = 11)
  )
  
  p1 <- p1 + title_theme
  p2 <- p2 + title_theme
  p3 <- p3 + title_theme
  p4 <- p4 + title_theme
  
  top_row <- p1 + p2 + plot_layout(ncol = 2, widths = c(1, 1))
  bottom_row <- p3 + p4 + plot_layout(ncol = 2, widths = c(1, 1))
  
  final_plot <- top_row / bottom_row + plot_layout(heights = c(1, 1), guides = "collect") & 
    theme(plot.margin = margin(5, 5, 5, 5))
  
  final_plot <- final_plot &
    plot_annotation(
      theme = theme(
        plot.margin = margin(20, 20, 20, 20)
      ),
      caption = NULL,
      tag_levels = NULL
    ) &
    labs(x = xlab, y = ylab)
  
  final_plot <- final_plot + plot_layout(heights = c(1, 1 + vspace))
  
  return(final_plot)
}
