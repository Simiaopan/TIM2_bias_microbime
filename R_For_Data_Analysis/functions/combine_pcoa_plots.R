#-------------------------------------------------------------------------------
# combine plots from different dists
#-------------------------------------------------------------------------------
combine_pcoa_plots <- function(plots.list, ncol = 2) {
  require(cowplot)
  
  #remove all legends and only keep the one of first plot
  plots.nolegend <- lapply(plots.list, function(p) p + theme(legend.position = "none"))
  legend <- cowplot::get_legend(plots.list[[1]] + theme(legend.position = "right"))
  
  # combine
  p.grid <- cowplot::plot_grid(plotlist = plots.nolegend, ncol = ncol)
  combined <- cowplot::plot_grid(p.grid, legend, rel_widths = c(0.85, 0.15))
  
  return(combined)
}
