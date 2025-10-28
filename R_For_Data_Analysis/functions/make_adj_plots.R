#-------------------------------------------------------------------------------
# make adjusted plots
#-------------------------------------------------------------------------------
make_adj_plots <- function(ps.obj, dists,
                           color_vec, shape_vec,
                           method = "PCoA",
                           ellipse = FALSE, hull = TRUE) {
  
  plots <- list()
  meta <- as(sample_data(ps.obj), "data.frame")
  
  for (dist in dists) {
    message(">>> Processing distance: ", dist)
    
    # Ordination
    ord <- ordinate(ps.obj, method = method, distance = dist)
    coords <- as.data.frame(ord$vectors[, 1:2])
    colnames(coords) <- c("Axis1", "Axis2")
    coords$Sample <- rownames(coords)
    
    # Safer merge with metadata
    meta.sub <- meta[match(coords$Sample, rownames(meta)), ]
    coords <- cbind(coords, meta.sub)
    
    # Debug checks
    message("Number of samples in ord: ", nrow(coords))
    message("Number of samples in meta: ", nrow(meta))
    message("Samples not matched: ", length(which(is.na(coords$Treat_ID))))
    if (any(is.na(coords$Treat_ID))) {
      message("Unmatched sample IDs: ", paste(coords$Sample[is.na(coords$Treat_ID)], collapse = ", "))
    }
    
    # Base plot
    p <- ggplot(coords, aes(x = Axis1, y = Axis2,
                            color = Treat_ID, shape = TimePoint,
                            group = Batch_ID)) +
      geom_path(alpha = 0.5, linewidth = 0.5) +   # connect by Batch_ID
      geom_point(alpha = 0.7, size = 2,
                 position = position_jitter(width = 0.01, height = 0.01)) +
      coord_fixed() +
      theme_bw()
    
    # Ellipse (per Treat_ID)
    if (ellipse) {
      p <- p + stat_ellipse(aes(fill = Treat_ID, group = Treat_ID), 
                            geom = "polygon",
                            alpha = 0.2, level = 0.90, color = NA)
    }
    
    # Hull (per Treat_ID)
    if (hull) {
      hull_df <- coords %>%
        group_by(Treat_ID) %>%
        slice(chull(Axis1, Axis2))
      p <- p + geom_polygon(data = hull_df,
                            aes(x = Axis1, y = Axis2, fill = Treat_ID, group = Treat_ID),
                            alpha = 0.2, color = NA)
    }
    
    # Colors and shapes
    if (!is.null(color_vec)) {
      p <- p + scale_color_manual(values = color_vec) +
        scale_fill_manual(values = color_vec)
    }
    if (!is.null(shape_vec)) {
      p <- p + scale_shape_manual(values = shape_vec)
    }
    
    plots[[dist]] <- p
  }
  
  return(plots)
}
