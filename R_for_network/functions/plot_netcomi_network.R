plot_netcomi_network <- function(layout_df,
                                 edge_mode = c("all",
                                               "candidate_connected",
                                               "candidate_candidate"),
                                 label_mode = c("none",
                                                "all",
                                                "candidate",
                                                "hub",
                                                "candidate_hub"),
                                 drop_isolated = FALSE,
                                 edge_alpha = 0.3,
                                 background_node_alpha = 0.35,
                                 candidate_node_alpha = 1,
                                 node_size_range = c(1, 6),
                                 edge_width = 0.25,
                                 label_size = 2,
                                 label_segment_size = 0.3,
                                 label_segment_color = "grey40",
                                 label_box_padding = 0.6,
                                 label_point_padding = 0.05,
                                 label_segment_alpha = 1,
                                 min_segment_length = 0) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Need package 'igraph'.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need package 'ggplot2'.")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Need package 'ggrepel'.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need package 'dplyr'.")
  
  edge_mode <- match.arg(edge_mode)
  label_mode <- match.arg(label_mode)
  
  g <- attr(layout_df, "graph")
  if (is.null(g)) {
    stop("layout_df does not contain graph attribute. Please create it with ggraph::create_layout().")
  }
  
  node_df <- as.data.frame(layout_df) |>
    dplyr::select(
      name, x, y,
      dplyr::everything()
    ) |>
    dplyr::mutate(
      candidate = as.logical(candidate),
      hub = as.logical(hub),
      candidate_str = ifelse(candidate, "Candidate", "Background")
    )
  
  edge_df <- igraph::as_data_frame(g, what = "edges")
  
  if (!"abs_asso" %in% colnames(edge_df)) {
    if (!"asso" %in% colnames(edge_df)) {
      stop("Edge table lacks both 'abs_asso' and 'asso'.")
    }
    edge_df$abs_asso <- abs(edge_df$asso)
  }
  
  if (!"sign" %in% colnames(edge_df)) {
    if (!"asso" %in% colnames(edge_df)) {
      stop("Edge table lacks both 'sign' and 'asso'.")
    }
    edge_df$sign <- ifelse(edge_df$asso > 0, "positive", "negative")
  }
  
  edge_df$sign <- as.character(edge_df$sign)
  
  candidate_nodes <- node_df$name[node_df$candidate]
  
  if (edge_mode == "all") {
    edge_keep <- rep(TRUE, nrow(edge_df))
  } else if (edge_mode == "candidate_connected") {
    edge_keep <- edge_df$from %in% candidate_nodes |
      edge_df$to %in% candidate_nodes
  } else if (edge_mode == "candidate_candidate") {
    edge_keep <- edge_df$from %in% candidate_nodes &
      edge_df$to %in% candidate_nodes
  }
  
  edge_df2 <- edge_df[edge_keep, , drop = FALSE]
  
  coord_from <- node_df |>
    dplyr::select(from = name, x, y)
  
  coord_to <- node_df |>
    dplyr::select(to = name, xend = x, yend = y)
  
  edge_plot_df <- edge_df2 |>
    dplyr::left_join(coord_from, by = "from") |>
    dplyr::left_join(coord_to, by = "to")
  
  if (any(is.na(edge_plot_df$x)) ||
      any(is.na(edge_plot_df$y)) ||
      any(is.na(edge_plot_df$xend)) ||
      any(is.na(edge_plot_df$yend))) {
    stop("Some edge endpoints do not have coordinates in layout_df.")
  }
  
  if (drop_isolated) {
    connected_nodes <- unique(c(edge_plot_df$from, edge_plot_df$to))
    node_df <- node_df |>
      dplyr::filter(name %in% connected_nodes)
  }
  
  label_df <- node_df
  
  if (label_mode == "none") {
    label_df <- label_df[0, , drop = FALSE]
  } else if (label_mode == "candidate") {
    label_df <- label_df |>
      dplyr::filter(candidate)
  } else if (label_mode == "hub") {
    label_df <- label_df |>
      dplyr::filter(hub)
  } else if (label_mode == "candidate_hub") {
    label_df <- label_df |>
      dplyr::filter(candidate & hub)
  }
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edge_plot_df,
      ggplot2::aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        colour = sign
      ),
      linewidth = edge_width,
      alpha = edge_alpha,
      show.legend = TRUE
    ) +
    ggplot2::scale_colour_manual(
      name = "Association",
      values = c(
        positive = "#D55E00",
        negative = "#0072B2"
      ),
      na.value = "grey70"
    ) +
    ggplot2::geom_point(
      data = node_df,
      ggplot2::aes(
        x = x,
        y = y,
        size = degree,
        fill = candidate_str,
        alpha = candidate_str
      ),
      shape = 21,
      colour = "black",
      stroke = 0.25
    ) +
    ggplot2::scale_fill_manual(
      name = "Node type",
      values = c(
        "Background" = "grey80",
        "Candidate"  = "#E64B35"
      )
    ) +
    ggplot2::scale_alpha_manual(
      name = "Node type",
      values = c(
        "Background" = background_node_alpha,
        "Candidate"  = candidate_node_alpha
      )
    ) +
    ggplot2::scale_size_continuous(
      name = "Degree",
      range = node_size_range
    ) +
    ggplot2::coord_equal() +
    ggplot2::theme_void()
  
  if (nrow(label_df) > 0) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        ggplot2::aes(
          x = x,
          y = y,
          label = name
        ),
        inherit.aes = FALSE,
        size = label_size,
        max.overlaps = Inf,
        box.padding = label_box_padding,
        point.padding = label_point_padding,
        segment.size = label_segment_size,
        segment.color = label_segment_color,
        segment.alpha = label_segment_alpha,
        min.segment.length = min_segment_length
      )
  }
  
  return(p)
}