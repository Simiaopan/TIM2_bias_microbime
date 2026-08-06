make_netcomi_igraph_one <- function(net_res, subset, candidate_features = NULL) {
  
  if (!requireNamespace("igraph", quietly = TRUE)) stop("Need igraph.")
  
  edge_df <- net_res$constructed[[subset]]$edgelist1
  
  cent <- net_res$analyzed[[subset]]$centralities
  clus <- net_res$analyzed[[subset]]$clustering$clust1
  hubs <- net_res$analyzed[[subset]]$hubs$hubs1
  
  nodes <- data.frame(
    name = names(cent$degree1),
    degree = as.numeric(cent$degree1),
    betweenness = as.numeric(cent$between1),
    closeness = as.numeric(cent$close1),
    eigenvector = as.numeric(cent$eigenv1),
    module = as.factor(clus[names(cent$degree1)]),
    hub = names(cent$degree1) %in% hubs,
    stringsAsFactors = FALSE
  )
  
  if (is.list(candidate_features)) {
    cand <- candidate_features[[subset]]
  } else {
    cand <- candidate_features
  }
  
  nodes$candidate <- nodes$name %in% cand
  
  g <- igraph::graph_from_data_frame(
    d = edge_df,
    vertices = nodes,
    directed = FALSE
  )
  
  igraph::E(g)$sign <- ifelse(igraph::E(g)$asso > 0, "positive", "negative")
  igraph::E(g)$abs_asso <- abs(igraph::E(g)$asso)
  
  g
}
