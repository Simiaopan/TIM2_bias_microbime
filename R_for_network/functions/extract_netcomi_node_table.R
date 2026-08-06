extract_netcomi_node_table <- function(g_list) {
  
  res <- lapply(names(g_list), function(subset) {
    
    g <- g_list[[subset]]
    
    data.frame(
      ASV = igraph::V(g)$name,
      degree = igraph::V(g)$degree,
      betweenness = igraph::V(g)$betweenness,
      closeness = igraph::V(g)$closeness,
      eigenvector = igraph::V(g)$eigenvector,
      module = igraph::V(g)$module,
      hub = igraph::V(g)$hub,
      candidate = igraph::V(g)$candidate,
      stringsAsFactors = FALSE
    )
  })
  
  names(res) <- names(g_list)
  
  res
}