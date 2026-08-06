extract_netcomi_network_table <- function(g_list) {
  
  res <- lapply(names(g_list), function(subset) {
    
    g <- g_list[[subset]]
    
    edge_sign <- igraph::E(g)$sign
    
    data.frame(
      n_nodes = igraph::vcount(g),
      n_edges = igraph::ecount(g),
      density = igraph::edge_density(g),
      average_degree = mean(igraph::V(g)$degree),
      n_modules = length(unique(igraph::V(g)$module)),
      n_hubs = sum(igraph::V(g)$hub),
      n_candidates = sum(igraph::V(g)$candidate),
      positive_edges = sum(edge_sign == "positive"),
      negative_edges = sum(edge_sign == "negative")
    )
  })
  
  names(res) <- names(g_list)
  
  res
}