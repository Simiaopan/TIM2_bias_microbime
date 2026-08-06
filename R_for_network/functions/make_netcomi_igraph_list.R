make_netcomi_igraph_list <- function(net_res,
                                     subsets = c("all", "t1", "t1ts"),
                                     candidate_features = NULL) {
  
  g_list <- lapply(subsets, function(subset) {
    make_netcomi_igraph_one(
      net_res = net_res,
      subset = subset,
      candidate_features = candidate_features
    )
  })
  
  names(g_list) <- subsets
  
  g_list
}