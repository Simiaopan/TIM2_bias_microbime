#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("../R_For_Data_Analysis/rdata/data.ls_for_LMM.Rdata")
load("../R_For_Data_Analysis/rdata/LMM_LinDA_MaAsLin2_match.Rdata")
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(ggraph)
  library(glue)
  library(NetCoMi)
})
 set.seed(20021)
#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/build_netcomi_networks.R")
source("functions/make_netcomi_igraph_one.R")
source("functions/make_netcomi_igraph_list.R")
source("functions/extract_netcomi_node_table.R")
source("functions/extract_netcomi_network_table.R")
source("functions/analyze_candidate_node_metrics.R")
source("functions/plot_netcomi_node_metric_boxplot.R")
source("functions/plot_netcomi_network.R")
 
#-------------------------------------------------------------------------------
# Run NetCoMi
#-------------------------------------------------------------------------------
net_res_0.3 <- build_netcomi_networks(
  ps_list = data.ls_for_LMM,
  candidate_features = LinDA_match$validated_vectors,
  n_background = 400,
  measure = "spearman",
  normMethod = "clr",
  zeroMethod = "pseudo",
  sparsMethod = "threshold",
  thresh = 0.3,
  graphlet = FALSE,
  gcmHeat = FALSE
)

# sapply(net_res$constructed, function(x) {
#   igraph::ecount(x$graph)
# })
# 
# sapply(net_res$constructed, function(x) {
#   igraph::vcount(x$graph)
# })
# 
# sapply(net_res$constructed, function(x) {
#   igraph::edge_density(x$graph)
# })

#-------------------------------------------------------------------------------
# Extract NetCoMi results for plot and analysis
#-------------------------------------------------------------------------------
netcomi_igraph_list_0.3 <- make_netcomi_igraph_list(
  net_res = net_res_0.3,
  subsets = c("all", "t1", "t1ts"),
  candidate_features = LinDA_match$validated_vectors
)

netcomi_node_table_0.3 <- extract_netcomi_node_table(
  netcomi_igraph_list_0.3
)

netcomi_network_table_0.3 <- extract_netcomi_network_table(
  netcomi_igraph_list_0.3
)

#-------------------------------------------------------------------------------
# Analyze candidate vs background nodes metrics
#-------------------------------------------------------------------------------
node_metric_tests_0.3 <- analyze_candidate_node_metrics(
  node_table_list = netcomi_node_table_0.3,
  metrics = c("degree", "betweenness", "closeness", "eigenvector"),
  group_col = "candidate",
  method = "wilcox",
  p_adjust_method = "BH"
)

# make box plots
netcomi_box_plots_0.3 <- plot_netcomi_node_metric_boxplot(
  node_table_list = netcomi_node_table_0.3,
  stats_table_list = node_metric_tests_0.3,
  metrics = c("degree", "betweenness", "closeness", "eigenvector"),
  mean_alpha = 0,
  box_width = 0.18,
  fill_cols = c("Background" = "grey80", "Candidate" = "#E64B35"),
  aspect_ratio = 1.1,
  point_size = 1, point_alpha = 0.2
)

netcomi_box_plots_0.3$all
#-------------------------------------------------------------------------------
# Make network plots
#-------------------------------------------------------------------------------
netcomi_layout_list_0.3 <- lapply(
  netcomi_igraph_list_0.3,
  function(g) {
    ggraph::create_layout(g, layout = "fr")
  }
)

netcomi_p_all_0.3_type1 <- plot_netcomi_network(
  layout_df = netcomi_layout_list_0.3$all,
  edge_alpha = 0.1,
  background_node_alpha = 0.15,
  label_segment_color = "grey20",
  label_segment_size = 0.3,
  candidate_node_alpha = 0.7,
  edge_mode = "candidate_connected",
  label_mode = "candidate_hub",
  label_size = 1.8,
  min_segment_length = 0,
  label_box_padding = 0.5
)

netcomi_p_all_0.3_type2 <- plot_netcomi_network(
  layout_df = netcomi_layout_list_0.3$all,
  edge_alpha = 0.2,
  background_node_alpha = 0.15,
  label_segment_color = "grey20",
  label_segment_size = 0.3,
  candidate_node_alpha = 0.7,
  edge_mode = "candidate_candidate",
  label_mode = "candidate_hub",
  label_size = 1.8,
  min_segment_length = 0,
  label_box_padding = 0.5
)

# ggplot2::ggsave(
#   filename = "ggraph_high_res_network5.pdf",
#   plot = netcomi_p_all_0.3_type1,
#   width = 10,
#   height = 10,
#   units = "in",
#   device = "pdf"
# )

save(list = c("netcomi_igraph_list_0.3", "netcomi_node_table_0.3", "netcomi_network_table_0.3"),
     file = "rdata/netcomi_extract.Rdata")
save(list = c("node_metric_tests_0.3"),
     file = "rdata/netcomi_stat.Rdata")
save(list = c("netcomi_box_plots_0.3","netcomi_p_all_0.3_type1","netcomi_p_all_0.3_type2"),
     file = "rdata/netcomi_plots.Rdata")

rm(list = ls())
gc()
