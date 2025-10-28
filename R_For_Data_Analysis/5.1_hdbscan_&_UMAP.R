#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("rdata/prm.Rdata")
load("rdata/aes.Rdata")
load("rdata/RPCA_90.Rdata") 
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)      
  library(dbscan)
  library(RColorBrewer)
})
seed <- prm.ls$general$seed_all
set.seed(seed)
#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/plot_rpca_biplot.R")
source("functions/run_hdbscan_rpca.R")
source("functions/plot_rpca_biplot_cluster.R")
source("functions/add_hclust_auto_layers_to_RPCA.R")
source("functions/plot_umap_cluster.R")
source("functions/run_umap_hdbscan_rpca.R")
source("functions/run_umap_visual_rpca.R")
source("functions/plot_umap_visual.R")
#-------------------------------------------------------------------------------
# User options
#-------------------------------------------------------------------------------
save_file_5.1 <- "rdata/hdbscan_&_UMAP_90.Rdata"
#norm.method <- prm.ls$pcoa$norm_m
dists <- prm.ls$pcoa$distance
dataset_names = c("_all","dtt1","t1","dtt1ts", 
                  "t1_T0","t1_T24","t1_T48","t1_T72")


#-------------------------------------------------------------------------------
# HBDSCAN on RPCA biplot
#-------------------------------------------------------------------------------

RPCA <- RPCA_filt
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "_all",
                         n_pcs = 10,minPts = 4,plot = FALSE)
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "dtt1",
                         n_pcs = 10,minPts = 4,plot = FALSE)
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "t1",
                              n_pcs = 15,minPts = 8,plot = FALSE)
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "dtt1ts",
                         n_pcs = 10,minPts = 3,plot = FALSE)
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "t1_T0",
                             n_pcs = 10,minPts = 4,plot = FALSE)
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "t1_T24",
                         n_pcs = 10,minPts = 4,plot = FALSE)
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "t1_T48",
                         n_pcs = 10,minPts = 4,plot = FALSE)
RPCA <- run_hdbscan_rpca(rpca_results = RPCA,dataset_name = "t1_T72",
                            n_pcs = 10,minPts = 4,plot = FALSE)

# ----------------------- UMAP + HDBSCAN ---------------------------------------
UMAP_ls <- list()
UMAP_ls$`_all` <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "_all",
                                     n_pcs = 8, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 20, 
                                     min_dist = 0.15, minPts = 6) 
UMAP_ls$t1 <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "t1",
                                   n_pcs = 10, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 20, 
                                   min_dist = 0.15, minPts = 6) # min cluster size
UMAP_ls$dtt1 <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "dtt1",
                                         n_pcs = 8, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 20, 
                                         min_dist = 0.15, minPts = 6) 
UMAP_ls$dtt1ts <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "dtt1ts",
                                       n_pcs = 8, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 10, 
                                       min_dist = 0.1, minPts = 4)
UMAP_ls$t1_T0 <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T0",
                                         n_pcs = 8, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 15, 
                                         min_dist = 0.15, minPts = 4)
UMAP_ls$t1_T24 <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T24",
                                        n_pcs = 8, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 15, 
                                        min_dist = 0.15, minPts = 4)
UMAP_ls$t1_T48 <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T48",
                                        n_pcs = 8, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 15, 
                                        min_dist = 0.15, minPts = 4)
UMAP_ls$t1_T72 <-  run_umap_hdbscan_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T72",
                                        n_pcs = 8, n_umap = 6, n_hdbscan_dim = 6, n_neighbors = 15, 
                                        min_dist = 0.15, minPts = 4)

#-------------------------------------------------------------------------------
# Add hierarchical clustering layers to RPCA biplot tables and UMAP
#-------------------------------------------------------------------------------
# RPCA
RPCA <- add_hclust_to_HDBSCAN(RPCA, pcs = paste0("PC", 1:15),
                              min_super = 2, min_sub = 7)
# UMAP
UMAP_ls <- add_hclust_to_HDBSCAN(UMAP_ls, pcs = paste0("PC", 1:15),
                              min_super = 2, min_sub = 7)


# -------------------- plot RPCA biplot ----------------------------------
# color 
my_colors <- brewer.pal(10, "Paired")
lv <- levels(RPCA[["dtt1"]]$biplot$Cluster)
n_clusters <- sum(lv != "Noise")
cluster_colors_vec <- my_colors[1:n_clusters]
if ("Noise" %in% lv) {
  names(cluster_colors_vec) <- lv[lv != "Noise"]
  cluster_colors_vec <- c(cluster_colors_vec, Noise = "grey70")
} else {
  names(cluster_colors_vec) <- lv
}
cluster_colors_vec

# make plots
p_rpca_cluster_t1 <- plot_rpca_biplot_cluster(
  rpca_results = RPCA, dataset_name = "t1",
  color_col = "Project", shape_col  = "TimePoint",
  add_species = FALSE, add_ellipse = FALSE, add_cluster = FALSE,# 不用 cluster 上色
  show_cluster_ellipse = FALSE,    # 只圈出 cluster 边界
  color_vec = aes.ls$col_alt$Project, shape_vec = aes.ls$shape$TimePoint)
p_rpca_cluster_t1

p_rpca_cluster_dtt1ts <- plot_rpca_biplot_cluster(
  rpca_results = RPCA, dataset_name = "dtt1ts",
  color_col = "Project", shape_col  = "TimePoint",
  add_species = FALSE, add_ellipse = FALSE, add_cluster = FALSE,# 不用 cluster 上色
  show_cluster_ellipse = FALSE,    # 只圈出 cluster 边界
  color_vec = aes.ls$col_alt$Project, shape_vec = aes.ls$shape$TimePoint)
p_rpca_cluster_dtt1ts

p_rpca_cluster_all <- plot_rpca_biplot_cluster(
  rpca_results = RPCA,
  dataset_name = "_all",
  color_col = "Project",          
  shape_col  = "TimePoint",
  add_species = FALSE,            
  add_ellipse = FALSE,          
  add_cluster = FALSE,            
  show_cluster_ellipse = TRUE,    
  color_vec = aes.ls$col_alt$Project,
  shape_vec = aes.ls$shape$TimePoint
)
p_rpca_cluster_all

p_rpca_cluster_t1_T0 <- plot_rpca_biplot_cluster(
  rpca_results = RPCA, dataset_name = "t1_T0",
  color_col = "Project", shape_col  = NULL,
  add_species = FALSE, top_n_species = 20,
  add_ellipse = FALSE, add_cluster = FALSE,
  show_cluster_ellipse = TRUE,   
  color_vec = aes.ls$col_alt$Project, shape_vec = aes.ls$shape$TimePoint)
p_rpca_cluster_t1_T0

p_rpca_cluster_t1_T24 <- plot_rpca_biplot_cluster(
  rpca_results = RPCA, dataset_name = "t1_T24",
  color_col = "Project", shape_col  = NULL,
  add_species = FALSE, add_ellipse = FALSE, add_cluster = FALSE,
  show_cluster_ellipse = TRUE,    
  color_vec = aes.ls$col_alt$Project, shape_vec = aes.ls$shape$TimePoint)
p_rpca_cluster_t1_T24

p_rpca_cluster_t1_T48 <- plot_rpca_biplot_cluster(
  rpca_results = RPCA, dataset_name = "t1_T48",
  color_col = "Project", shape_col  = NULL,
  add_species = FALSE, add_ellipse = FALSE, add_cluster = FALSE,# 不用 cluster 上色
  show_cluster_ellipse = TRUE,    
  color_vec = aes.ls$col_alt$Project, shape_vec = aes.ls$shape$TimePoint)
p_rpca_cluster_t1_T48

p_rpca_cluster_t1_T72 <- plot_rpca_biplot_cluster(
  rpca_results = RPCA, dataset_name = "t1_T72",
  color_col = "Project", shape_col  = NULL,
  add_species = FALSE, add_ellipse = FALSE, add_cluster = FALSE,# 不用 cluster 上色
  show_cluster_ellipse = TRUE,    
  color_vec = aes.ls$col_alt$Project, shape_vec = aes.ls$shape$TimePoint)
p_rpca_cluster_t1_T72

# save(list = c("RPCA", "p_rpca_cluster_dtt1", "p_rpca_cluster_all","p_rpca_cluster_dtt1_T0",
#               "p_rpca_cluster_dtt1_T24", "p_rpca_cluster_dtt1_T48", "p_rpca_cluster_dtt1_T72"),
#      file = "rdata/hdbscan_rpca_result.Rdata")

#----------------------- plot UMAP ---------------------------------------------
UMAP_v <- list()
UMAP_v$`_all` <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "_all",
                                      n_pcs = 8, n_umap = 6,n_neighbors = 20, min_dist = 0.15)
UMAP_v$t1 <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "t1",
                                  n_pcs = 10, n_umap = 6,n_neighbors = 20, min_dist = 0.15)
UMAP_v$dtt1 <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "dtt1",
                                    n_pcs = 8, n_umap = 6,n_neighbors = 20, min_dist = 0.15)
UMAP_v$dtt1ts <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "dtt1ts",
                                      n_pcs = 8, n_umap = 6,n_neighbors = 10, min_dist = 0.1)
UMAP_v$t1_T0 <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T0",
                                     n_pcs = 8, n_umap = 6,n_neighbors = 15, min_dist = 0.15)
UMAP_v$t1_T24 <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T24",
                                      n_pcs = 8, n_umap = 6,n_neighbors = 15, min_dist = 0.15)
UMAP_v$t1_T48 <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T48",
                                      n_pcs = 8, n_umap = 6,n_neighbors = 15, min_dist = 0.15)
UMAP_v$t1_T72 <- run_umap_visual_rpca(rpca_results = RPCA_filt, dataset_name = "t1_T72",
                                      n_pcs = 8, n_umap = 6,n_neighbors = 15, min_dist = 0.15)

UMAP_NAMES <- names(UMAP_v)

UMAP_plots <- lapply(UMAP_NAMES, function(ds) {
  p <- plot_umap_visual (rpca_results = UMAP_v, dataset_name = ds,
                         color_col = "Project", shape_col  = "TimePoint",
                         color_vec = aes.ls$col_alt$Project,
                         shape_vec = aes.ls$shape$TimePoint,
                         title_prefix = NULL)
  print(p)
  return(p)
})
names(UMAP_plots) <- UMAP_NAMES

# UMAP_NAMES <- names(UMAP_ls)
# 
# UMAP_plots <- lapply(UMAP_NAMES, function(ds) {
#   p <- plot_umap_cluster(rpca_results = UMAP_ls, dataset_name = ds,
#                          color_col = "Project", shape_col  = "TimePoint",
#                          color_vec = aes.ls$col_alt$Project,
#                          shape_vec = aes.ls$shape$TimePoint)
#   return(p)
# })
# names(UMAP_plots) <- UMAP_NAMES
UMAP_plots$t1


save(list = c("RPCA", "UMAP_ls", "UMAP_v", "UMAP_plots"), 
     file = save_file_5.1)

rm(list = ls())
gc()




