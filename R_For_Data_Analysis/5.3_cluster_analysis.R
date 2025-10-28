#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("rdata/prm.Rdata")
load("rdata/rdata/hdbscan_&_UMAP_filt.Rdata")
load("rdata/aes.Rdata")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggalluvial)
  library(stringr)
  library(rlang)
  library(vegan)
})

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/build_cluster_timeline.R")
source("functions/plot_cluster_sankey.R")
source("functions/plot_cluster_hierarchy.R")
source("functions/analyze_meta_effects_cluster.R")
source("functions/classify_batch_behavior.R")
source("functions/plot_behavior_sankey.R") 
source("functions/classify_cluster_behavior.R")
#-------------------------------------------------------------------------------
# extract data
#-------------------------------------------------------------------------------
data_set <- UMAP_ls$t1_T72
biplot_table <- data_set$biplot
# set RUN as factor
as.factor(biplot_table$RUN)

# rearrange the levels of Cluster
cl_chr <- trimws(as.character(biplot_table$Cluster)) 
lvls   <- unique(cl_chr)
is_num <- grepl("^[0-9]+$", lvls)
num_labels <- lvls[is_num]
num_labels <- num_labels[order(as.numeric(num_labels))]
non_num_labels <- lvls[!is_num]
nn_norm <- tolower(non_num_labels)
tail_order <- c("outlier", "unassigned", "unknown", "noise")
tail_rank  <- match(nn_norm, tail_order, nomatch = 0)
non_special <- non_num_labels[tail_rank == 0]
special    <- non_num_labels[tail_rank > 0]
special    <- special[order(tail_rank[tail_rank > 0])]
new_levels <- c(num_labels, non_special, special)
biplot_table$Cluster <- factor(cl_chr, levels = new_levels, ordered = TRUE)
levels(biplot_table$Cluster)

#-------------------------------------------------------------------------------
# User options
#-------------------------------------------------------------------------------
seed <- prm.ls$general$seed_all
set.seed(seed)

time_colors <- c(T0 = "#264653",T24 = "#2A9D8F",T48 = "#E9C46A",T72 = "#E76F51")

cluster_colors <- c(
  "#264653", "#2A9D8F", "#8AB17D", "#E9C46A", "#F4A261",
  "#E76F51", "#6D597A", "#B56576", "#E56B6F", "#EAAC8B",
  "#457B9D", "#A8DADC", "#EAE7DC", "#E63946", "#9B2226",
  "#003049", "#669BBC", "#EAE2B7", "#F77F00", "#D62828",
  "#D4A373", "#B07C59", "#D1CBC1", "#A7C6C9"
)

RUN_colors <- c(
  "#4A6FA5", "#73A9AD", "#A7C6C9", "#D1CBC1", "#E8D8C3",
  "#D4A373", "#B07C59", "#855D44", "#5B4C43", "#3A3A3A",
  "#4E6E81", "#7B9E89", "#A5BFA0", "#C5C1A9", "#E0D8C3",
  "#C9986C", "#A06A4E", "#73564B", "#555555", "#2F2F2F"
)


# cluster_names <- as.character(sort(unique(cluster_behavior_df$Cluster)))
# cluster_colors <- setNames(cluster_colors[seq_along(cluster_names)], cluster_names)

#-------------------------------------------------------------------------------
# Summarize cluster change along time 
#-------------------------------------------------------------------------------
cluster_timeline <- build_cluster_timeline(biplot_table, batch_col = "Batch_ID", 
                       time_col = "TimePoint", cluster_col = "Cluster",
                       meta_cols = c("Project","RUN"))

time_cols <- grep("^T", names(cluster_timeline$wide), value = TRUE)

cluster_timeline$wide %>%
  rowwise() %>%
  mutate(stable = length(unique(c_across(all_of(time_cols)))) == 1) %>%
  ungroup() %>%
  summarise(prop_stable = mean(stable))


#-------------------------------------------------------------------------------
# Plot Sankey diagram showing cluster transitions over time
#-------------------------------------------------------------------------------
plot.cluster.sankey <- plot_cluster_sankey(cluster_timeline$long, 
                                           batch_col = "Batch_ID")
plot.cluster.sankey


#------------ cluster-project summary --------------
plot_cluster_hierarchy(
  biplot_tbl = biplot_table,
  cluster_long = cluster_timeline$long,
  fill_col = "Project",
  hc          = data_set$hclust_cluster,
  fill_colors = aes.ls[["col_alt"]][["Project"]]
)

#------------ cluster-RUN summary --------------
plot_cluster_hierarchy(
  biplot_tbl = biplot_table,
  cluster_long = cluster_timeline$long,
  fill_col = "RUN",
  fill_colors = RUN_colors,
  hc          = data_set$hclust_cluster,
  title = "Sample distribution by RUN"
)

#------------ cluster-time summary --------------
plot_cluster_hierarchy(
  biplot_tbl = biplot_table,
  cluster_long = cluster_timeline$long,
  fill_col = "TimePoint",
  fill_colors = time_colors,
  hc          = data_set$hclust_cluster,
  title = "Sample distribution by Time"
)

1#------------run LASSO multinomial logistic regression -------------------
# TimePoint, Project, RUN
res_lasso <- analyze_meta_effects_cluster(biplot_table, min_cluster_size = 8)

#------------ check bath_ID behavior in cluster during time --------------
behavior_df <- classify_batch_behavior(cluster_timeline$wide)

# summarize and plot the results
behavior_summary <- behavior_df %>%
  as.data.frame(stringsAsFactors = FALSE) %>% 
  dplyr::select(where(~!is.list(.))) %>%                     
  mutate(across(everything(), as.character)) %>%             
  dplyr::count(behavior, Project, name = "n") %>%            
  group_by(behavior) %>%
  mutate(
    total_n = sum(n),
    prop = n / total_n
  ) %>%
  ungroup()


ggplot(behavior_summary, aes(x = behavior, y = n, fill = Project)) +
  geom_col(color = "grey30") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Batch behavior composition by Project",
    x = "Behavior category",
    y = "Number of Batch_IDs",
    fill = "Project"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1))

# sankey for special Batch_ID
p <- plot_behavior_sankey_T0_T72(
  behavior_df = behavior_df,
  cluster_long = cluster_timeline$long,
  behaviors_to_include = c("Converging", "Dynamic transition"),
  title = "Cluster transitions (T0 → T72): Converging & Dynamic transition"
)

#------- check cluster behavior over time ------------------------------

cluster_behavior_df <- classify_cluster_behavior(cluster_timeline$long, 
                                                 epsilon = 0)


behavior_summary <- cluster_behavior_df %>%
  group_by(behavior) %>%
  summarise(n_cluster = n(), .groups = "drop")

cluster_behavior_summary <- cluster_behavior_df %>%
  group_by(behavior, Cluster) %>%
  summarise(n = 1, .groups = "drop")  

ggplot(cluster_behavior_summary,
       aes(x = behavior, fill = Cluster)) +
  geom_bar(position = "stack", color = "grey30") +
  scale_fill_manual(values = cluster_colors) + 
  theme_minimal(base_size = 14) +
  labs(
    title = "Cluster behavior distribution",
    x = "Behavior type",
    y = "Number of clusters",
    fill = "Cluster"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

