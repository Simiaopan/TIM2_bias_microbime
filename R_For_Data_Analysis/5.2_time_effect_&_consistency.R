#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("rdata/prm.Rdata")
load("rdata/hdbscan_&_UMAP_filt.Rdata") # change when process different datasets
load("rdata/aes.Rdata")

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(tidyr)
  library(emmeans)
  library(ggplot2)
  library(car)
  library(pheatmap)
  library(vegan)
})

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/run_rpca_lmm.R")
source("functions/visualize_adonis2.R")
source("functions/test_time_consistency.R")
source("functions/plot_pc_time_trends.R")
source("functions/extract_manova_results.R")
source("functions/visualize_manova_eta2.R")
source("functions/check_pc_lmm.R")
#-------------------------------------------------------------------------------
# user options
#-------------------------------------------------------------------------------
data_set <- RPCA$t1
biplot_table <- data_set$biplot
save_file <- "rdata/time_effect_90_t1.Rdata"
color_vec <- c( "#E76F51", "#E9C46A", "#2A9D8F")


# ------- LMM for each selected-rpca-axis --------------------------------------
res_time_numeric <- run_rpca_lmm(biplot_df = biplot_table, time_var = "Time",
                                 project_var = "Project", random_var = "Batch_ID",
                                 n_pcs = 5)

# ------- adonis2 for time and time*project ------------------------------------
dist_rpca <- dist(biplot_table[, c("PC1", "PC2", "PC3", "PC4","PC5")])
adonis_res<- adonis2(dist_rpca ~ TimePoint + Project + TimePoint:Project, 
                      data = biplot_table,permutations = 999,by = "terms")

adonis_p <- visualize_adonis2(adonis_res, color_vec = color_vec)

adonis_results <- list(adonis_res, adonis_p)

# ------- test time consistency: LMM + MANOVA + cosine similarity --------------
res_time_consistency <- test_time_consistency(biplot_table, lmm_result = res_time_numeric,
            time_var = "Time", project_var = "Project", random_var = "Batch_ID",
            pcs_for_manova = 1:5)

# plot LMM-time-proj in each PC
PC_time_plot <- plot_pc_time_trends  (res_time_consistency)
PC_time_plot # nSignif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ others: ns"

# visualize manova results
manova_results <- extract_manova_results(res_time_consistency$manova)
    # chek the data.frame before move to the next step!
manova_plot <- visualize_manova_eta2(pillai_df, color_vec = color_vec)

# check pc-lmm models
pcs_to_check <- c("PC1", "PC2", "PC3", "PC4", "PC5")

pc_check_plot <- check_pc_lmm(data = biplot_table,pc_names = pcs_to_check,
                  time_var = "Time",project_var = "Project",random_var = "Batch_ID")
pc_check_plot



save(list = c("res_time_numeric", "adonis_results", "res_time_consistency", 
              "manova_plot", "PC_time_plot", "pc_check_plot"),
     file = save_file)

rm(list = ls())
gc()
