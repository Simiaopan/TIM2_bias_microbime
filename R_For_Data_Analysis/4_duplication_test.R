#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("rdata/prm.Rdata")
# load("rdata/data.Rdata") # meta and phyloseq
# load("output/beta_plots.Rdata")
load("rdata/RPCA_90.Rdata")  #change when process new datasets
load("rdata/aes.Rdata")
suppressPackageStartupMessages({
  library(ggplot2)
  library(phyloseq)
  library(lme4)
  library(lmerTest)
  library(dplyr)
  library(merTools)
})

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/extract_rpca_coords_meta.R")
source("functions/dists_within_out_dupli.R")
source("functions/plot_violin_dupli.R")
source("functions/run_lmm_dupli.R")
source("functions/plot_lmm_clr_dupli.R")
source("functions/duplication_permutation_test.R")
# source("functions/plot_duplication_permutation.R")
#-------------------------------------------------------------------------------
# User options
#-------------------------------------------------------------------------------
seed <- prm.ls$general$seed_all
set.seed(seed)


#-------------------------------------------------------------------------------
# extract data
#-------------------------------------------------------------------------------
rpca_dtt1 <- extract_rpca_coords_meta(RPCA_filt, "dtt1")
coords_rpca_dtt1 <- rpca_dtt1$coords_all
meta_rpca_dtt1 <- rpca_dtt1$meta

# coords_dtt1_clr <- pcoa_clr_dtt1$coords_all
# meta_dtt1 <- as(sample_data(data.ls$Ps$ASV$CLR$dtt1), "data.frame")

#-------------------------------------------------------------------------------
# run duplication check, and make violin plot
#-------------------------------------------------------------------------------
comp_dupli_rpca <- dists_within_out_dupli(coords_rpca_dtt1, meta_rpca_dtt1,
                                         treat_col = "Treat_ID",
                                         batch_col = "Batch_ID",
                                         tp_col = "TimePoint",time_col = "Time")

# point mode
plot.dupli.rpca <- plot_violin_dupli(comp_dupli_rpca, aes.ls, mode = "point",
                    show_box = TRUE, time_levels = c("T0", "T24", "T48", "T72"))

plot.dupli.rpca
#-------------------------------------------------------------------------------
# run LMM on within/out duplication pairs
#-------------------------------------------------------------------------------
# lmm_clr_dupli <- run_lmm_dupli(comp_dupli_clr, do_fdr = TRUE)
# lmm_in_out_table <- lmm_clr_dupli$fixed

lmm_rpca_dupli <- run_lmm_dupli(comp_dupli_rpca$distances, do_fdr = TRUE)
lmm_in_out_table_rpca <- lmm_rpca_dupli$fixed

# plot data
# plot_lmm_clr_dupli(lmm_rpca_dupli)

ggplot(lmm_rpca_dupli$pred, aes(x = Time, y = fit, color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  scale_color_manual(values = c("Within" = "#E76F51", "Between" = "#2A9D8F")) +
  scale_fill_manual(values = c("Within" = "#E76F51", "Between" = "#2A9D8F")) +
  scale_x_continuous(breaks = c(0, 24, 48, 72), limits = c(0, 72)) +
  scale_y_continuous(limits = c(0, 0.35), expand = c(0, 0)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black"),
        aspect.ratio = 0.9) +
  labs(x = "Time (h)", y = "Predicted Distance", color = "Group", fill = "Group")


#-------------------------------------------------------------------------------
# run duplication permutation
#-------------------------------------------------------------------------------
# permutaion_clr_sample <- duplication_permutation_test(coords_dtt1_clr, meta_dtt1,
#                            mode = "sample", n_iter = 300)
# 
# permutation_clr_centroid <- duplication_permutation_test(coords_dtt1_clr, meta_dtt1,
#                                                mode = "centroid",
#                                                n_iter = 300)
# 
# plot.permutation.clr.centroid <- plot_duplication_permutation(check_clr_centroid,
#                                                       title = "Stability (centroid mode)")
# 
# plot_violin_dupli(comp_dupli_rpca, aes.ls, mode = "point",
#                time_levels = c("T0","T24","T48","T72"))




