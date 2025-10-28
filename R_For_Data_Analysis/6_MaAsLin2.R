#===============================================================================
#  MaAsLin2-based microbiome LMM analysis
#===============================================================================

load("rdata/prm.Rdata")
load("rdata/filt_data.Rdata")

seed = prm.ls$general$seed_all
set.seed(seed)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ggplot2)
  library(glue)
  library(Maaslin2)
  library(lmtest)
})

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/check_abundance_distribution.R")
source("functions/filter_microbiome_data_long.R")
source("functions/run_MaAsLin2_LMM.R")   
source("functions/analyze_time_trends.R")
source("functions/merge_all_sub_results.R")
source("functions/diagnose_MaAsLin2_fit.R")
source("functions/plot_volcano_M.R")
#-------------------------------------------------------------------------------
# extract data
#-------------------------------------------------------------------------------
subdata_raw <- data.ls_filt$RAW$t1
subdata_clr <- data.ls_filt$CLR$t1

# RAW data (ASV)
ps_raw_asv <- data.ls_filt$RAW$`_all`
asv_raw <- psmelt(ps_raw_asv) %>%
  dplyr::select(OTU, Abundance, Genus, Sample, Species, 
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data (ASV)
ps_clr_asv <- data.ls_filt$CLR$`_all`
asv_clr <- psmelt(ps_clr_asv) %>%
  dplyr::select(OTU, Abundance, Genus, Sample, Species, 
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# RAW data (genus)
ps_raw_genus <- tax_glom(ps_raw_asv, taxrank = "Genus")
genus_raw <- psmelt(ps_raw_genus) %>%
  dplyr::select(Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data (genus)
ps_clr_genus <- tax_glom(ps_clr_asv, taxrank = "Genus")
genus_clr <- psmelt(ps_clr_genus) %>%
  dplyr::select(Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# RAW data sub (ASV)
ps_raw_asv_sub <- subdata_raw
asv_raw_sub <- psmelt(ps_raw_asv_sub) %>%
  dplyr::select(OTU, Abundance, Genus, Sample, Species, 
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data sub (ASV)
ps_clr_asv_sub <- subdata_clr
asv_clr_sub <- psmelt(ps_clr_asv_sub) %>%
  dplyr::select(OTU, Abundance, Genus, Sample, Species, 
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# RAW data sub (genus)
ps_raw_genus_sub <- tax_glom(subdata_raw, taxrank = "Genus")
genus_raw_sub <- psmelt(ps_raw_genus_sub) %>%
  dplyr::select(Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data sub (genus)
ps_clr_genus_sub <- tax_glom(subdata_clr, taxrank = "Genus")
genus_clr_sub <- psmelt(ps_clr_genus_sub) %>%
  dplyr::select(Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

#-------------------------------------------------------------------------------
# check and filter
#-------------------------------------------------------------------------------
asv_raw_filtered_list <- filter_microbiome_data_long(asv_raw, level = "OTU",
                                                     abundance_col = "Abundance",
                                                     sample_col = "Sample",
                                                     prev_cutoff = 0, mean_cutoff = 0)
asv_raw <- asv_raw_filtered_list$filtered_df

genus_raw_filtered_list <- filter_microbiome_data_long(genus_raw, level = "Genus",
                                                       abundance_col = "Abundance",
                                                       sample_col = "Sample",
                                                       prev_cutoff = 0, mean_cutoff = 0)
genus_raw <- genus_raw_filtered_list$filtered_df

asv_raw_sub_filtered_list <- filter_microbiome_data_long(asv_raw_sub, level = "OTU",
                                                         abundance_col = "Abundance",
                                                         sample_col = "Sample",
                                                         prev_cutoff = 0, mean_cutoff = 0)
asv_raw_sub <- asv_raw_sub_filtered_list$filtered_df

genus_raw_sub_filtered_list <- filter_microbiome_data_long(genus_raw_sub, level = "Genus",
                                                           abundance_col = "Abundance",
                                                           sample_col = "Sample",
                                                           prev_cutoff = 0, mean_cutoff = 0)
genus_raw_sub <- genus_raw_sub_filtered_list$filtered_df

#===============================================================================
#------------------------ run MaAsLin2 LMM analysis ----------------------------
#===============================================================================
# ASV (all)
inp_asv <- prep_maaslin_inputs(asv_raw, "OTU")
res_asv_time_proj <- run_MaAsLin2_LMM(
  inp_asv$feat, inp_asv$meta,
  output_dir = "results/maaslin2_asv_all",
  fixed_effects = c("Time_c"),
  random_effects = "Batch_ID"
)
res_asv_time_proj_all <- res_asv_time_proj$result

# GENUS (all)
inp_genus <- prep_maaslin_inputs(genus_raw, "Genus")
res_genus_time_proj <- run_MaAsLin2_LMM(
  inp_genus$feat, inp_genus$meta,
  output_dir = "results/maaslin2_genus_all",
  fixed_effects = c("Time_c"),
  random_effects = "Batch_ID"
)
res_genus_time_proj_all <- res_genus_time_proj$result

# ASV (subset)
inp_asv_sub <- prep_maaslin_inputs(asv_raw_sub, "OTU")
res_asv_time_proj_sub <- run_MaAsLin2_LMM(
  inp_asv_sub$feat, inp_asv_sub$meta,
  output_dir = "results/maaslin2_asv_sub",
  fixed_effects = c("Time_c"),
  random_effects = "Batch_ID"
)
res_asv_time_proj_sub_all <- res_asv_time_proj_sub$result

# GENUS (subset)
inp_genus_sub <- prep_maaslin_inputs(genus_raw_sub, "Genus")
res_genus_time_proj_sub <- run_MaAsLin2_LMM(
  inp_genus_sub$feat, inp_genus_sub$meta,
  output_dir = "results/maaslin2_genus_sub",
  fixed_effects = c("Time_c"),
  random_effects = "Batch_ID"
)
res_genus_time_proj_sub_all <- res_genus_time_proj_sub$result


#-------------------------------------------------------------------------------
# add more info to evaluate
#-------------------------------------------------------------------------------
out_genus_sub <- analyze_time_trends(res_df  = res_genus_time_proj_sub_all, 
                                     sub_df  = genus_raw_sub, feature_col = "Genus")
genus_summary_sub_M <- out_genus_sub$merged

out_asv_sub <- analyze_time_trends(res_df  = res_asv_time_proj_sub_all,
                                   sub_df  = asv_raw_sub,feature_col = "OTU")
asv_summary_sub_M <- out_asv_sub$merged

out_genus_all <- analyze_time_trends(res_df  = res_genus_time_proj_all,  
                                     sub_df  = genus_raw_sub, feature_col = "Genus")
genus_summary_all_M <- out_genus_all$merged

out_asv_all <- analyze_time_trends(res_df  = res_asv_time_proj_sub_all,  
                                     sub_df  = asv_raw_sub, feature_col = "OTU")
asv_summary_all_M <- out_asv_all$merged

# merge sub and all
merged_genus_M <- merge_all_sub_results(res_all = genus_summary_all_M,
                  res_sub = genus_summary_sub_M, feature_col = "feature")

merged_asv_M <- merge_all_sub_results(res_all = asv_summary_all_M,
                  res_sub = asv_summary_sub_M, feature_col = "feature")

#---------------- diagnose and add info ----------------------------------------
# asv(all)
res_diag_asv <- diagnose_MaAsLin2_fit(fit = res_asv_time_proj$fit,
                merged_df = merged_asv_M$all_with_both, feature_col = "feature")

# genus(all)
res_diag_genus <- diagnose_MaAsLin2_fit(fit = res_genus_time_proj$fit,
                  merged_df = merged_genus_M$all_with_both, feature_col = "feature")

# asv(sub)
res_diag_genus_sub <- diagnose_MaAsLin2_fit(fit = res_asv_time_proj_sub$fit,
                  merged_df = merged_asv_M$sub_with_both, feature_col = "feature")
# genus(sub)
res_diag_genus_sub <- diagnose_MaAsLin2_fit(fit = res_genus_time_proj_sub$fit,
                                        merged_df = merged_genus_M$sub_with_both, feature_col = "feature")
#-------------------------- plot results ---------------------------------------
volcano_all_asv <- plot_volcano_M(res_diag_asv$merged_out)
volcano_genus_asv <- plot_volcano_M(res_diag_genus$merged_out)


feature_asv_all <- sample_remarkable_features(res_diag_asv$merged_out, 
                                              n = 6, m = 6, k = 0.00039)
feature_genus_all <- sample_remarkable_features(res_diag_asv$merged_out, 
                                                n = 10, m = 10, k = 0.001)

plot_feature_diagnostics(fit = res_asv_time_proj$fit,
                         features = feature_asv_all$Remarkable_Suspect_Sampling,
                         ncol = 3)
plot_feature_diagnostics(fit = res_asv_time_proj$fit,
                         features = feature_asv_all$Remarkable_Reliable_Sampling,
                         ncol = 3)
