#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------

load("rdata/prm.Rdata")
#load("rdata/data.Rdata")
load("rdata/filt_data.Rdata")

seed = prm.ls$general$seed_all
set.seed(seed)

suppressPackageStartupMessages({
  library(phyloseq)
  library(lme4)
  library(dplyr)
  library(glmmTMB)
  library(broom.mixed)
  library(ggplot2)
  library(lmerTest)
  library(glue)
  library(DHARMa)
})

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/check_abundance_distribution.R")
source("functions/filter_microbiome_data_long.R")
source("functions/run_glmmTMB_LMM_plus.R")
source("functions/run_glmmTMB_LMM.R")
source("functions/plot_GLMM_LMM_volcano.R")
#-------------------------------------------------------------------------------
# extract data
#-------------------------------------------------------------------------------
subdata_raw <- data.ls_filt$RAW$t1
subdata_clr <- data.ls_filt$CLR$t1
# RAW data (ASV)
# ps_raw_asv <- data.ls$Ps$ASV$RAW$`_all`
ps_raw_asv <- data.ls_filt$RAW$`_all`
asv_raw <- ps_raw_asv %>%
  psmelt(.) %>%
  dplyr::select(.,OTU, Abundance, Genus, Sample, Species, 
         Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data (ASV)
#ps_asv_clr <- data.ls$Ps$ASV$CLR$`_all`
ps_clr_asv <- data.ls_filt$CLR$`_all`
asv_clr <- ps_clr_asv %>%
  psmelt() %>%
  dplyr::select(.,OTU, Abundance, Genus, Sample, Species, 
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)


# RAW data (genus)
ps_raw_genus <- tax_glom(ps_raw_asv, taxrank = "Genus")
genus_raw <- ps_raw_genus %>%
  psmelt(.) %>%
  dplyr::select(.,Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data (genus)
ps_clr_genus <- tax_glom(ps_clr_asv, taxrank = "Genus")
genus_clr <- ps_clr_genus %>%
  psmelt() %>%
  dplyr::select(.,Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# RAW data sub (ASV)
ps_raw_asv_sub <- subdata_raw
asv_raw_sub <- ps_raw_asv_sub %>%
  psmelt(.) %>%
  dplyr::select(.,OTU, Abundance, Genus, Sample, Species, 
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data sub (ASV)
ps_clr_asv_sub <- subdata_clr
asv_clr_sub <- ps_clr_asv_sub %>%
  psmelt() %>%
  dplyr::select(.,OTU, Abundance, Genus, Sample, Species, 
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# RAW data sub (genus)
ps_raw_genus_sub <- tax_glom(subdata_raw, taxrank = "Genus")
genus_raw_sub <- ps_raw_genus_sub %>%
  psmelt(.) %>%
  dplyr::select(.,Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

# clr data sub (genus)
ps_clr_genus_sub <- tax_glom(subdata_clr, taxrank = "Genus")
genus_clr_sub <- ps_clr_genus_sub %>%
  psmelt() %>%
  dplyr::select(.,Genus, Abundance, Sample,
                Time, Project, RUN, TimePoint, TimePattern, Batch_ID, Treat_ID)

#-------------------------------------------------------------------------------
# check and filter
#-------------------------------------------------------------------------------
# distribution_check_raw <- check_abundance_distribution(asv_raw, level = "OTU", 
#                                                        value_col = "Abundance")

# asv_raw_filtered_list <- filter_microbiome_data_long(asv_raw, level = "OTU", 
#                             abundance_col = "Abundance", sample_col = "Sample",
#                             prev_cutoff = 0, mean_cutoff = 1e-4)
# asv_raw_filtered <- asv_raw_filtered_list$filtered_df

# distribution_check_filtered <- check_abundance_distribution(asv_raw_filtered, 
#                                       level = "OTU", value_col = "RelAbundance")

# asv_clr_filtered <- asv_clr %>%
#   filter(OTU %in% asv_raw_filtered_list$kept_ids)

#all_asv
asv_raw_filtered_list <- filter_microbiome_data_long(asv_raw, level = "OTU",
                            abundance_col = "Abundance", sample_col = "Sample",
                            prev_cutoff = 0, mean_cutoff = 0)
asv_raw <- asv_raw_filtered_list$filtered_df

# all_genus
genus_raw_filtered_list <- filter_microbiome_data_long(genus_raw, level = "Genus", 
                                                     abundance_col = "Abundance", sample_col = "Sample",
                                                     prev_cutoff = 0, mean_cutoff = 0)
genus_raw <- genus_raw_filtered_list$filtered_df

# sub_asv
asv_raw_sub_filtered_list <- filter_microbiome_data_long(asv_raw_sub, level = "OTU", 
                                                     abundance_col = "Abundance", sample_col = "Sample",
                                                     prev_cutoff = 0, mean_cutoff = 0)
asv_raw_sub <- asv_raw_sub_filtered_list$filtered_df

# sub_genus
genus_raw_sub_filtered_list <- filter_microbiome_data_long(genus_raw_sub, level = "Genus", 
                                                         abundance_col = "Abundance", sample_col = "Sample",
                                                         prev_cutoff = 0, mean_cutoff = 0)
genus_raw_sub <- asv_raw_sub_filtered_list$filtered_df

#-------------------------------------------------------------------------------
# run GLMM and CLR-LMM
#-------------------------------------------------------------------------------
# Project for intercept, Batch_ID for intercept
# asv_all
res_asv_inters_Time <- run_glmmTMB_LMM_plus(asv_raw, asv_clr, level = "OTU",
                                            random_structure = "intercept")
#res_all_inters_Time_50 <- run_glmmTMB_LMM(asv_raw_filtered, asv_clr_filtered)
res_asv_inters_Time_merged <- res_asv_inters_Time$merged

# genus_all 
res_genus_inters_time <- run_glmmTMB_LMM_plus(genus_raw,genus_clr, level = "Genus",
                                              random_structure = "intercept")
res_genus_inters_time_merged <- res_genus_inters_time$merged

# asv_sub
res_asv_inters_Time_sub <- run_glmmTMB_LMM_plus(asv_raw_sub, asv_clr_sub,level = "OTU",
                                                random_structure = "intercept")
res_asv_inters_Time_sub_merged <- res_asv_inters_Time_sub$merged

#genus_sub
res_genus_inters_Time_sub <- run_glmmTMB_LMM_plus(genus_raw_sub,genus_clr_sub, level = "Genus",
                                                  random_structure = "intercept")
res_genus_inters_Time_sub_merged <- res_genus_inters_Time_sub$merged
#-------------------------------------------------------------------------------
# make volcano plot
#-------------------------------------------------------------------------------
#  all_asv
p_volcano_all_asv <- plot_GLMM_LMM_volcano(res_merged = res_all_inters_Time_merged,
  abund_df   = asv_raw_filtered,
  key_col    = "OTU", abund_col  = "RelAbundance"
)
p_volcano_all_asv

#  all_genus
p_volcano_all_genus <- plot_GLMM_LMM_volcano(res_merged = res_genus_inters_time_merged,
                                           abund_df   = genus_raw,
                                           key_col    = "Genus", abund_col  = "RelAbundance")
p_volcano_all_genus

#  dtt1_asv
p_volcano_dtt1_asv <- plot_GLMM_LMM_volcano(res_merged = res_all_inters_Time_dtt1_merged,
                                           abund_df   = asv_raw_dtt1,
                                           key_col    = "OTU", abund_col  = "RelAbundance"
)
p_volcano_dtt1_asv

