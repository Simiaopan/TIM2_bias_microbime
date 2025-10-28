#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------

load("rdata/prm.Rdata")
load("rdata/data.Rdata")

seed = prm.ls$general$seed_all
set.seed(seed)

suppressPackageStartupMessages({
  library(phyloseq)
  library(qiime2R)
  library(tidyverse)
})
#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/add_relative_abundance.R")
source("functions/filter_asv_time_project.R")
source("functions/check_abundance_distribution.R")
source("functions/filter_core_asv.R")
source("functions/plot_asv_distribution_T.R")
source("functions/filter_phyloseq.R")
#-------------------------------------------------------------------------------
# User options
#-------------------------------------------------------------------------------
data_set_vec <- c("_all","t1","dtt1","dtt1ts", "dtt1_T0","dtt1_T24","dtt1_T48",
                  "dtt1_T72","t1_T0","t1_T24", "t1_T48", "t1_T72")
save_file <- "rdata/filt_data.Rdata"
#-------------------------------------------------------------------------------
# extract data and add relative abundance
#-------------------------------------------------------------------------------
# RAW data (ASV) for filtering
ps_raw_asv <- data.ls$Ps$ASV$RAW$`_all`
asv_raw <- ps_raw_asv %>%
  psmelt(.) %>%
  select(.,OTU, Abundance, Sample, Project, RUN, Genus,Species, 
         Time, TimePoint, TimePattern, Batch_ID, Treat_ID)
asv_raw <- add_relative_abundance(asv_raw)

# create new data.ls
data.ls_new1 <- list()
data.ls_new1$RAW <- data.ls$Ps$ASV$RAW[data_set_vec]
data.ls_new2 <- list()
data.ls_new2$CLR <- data.ls$Ps$ASV$CLR[data_set_vec]
data.ls_new <- c(data.ls_new1, data.ls_new2)
#-------------------------------------------------------------------------------
# filter asv by Time point, project
# check distrubution before and after filtering
#-------------------------------------------------------------------------------
check_abundance_distribution(asv_raw, 
                             level = "OTU", value_col = "RelAbundance")

filtered_list <- filter_asv_time_project(df = asv_raw, prev_cutoff_time = 0.25, 
                mean_cutoff = 1e-4, min_projects = 3, min_timepoints_per_project = 2)
asv_raw_filtered <- filtered_list$filtered_df

check_abundance_distribution(asv_raw_filtered, 
                             level = "OTU", value_col = "RelAbundance")

#-------------------------------------------------------------------------------
# filter phyloseq list based on the last step
#-------------------------------------------------------------------------------
keep_asv_robust <- unique(asv_raw_filtered$OTU)


data.ls_filt <- lapply(data.ls_new, function(sublist) {
  if (inherits(sublist, "phyloseq")) {
    prune_taxa(taxa_names(sublist) %in% keep_asv_robust, sublist)
  } else if (is.list(sublist)) {
    lapply(sublist, function(ps) {
      prune_taxa(taxa_names(ps) %in% keep_asv, ps)
    })
  } else {
    sublist
  }
})

# simplely check
sapply(data.ls_filt$RAW, ntaxa)
sapply(data.ls_filt$CLR, ntaxa)
length(keep_asv)

# ------------- core_asv: 80% 50% 25%----------------------------------
# get asv table
core_asv_90_list <- filter_core_asv(asv_raw, 0.9, 0, 3)
core_asv_80_list <- filter_core_asv(asv_raw, 0.8, 0, 3)
core_asv_50_list <- filter_core_asv(asv_raw, 0.5, 0, 3)
core_asv_25_list <- filter_core_asv(asv_raw, 0.25,0, 3)

check_abundance_distribution(core_asv_80_list$filtered_table, 
                             level = "OTU", value_col = "RelAbundance")

# plot changing with time
dt_90 <- plot_core_asv_T(core_asv_90_list$filtered_table, subject_col = "Batch_ID")
dt_80 <- plot_core_asv_T(core_asv_80_list$filtered_table, subject_col = "Batch_ID")
dt_50 <- plot_core_asv_T(core_asv_50_list$filtered_table, subject_col = "Batch_ID")
dt_25 <- plot_core_asv_T(core_asv_25_list$filtered_table, subject_col = "Batch_ID")
p90 <- dt_90$plot
p80 <- dt_80$plot
p50 <- dt_50$plot
p25 <- dt_25$plot

source("functions/combine_core_asv_plots.R")
p_combined <- combine_core_asv_plots(p90, p80, p50, p25)
p_combined
# filter ps
data.ls_80 <- filter_phyloseq(data.ls_new, core_asv_80_list$keep_asv)
data.ls_90 <- filter_phyloseq(data.ls_new, core_asv_90_list$keep_asv)

save(list = c("filtered_list", "data.ls_filt"), file = save_file)
save(list = c("data.ls_80"), file = "rdata/data.ls_80.Rdata")
save(list = c("data.ls_90"), file = "rdata/data.ls_90.Rdata")


rm(list = ls())
gc()
