#-------------------------------------------------------------------------------
# Set environment 
#-------------------------------------------------------------------------------
load("rdata/prm.Rdata")
#load("rdata/data.Rdata")
load("rdata/data.ls_90.Rdata")

seed = prm.ls$general$seed_all
set.seed(seed)

suppressPackageStartupMessages({
  library(phyloseq)
  library(qiime2R)
  library(dplyr)
})

#-------------------------------------------------------------------------------
# Load custom functions
#-------------------------------------------------------------------------------
source("functions/export_for_deicode.R")
source("functions/merge_rpca_meta.R")

#-------------------------------------------------------------------------------
# User options
#-------------------------------------------------------------------------------
save_file <- "rdata/RPCA_filt_t1.Rdata" # change it when process new datasets
input.dir <- prm.ls$RPCA$input_dir
output.dir <- prm.ls$RPCA$output_dir # change it when process new datasets
dataset.names <- prm.ls$RPCA$dataset_names
data.ls <- data.ls_90 # change it when process new datasets
#-------------------------------------------------------------------------------
# extract data
#-------------------------------------------------------------------------------
#ps.list <- data.ls$Ps$ASV$RAW[dataset.names]
ps.list <- data.ls$RAW[dataset.names]
#-------------------------------------------------------------------------------
# Prepare data for RPCA (deicode)
#-------------------------------------------------------------------------------
export_for_deicode(ps.list, outdir = input.dir)

#-------------------------------------------------------------------------------
# Run RPCA (deicode) in bash
#-------------------------------------------------------------------------------
print("Run RPCA need to be done in bash with qiime2, script: run_deicode_rcpca.sh")

#-------------------------------------------------------------------------------
# merge rPCA data with meta
#-------------------------------------------------------------------------------
# RPCA <- list()
# RPCA <- merge_rpca_meta(dataset.names,indir = "data/rpca_output",
#                         metadir = "data/rpca_input")
# 
# RPCA[["dtt1"]]$proportion

RPCA_filt <- list()
RPCA_filt <- merge_rpca_meta(dataset.names,indir = output.dir,
                        metadir = input.dir)

RPCA_filt[["t1"]]$proportion


#----------------- save data ----------------------
# save(list = c("RPCA"), file = "rdata/RPCA.Rdata")

save(list = c("RPCA_filt"), file = save_file)

rm(list = ls())
gc()
