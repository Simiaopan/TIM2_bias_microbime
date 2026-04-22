#===============================================================================
# run_MaAsLin2_LMM
#===============================================================================

run_MaAsLin2_LMM <- function(feature_table, metadata,
                             output_dir = "output/MaAsLin2",
                             fixed_effects = c("Time_c"),
                             random_effects = "Batch_ID",
                             normalization = "TSS",
                             transform = "LOG",
                             min_abundance = 1e-4,
                             min_prevalence = 0.1,
                             cores = 7,
                             time_diff = 72,
                             save_output = TRUE) {
  
  suppressPackageStartupMessages(library(Maaslin2))
  suppressPackageStartupMessages(library(dplyr))
  

