#===============================================================================
# run_MaAsLin2_LMM
#-------------------------------------------------------------------------------
# A clean, flexible MaAsLin2 wrapper
# Supports:
#   - Automatic path creation
#   - Optional custom output folder
#   - Automatic logging
#   - R-returned result table
#   - Default save path: "output/MaAsLin2"
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
  
  #-------------------------------------------------------------------------------
  # 1️⃣ Handle and normalize output path
  #-------------------------------------------------------------------------------
  output_dir <- gsub("\\\\", "/", output_dir)               # unify path separator
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # unique subfolder for this run (timestamped)
  subdir <- file.path(output_dir,
                      paste0("MaAsLin2_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(subdir, showWarnings = FALSE, recursive = TRUE)
  
  log_file <- file.path(subdir, "maaslin2.log")
  
  #-------------------------------------------------------------------------------
  # 2️⃣ Run MaAsLin2 model
  #-------------------------------------------------------------------------------
  fit <- Maaslin2(
    input_data = feature_table,
    input_metadata = metadata,
    output = subdir,
    normalization = normalization,
    transform = transform,
    analysis_method = "LM",
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    min_abundance = min_abundance,
    min_prevalence = min_prevalence,
    correction = "BH",
    standardize = FALSE,
    cores = cores
  )
  
  #-------------------------------------------------------------------------------
  # 3️⃣ Load and enrich results
  #-------------------------------------------------------------------------------
  res_path <- file.path(subdir, "all_results.tsv")
  if (!file.exists(res_path)) {
    stop(glue::glue("❌ No result file found at {res_path}. Check MaAsLin2 logs."))
  }
  
  res <- read.delim(res_path) %>%
    filter(metadata %in% fixed_effects | grepl(":", metadata)) %>%
    mutate(
      FDR_lmm = qval,
      est_lmm = coef,
      FC72 = exp(coef * time_diff),
      FC72_low = exp((coef - stderr) * time_diff),
      FC72_high = exp((coef + stderr) * time_diff)
    )
  
  #-------------------------------------------------------------------------------
  # 4️⃣ Optionally keep or delete folder
  #-------------------------------------------------------------------------------
  if (!save_output) {
    unlink(subdir, recursive = TRUE, force = TRUE)
    message("🧹 Temporary MaAsLin2 output deleted after import.")
  } else {
    message(glue::glue("✅ MaAsLin2 output saved at: {subdir}"))
    message(glue::glue("Log file: {log_file}"))
  }
  
  #-------------------------------------------------------------------------------
  # 5️⃣ Return result and model object
  #-------------------------------------------------------------------------------
  return(list(result = res, fit = fit, log_path = log_file, output_dir = subdir))
}
