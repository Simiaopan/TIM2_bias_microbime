#-------------------------------------------------------------------------------
# LMM for each selected-rpca-axis
#-------------------------------------------------------------------------------
run_rpca_lmm <- function(biplot_df, 
                         pc_prefix = "PC", 
                         n_pcs = 15, 
                         time_var = "Time",         # numeric
                         project_var = "Project", 
                         random_var = "Batch_ID",
                         fdr_method = "BH") {
  
  # --- load + explicit namespace binding to prevent pollution ---
  suppressPackageStartupMessages({
    library(lme4)
    library(lmerTest)
    library(dplyr)
    library(tidyr)
  })
  
  # explicit namespace binding
  select <- dplyr::select
  filter <- dplyr::filter
  mutate <- dplyr::mutate
  pull   <- dplyr::pull
  group_by <- dplyr::group_by
  ungroup <- dplyr::ungroup
  summarise <- dplyr::summarise
  bind_rows <- dplyr::bind_rows
  
  message(">>> Running Linear Mixed Models (LMM) for each selected RPCA axis ...")
  
  #-------------------------------------------------------------
  # select axis
  #-------------------------------------------------------------
  pc_cols <- paste0(pc_prefix, 1:n_pcs)
  pcs_exist <- pc_cols[pc_cols %in% colnames(biplot_df)]
  if (length(pcs_exist) == 0) stop("No PC columns found in dataset!")
  
  results <- list()
  
  message("------------------------------------------------------")
  message("Running LMMs on ", length(pcs_exist), " PCs ...")
  message("Using dplyr::select from version: ", as.character(utils::packageVersion("dplyr")))
  message("------------------------------------------------------")
  
  #-------------------------------------------------------------
  # modeling for every selected axis
  #-------------------------------------------------------------
  for (pc in pcs_exist) {
    formula_str <- paste0(pc, " ~ ", time_var, " * ", project_var, " + (1 | ", random_var, ")")
    message(">>> Fitting model for ", pc, ": ", formula_str)
    
    model <- lmer(as.formula(formula_str), data = biplot_df, REML = TRUE)
    anova_res <- anova(model, type = 3)  # Type III ANOVA
    
    # --- safe extraction of random effect variance ---
    vc_df <- as.data.frame(VarCorr(model))
    col_vcov <- intersect(names(vc_df), c("vcov", "Variance", "vcov.", "vcovariance"))
    if (length(col_vcov) == 0) {
      warning(pc, ": cannot find a variance column in VarCorr(model), using NA.")
      rand_var <- NA_real_
    } else {
      rand_var <- vc_df %>% filter(grp == random_var)
      rand_var <- if (nrow(rand_var) == 0) NA_real_ else rand_var[[col_vcov[1]]][1]
    }
    
    fixed_pvals <- data.frame(
      PC = pc,
      Effect = rownames(anova_res),
      Fvalue = anova_res$`F value`,
      Pvalue = anova_res$`Pr(>F)`,
      RandomVar = rand_var
    )
    results[[pc]] <- fixed_pvals
  }
  
  #-------------------------------------------------------------
  # merge + adjust
  #-------------------------------------------------------------
  summary_df <- bind_rows(results)
  
  summary_df <- summary_df %>%
    group_by(Effect) %>%
    mutate(P_adj = p.adjust(Pvalue, method = fdr_method),
           Signif = P_adj < 0.05) %>%
    ungroup()
  
  #-------------------------------------------------------------
  # summarize the result
  #-------------------------------------------------------------
  summary_stats <- summary_df %>%
    group_by(Effect) %>%
    summarise(
      n_tests = n(),
      n_sig = sum(Signif),
      prop_sig = round(mean(Signif), 3),
      min_P = signif(min(P_adj), 3),
      max_P = signif(max(P_adj), 3),
      .groups = "drop"
    )
  
  message("------------------------------------------------------")
  message("FDR correction done (method = ", fdr_method, ")")
  message("Summary of significant PCs per effect:")
  print(summary_stats)
  message("------------------------------------------------------")
  
  return(list(
    table = summary_df,
    summary = summary_stats
  ))
}
