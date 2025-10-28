#-------------------------------------------------------------------------------
# test time consistency: LMM + MANOVA + cosine similarity
#-------------------------------------------------------------------------------
test_time_consistency <- function(biplot_df, 
                                  lmm_result,
                                  time_var = "Time",
                                  time_factor_var = "TimePoint",
                                  project_var = "Project",
                                  random_var = "Batch_ID",
                                  pcs_for_manova = 1:5,
                                  pc_prefix = "PC") {
  # ---- Dependencies ----
  suppressPackageStartupMessages({
    library(lme4)
    library(lmerTest)
    library(emmeans)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(car)
    library(pheatmap)
  })
  
  message("------------------------------------------------------")
  message("STEP 1: Extract significant principal components")
  message("------------------------------------------------------")
  
  sig_pcs <- lmm_result$table %>%
    dplyr::filter(Effect %in% c(time_var, paste0(time_var, ":", project_var))) %>%
    dplyr::filter(P_adj < 0.05) %>%
    dplyr::distinct(PC) %>%
    dplyr::pull(PC)
  
  message("Significant time-related PCs: ", paste(sig_pcs, collapse = ", "))
  
  results_slopes <- list()
  
  message("------------------------------------------------------")
  message("STEP 2: Extract time slopes (emtrends) for each Project")
  message("------------------------------------------------------")
  
  for (pc in sig_pcs) {
    message("  → ", pc)
    fmla <- as.formula(paste0(pc, " ~ ", time_var, " * ", project_var, " + (1|", random_var, ")"))
    fit <- lmer(fmla, data = biplot_df, REML = TRUE)
    
    em <- emmeans::emtrends(fit, ~ Project, var = time_var)
    slopes <- summary(em, infer = c(TRUE, TRUE))
    slopes$PC <- pc
    results_slopes[[pc]] <- slopes
  }
  
  slopes_df <- dplyr::bind_rows(results_slopes)
  
  # ---- Compute direction consistency ----
  dir_summary <- slopes_df %>%
    dplyr::group_by(PC) %>%
    dplyr::summarise(
      mean_slope = mean(Time.trend),
      sd_slope = sd(Time.trend),
      prop_positive = mean(Time.trend > 0),
      prop_significant = mean(lower.CL * upper.CL > 0),
      .groups = "drop"
    )
  
  message("------------------------------------------------------")
  message("STEP 3: MANOVA across multiple PCs to test overall time effect")
  message("------------------------------------------------------")
  
  pcs <- paste0(pc_prefix, pcs_for_manova)
  pcs_exist <- pcs[pcs %in% colnames(biplot_df)]
  Y <- as.matrix(biplot_df[, pcs_exist])
  
  fit_mlm <- stats::lm(Y ~ biplot_df[[time_var]] * biplot_df[[project_var]])
  manova_res <- summary(car::Anova(fit_mlm, type = 3), multivariate = TRUE)
  
  message("→ MANOVA results summary:")
  print(manova_res)
  
  #----------------------------------------
  # STEP 4: Multivariate Δ-vector similarity
  #----------------------------------------
  message("------------------------------------------------------")
  message("STEP 4: Multivariate Δ-vector similarity (trajectory direction similarity)")
  message("------------------------------------------------------")
  
  pcs_sel <- pcs_exist
  delta <- biplot_df %>%
    dplyr::group_by(!!rlang::sym(project_var), !!rlang::sym(time_factor_var)) %>% 
    dplyr::summarise(dplyr::across(dplyr::all_of(pcs_sel), mean), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = !!rlang::sym(time_factor_var),
                       values_from = dplyr::all_of(pcs_sel))
  
  time_levels <- sort(unique(biplot_df[[time_factor_var]]))
  start_tp <- time_levels[1]
  end_tp <- time_levels[length(time_levels)]
  
  for (pc in pcs_sel) {
    delta[[paste0("d_", pc)]] <- delta[[paste0(pc, "_", end_tp)]] - delta[[paste0(pc, "_", start_tp)]]
  }
  
  delta_mat <- as.matrix(delta %>% dplyr::select(dplyr::starts_with("d_")))
  
  cos_sim <- function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
  sim_mat <- outer(1:nrow(delta_mat), 1:nrow(delta_mat),
                   Vectorize(function(i, j) cos_sim(delta_mat[i, ], delta_mat[j, ])))
  rownames(sim_mat) <- delta[[project_var]]
  colnames(sim_mat) <- delta[[project_var]]
  
  p <- pheatmap::pheatmap(
    sim_mat,
    main = "Cosine similarity of time-trajectory vectors across projects",
    fontsize = 10,
    fontsize_row = 9,
    fontsize_col = 9,
    treeheight_row = 30,
    treeheight_col = 30
  )
  
  message("------------------------------------------------------")
  message("Analysis completed ✅")
  message("------------------------------------------------------")
  
  return(list(
    slopes = slopes_df,
    dir_summary = dir_summary,
    manova = manova_res,
    cos_sim = sim_mat,
    heatmap = p
  ))
}
