#------------------------------------------------------------------------------
# MaAsLin2 diagnostic and labeling function
#------------------------------------------------------------------------------
diagnose_MaAsLin2_fit <- function(
    fit,
    merged_df,
    feature_col = "feature",
    cor_threshold = 0.3,
    hetero_p = 0.05,
    rmse_q = 0.95,
    sample_size = 20000
) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(rlang)
  })
  
  message("▶ Step 1. Extract fitted/residual matrices and feature list")
  F <- fit$fitted
  R <- fit$residuals
  feat_fit <- unique(na.omit(fit$results$feature))
  feat_merged <- unique(na.omit(merged_df[[feature_col]]))
  
  if (is.null(colnames(F)) && length(feat_fit) == ncol(F)) {
    colnames(F) <- feat_fit
    colnames(R) <- feat_fit
    message("Automatically assigned feature names to F/R from fit$results$feature")
  }
  if (nrow(F) == length(feat_fit)) {
    F <- t(F); R <- t(R)
  }
  
  # Auto-align features between fitted and merged data
  common_features <- intersect(na.omit(colnames(F)), na.omit(feat_merged))
  if (length(common_features) < 5) {
    stop(sprintf("❌ Only %d shared features found; please check naming or filtering.", length(common_features)))
  } else if (length(common_features) < length(feat_merged)) {
    warning(sprintf("⚠️ %d shared features found, fewer than merged_df (%d); automatically aligning by intersection.",
                    length(common_features), length(feat_merged)))
  }
  
  F <- F[, common_features, drop = FALSE]
  R <- R[, common_features, drop = FALSE]
  merged_df <- merged_df %>% filter(.data[[feature_col]] %in% common_features)
  message("✅ Feature alignment complete: ", length(common_features), " features")
  
  # Step 2: Compute basic diagnostic metrics
  message("Step 2. Compute basic diagnostic metrics")
  mean_fit <- colMeans(F, na.rm = TRUE)
  sd_resid <- apply(R, 2, sd, na.rm = TRUE)
  rmse     <- sqrt(colMeans(R^2, na.rm = TRUE))
  
  safe_cor <- function(x, y) {
    res <- try(suppressWarnings(cor(abs(x), y, use = "pairwise")), silent = TRUE)
    if (inherits(res, "try-error") || is.na(res)) return(NA_real_) else return(res)
  }
  cor_abs <- sapply(seq_len(ncol(R)), function(j) safe_cor(R[, j], F[, j]))
  
  # Step 3: Heteroscedasticity test (Breusch–Pagan test + sampling)
  message("Step 3. Heteroscedasticity test (BP test + sampling)")
  max_feats <- min(300, ncol(R))
  sample_feats <- sample(seq_len(ncol(R)), max_feats)
  
  p_hetero <- rep(NA_real_, ncol(R))
  for (j in sample_feats) {
    res_j <- as.numeric(R[, j])
    fit_j <- as.numeric(F[, j])
    if (sum(is.finite(res_j + fit_j)) > 20 && sd(fit_j, na.rm = TRUE) > 0) {
      mod <- try(lm(res_j ~ fit_j), silent = TRUE)
      if (!inherits(mod, "try-error")) {
        bp <- try(lmtest::bptest(mod, varformula = ~ fit_j), silent = TRUE)
        if (!inherits(bp, "try-error")) p_hetero[j] <- bp$p.value
      }
    }
  }
  q_hetero <- p.adjust(p_hetero, "BH")
  
  # Step 4: Summarize diagnostic results
  message("Step 4. Summarize diagnostic result table")
  diag_df <- data.frame(
    feature = common_features,
    mean_fit, sd_resid, rmse, cor_abs,
    p_hetero, q_hetero
  ) %>%
    mutate(
      flag_rmse   = rmse > quantile(rmse, rmse_q, na.rm = TRUE),
      flag_cor    = abs(cor_abs) > cor_threshold,
      flag_hetero = !is.na(q_hetero) & q_hetero < hetero_p,
      suspicious  = flag_rmse + flag_cor + flag_hetero,
      suspect_flag = suspicious >= 2,
      hetero_flag  = flag_hetero & flag_cor
    ) %>%
    arrange(desc(suspicious), desc(rmse))
  
  # Step 5: Visualization
  message("Step 5. Plot diagnostics")
  g1 <- ggplot(diag_df, aes(mean_fit, sd_resid)) +
    geom_point(alpha = .6, size = 1) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "red") +
    labs(title = "Per-feature residual SD vs mean fitted",
         x = "Mean fitted value", y = "Residual SD") +
    theme_minimal()
  
  pooled <- data.frame(fitted = as.vector(F), residual = as.vector(R))
  n_sample <- min(sample_size, nrow(pooled))
  set.seed(42)
  pooled_sub <- pooled[sample(seq_len(nrow(pooled)), n_sample), ]
  
  g2 <- ggplot(pooled_sub, aes(fitted, residual)) +
    geom_bin2d(bins = 60) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = "red") +
    labs(title = "Pooled residuals vs fitted (sampled)",
         subtitle = paste("Sampled", n_sample, "points"),
         x = "Fitted value", y = "Residual") +
    theme_minimal()
  
  # Step 6: Merge model_quality back to merged_df — fix NSE issues
  message("Step 6. Merge into merged_df and mark model_quality")
  merged_out <- merged_df %>%
    left_join(diag_df %>% select(feature, suspect_flag),
              by = setNames("feature", feature_col)) %>%
    mutate(model_quality = case_when(
      is.na(suspect_flag) ~ "Unknown",
      suspect_flag        ~ "Suspect",
      TRUE                ~ "Reliable"
    ))
  
  # Step 7: Output summary statistics
  message("Step 7. Summary statistics")
  tab_quality <- table(merged_out$model_quality)
  print(tab_quality)
  
  invisible(list(
    diag_df = diag_df,
    merged_out = merged_out,
    g1 = g1,
    g2 = g2,
    suspect_features = diag_df %>% filter(suspect_flag)
  ))
}

