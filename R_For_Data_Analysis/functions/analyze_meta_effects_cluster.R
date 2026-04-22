##--------------------------------------------------------------------
## Function:
##   - Perform LASSO multinomial logistic regression to assess 
##     the independent contribution of meta variables to cluster structure
##   - Evaluate overall explanatory power of TimePoint / Project / RUN
##   - Output and visualize key results
## -------------------------------------------------------------------

analyze_meta_effects_cluster <- function(RPCA_biplot_table, 
                                         lambda_ratio = 1e-3, 
                                         min_cluster_size = 8,
                                         seed = 20021) {
  ## -------------------------
  ## Step 0. Load dependencies
  ## -------------------------
  suppressPackageStartupMessages({
    library(dplyr)
    library(glmnet)
    library(ggplot2)
    library(pheatmap)
    library(tidyr)
  })
  
  ## -------------------------
  ## Step 1. Data preparation and filtering
  ## -------------------------
  cat("=== Step 1. Data filtering ===\n")
  df <- RPCA_biplot_table %>%
    filter(Cluster != "Noise") %>%
    mutate(
      Cluster = factor(Cluster),
      TimePoint = factor(TimePoint),
      Project = factor(Project),
      RUN = factor(RUN)
    )
  
  cat("Original Cluster distribution:\n")
  print(table(df$Cluster))
  
  keep_clusters <- names(which(table(df$Cluster) >= min_cluster_size))
  df_sub <- df %>%
    filter(Cluster %in% keep_clusters) %>%
    mutate(Cluster = droplevels(Cluster))
  
  cat("Remaining number of clusters:", length(unique(df_sub$Cluster)), "\n")
  cat("Remaining number of samples:", nrow(df_sub), "\n\n")
  
  ## -------------------------
  ## Step 2. Build design matrix
  ## -------------------------
  x <- model.matrix(Cluster ~ TimePoint + Project + RUN, data = df_sub)[, -1]
  y <- df_sub$Cluster
  cat("=== Step 2. Design matrix built ===\n")
  print(table(y))
  
  ## -------------------------
  ## Step 3. Run LASSO multinomial logistic regression (CV)
  ## -------------------------
  cat("\n=== Step 3. LASSO cross-validation ===\n")
  set.seed(seed)
  cvfit <- cv.glmnet(
    x, y,
    family = "multinomial",
    type.multinomial = "grouped",
    alpha = 1,
    nfolds = 10,
    maxit = 200000,
    lambda.min.ratio = lambda_ratio
  )
  
  plot(cvfit)
  title("Cross-validation curve for multinomial LASSO", line = 2.5)
  cvfit_plot <- recordPlot()
  
  best_lambda <- cvfit$lambda.1se
  cat("Selected λ (1se):", best_lambda, "\n")
  
  ## -------------------------
  ## Step 4. Refit the model
  ## -------------------------
  model_lasso <- glmnet(
    x, y,
    family = "multinomial",
    type.multinomial = "grouped",
    alpha = 1,
    lambda = best_lambda
  )
  
  ## -------------------------
  ## Step 5. Extract coefficient matrix
  ## -------------------------
  coef_list <- coef(model_lasso)
  coef_df <- do.call(rbind, lapply(names(coef_list), function(cl) {
    co <- as.matrix(coef_list[[cl]])
    data.frame(Cluster = cl,
               Variable = rownames(co),
               Coef = co[, 1],
               stringsAsFactors = FALSE)
  })) %>%
    filter(Variable != "(Intercept)")
  
  ## -------------------------
  ## Step 6. Overall variable importance
  ## -------------------------
  importance_df <- coef_df %>%
    group_by(Variable) %>%
    summarise(TotalEffect = sum(abs(Coef))) %>%
    arrange(desc(TotalEffect))
  
  cat("\nTop variables explaining Cluster structure:\n")
  print(head(importance_df, 10))
  
  ## -------------------------
  ## Step 7. Visualize overall importance
  ## -------------------------
  importance_plot <- ggplot(importance_df, aes(x = reorder(Variable, TotalEffect), y = TotalEffect)) +
    geom_col(fill = "#2b83ba") +
    coord_flip() +
    theme_minimal(base_size = 13) +
    labs(x = "Variable", y = "Overall |β|",
         title = "Meta-variable importance in explaining Cluster structure")
  print(importance_plot)
  
  ## -------------------------
  ## Step 8. Coefficient heatmap
  ## -------------------------
  coef_matrix <- coef_df %>%
    pivot_wider(names_from = Cluster, values_from = Coef, values_fill = 0) %>%
    as.data.frame()
  rownames(coef_matrix) <- coef_matrix$Variable
  coef_matrix <- coef_matrix[, -1]
  
  pheatmap(as.matrix(coef_matrix),
           color = colorRampPalette(c("#4A6FA5", "white", "#E56B6F"))(100),
           main = "LASSO coefficients (Variable vs Cluster)",
           fontsize_row = 7, fontsize_col = 8, fontsize = 6,
           treeheight_row = 25, treeheight_col = 25)
  coef_heatmap <- recordPlot()   # ✅ Save
  
  ## -------------------------
  ## Step 9. Independent contributions: Project vs RUN
  ## -------------------------
  project_run_effect <- coef_df %>%
    filter(grepl("^Project|^RUN", Variable)) %>%
    group_by(Variable) %>%
    summarise(TotalEffect = sum(abs(Coef))) %>%
    arrange(desc(TotalEffect))
  
  cat("\nIndependent contributions (Project vs RUN):\n")
  print(project_run_effect)
  
  project_run_plot <- ggplot(project_run_effect, aes(x = reorder(Variable, TotalEffect), y = TotalEffect)) +
    geom_col(fill = "#d7191c") +
    coord_flip() +
    theme_minimal(base_size = 13) +
    labs(x = "Variable", y = "Total |β|",
         title = "Independent contributions of Project and RUN")
  print(project_run_plot)
  
  ## -------------------------
  ## Step 10–12. Overall meta-level contribution
  ## -------------------------
  meta_importance <- coef_df %>%
    mutate(Meta = case_when(
      grepl("^Project", Variable) ~ "Project",
      grepl("^RUN", Variable) ~ "RUN",
      grepl("^TimePoint", Variable) ~ "TimePoint",
      TRUE ~ "Other"
    )) %>%
    group_by(Meta) %>%
    summarise(TotalEffect = sum(abs(Coef))) %>%
    mutate(Percent = round(TotalEffect / sum(TotalEffect) * 100, 2)) %>%
    arrange(desc(TotalEffect))
  
  cat("\n=== Independent contribution of each meta variable ===\n")
  print(meta_importance)
  
  meta_plot <- ggplot(meta_importance, aes(x = reorder(Meta, TotalEffect), y = TotalEffect, fill = Meta)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("Project" = "#E69F00", "RUN" = "#56B4E9", "TimePoint" = "#009E73"))+
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), legend.position = "none") +
    labs(x = "Meta variable", y = "Total |β|", title = "Relative independent contributions of meta variables")
  
  print(meta_plot)
  
  
  run_residual_ratio <- meta_importance %>%
    filter(Meta == "RUN") %>%
    pull(Percent)
  
  cat("\nRUN residual explanatory power:", run_residual_ratio, "%\n")
  
  ## -------------------------
  ## Step 13. Save all plot objects
  ## -------------------------
  plots <- list(
    cvfit_plot = cvfit_plot,
    importance_plot = importance_plot,
    coef_heatmap = coef_heatmap,
    project_run_plot = project_run_plot,
    meta_plot = meta_plot
  )
  
  ## -------------------------
  ## Step 14. Compile and return all results
  ## -------------------------
  result <- list(
    filtered_data = df_sub,
    best_lambda = best_lambda,
    cvfit = cvfit,
    model_lasso = model_lasso,
    coef_df = coef_df,
    importance_df = importance_df,
    project_run_effect = project_run_effect,
    meta_importance = meta_importance,
    run_residual_ratio = run_residual_ratio,
    plots = plots
  )
  
  cat("\n=== Analysis complete ===\n")
  return(result)
}
