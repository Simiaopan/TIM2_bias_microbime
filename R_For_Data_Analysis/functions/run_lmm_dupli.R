#-------------------------------------------------------------------------------
# LMM for Within / Between groups + summary table + predictions (w/ CI)
#-------------------------------------------------------------------------------
run_lmm_dupli <- function(df, do_fdr = FALSE,
                          time_seq = seq(0, 72, by = 1),
                          ci_level = 0.95, n_sims = 1000) {
  # dependencies
  suppressPackageStartupMessages({
    library(lme4)
    library(lmerTest)
    library(dplyr)
  })
  has_merTools <- requireNamespace("merTools", quietly = TRUE)
  
  results   <- list()
  all_fixed <- list()
  pred_list <- list()
  
  for (grp in c("Within", "Between")) {
    message(">>> Processing group: ", grp)
    
    subdf <- df[df$Group == grp, , drop = FALSE]
    if (nrow(subdf) == 0) {
      warning("No rows found for group: ", grp)
      next
    }
    
    subdf <- subdf %>%
      dplyr::select(Distance, Time, PairID) %>%
      mutate(PairID = factor(PairID))
    
    message("   Rows: ", nrow(subdf),
            " | Unique PairIDs: ", length(unique(subdf$PairID)))
    
    # --- fit LMM ---
    model <- lmer(Distance ~ Time + (1|PairID),
                  data = subdf,
                  REML = TRUE,
                  control = lmerControl(check.nobs.vs.nRE = "ignore"))
    
    coef_tab <- as.data.frame(coef(summary(model)))
    coef_tab$Effect <- rownames(coef_tab)
    rownames(coef_tab) <- NULL
    coef_tab$Group <- grp
    
    results[[grp]]   <- model
    all_fixed[[grp]] <- coef_tab
    
    # --- predictions: prefer merTools; fallback to Wald approximation ---
    # Only plot overall average trend (fixed effects only), no random effects added
    newdata <- data.frame(Time = time_seq)
    pred_df <- NULL
    
    if (has_merTools) {
      pred_df <- tryCatch({
        preds <- merTools::predictInterval(model,
                                           newdata = newdata,
                                           level = ci_level, n.sims = n_sims,
                                           which = "fixed", include.resid.var = FALSE)
        # Ensure consistent column names: fit/lwr/upr
        nm <- names(preds)
        names(preds)[match(c("median","lower","upper"), nm, nomatch = 0)] <- c("fit","lwr","upr")[match(c("median","lower","upper"), nm, nomatch = 0) != 0]
        if (!all(c("fit","lwr","upr") %in% names(preds))) {
          # Some versions already have fit/lwr/upr; no action needed
        }
        preds$Time  <- newdata$Time
        preds$Group <- grp
        preds[, c("Time","fit","lwr","upr","Group")]
      }, error = function(e) NULL)
    }
    
    # Wald fallback (fixed effects only)
    if (is.null(pred_df)) {
      fe  <- lme4::fixef(model)                # (Intercept), Time
      V   <- as.matrix(stats::vcov(model))     # Variance-covariance matrix of fixed effects
      X   <- stats::model.matrix(~ Time, data = newdata)  # Design matrix (intercept + Time)
      fit <- as.numeric(X %*% fe)
      se  <- sqrt(diag(X %*% V %*% t(X)))
      z   <- qnorm((1 + ci_level)/2)
      lwr <- fit - z * se
      upr <- fit + z * se
      pred_df <- data.frame(Time = newdata$Time, fit = fit, lwr = lwr, upr = upr, Group = grp)
    }
    
    pred_list[[grp]] <- pred_df
  }
  
  # --- combine fixed effects (same logic as before) ---
  fixed_df <- bind_rows(all_fixed)
  
  if ("Std. Error" %in% names(fixed_df)) {
    names(fixed_df)[names(fixed_df) == "Std. Error"] <- "SE"
  }
  if ("t value" %in% names(fixed_df)) {
    names(fixed_df)[names(fixed_df) == "t value"] <- "t"
  }
  if ("Pr(>|t|)" %in% names(fixed_df)) {
    names(fixed_df)[names(fixed_df) == "Pr(>|t|)"] <- "p"
  }
  
  fixed_df <- fixed_df %>%
    dplyr::select(Group, Effect, Estimate, SE, t, p)
  
  if (do_fdr) {
    fixed_df <- fixed_df %>%
      dplyr::group_by(Effect) %>%
      dplyr::mutate(p_adj = p.adjust(p, method = "fdr")) %>%
      dplyr::ungroup()
  }
  
  # --- combine predictions (ensure non-empty and column consistency) ---
  pred_df <- bind_rows(pred_list)
  # Fallback: if empty, return placeholder to avoid plotting errors
  if (nrow(pred_df) == 0) {
    pred_df <- data.frame(Time = numeric(0), fit = numeric(0),
                          lwr = numeric(0), upr = numeric(0),
                          Group = character(0))
  }
  
  return(list(
    models = results,
    fixed  = fixed_df,
    pred   = pred_df
  ))
}
