#-------------------------------------------------------------------------------
# Duplication stability test (sample vs centroid mode)
#-------------------------------------------------------------------------------
duplication_stability_test <- function(coords_all,
                                  meta,
                                  treat_col   = "Treat_ID",
                                  batch_col   = "Batch_ID",
                                  tp_col      = "TimePoint",
                                  mode        = c("sample","centroid"),
                                  dist_method = "euclidean",
                                  n_iter      = 300) {
  mode <- match.arg(mode)
  
  # 1. Align coordinates with metadata
  stopifnot(all(rownames(coords_all) %in% rownames(meta)))
  meta <- meta[rownames(coords_all), , drop=FALSE]
  
  # distance matrix
  Dmat <- as.matrix(dist(coords_all, method=dist_method))
  ut   <- upper.tri(Dmat)
  
  treat <- meta[[treat_col]]
  batch <- meta[[batch_col]]
  tp    <- meta[[tp_col]]
  
  # 2. Identify real duplication pairs
  if (mode == "sample") {
    # same Treat_ID, same TimePoint, different Batch_ID
    mask <- ut & outer(treat,treat,"==") & outer(tp,tp,"==") & !outer(batch,batch,"==")
    d_within <- Dmat[mask]
    
  } else if (mode == "centroid") {
    # calculate centroids for (Treat_ID, Batch_ID, TimePoint)
    df_coords <- as.data.frame(coords_all)
    centroids <- aggregate(df_coords,
                           by = list(Treat=meta[[treat_col]],
                                     Batch=meta[[batch_col]],
                                     TP=meta[[tp_col]]),
                           FUN = mean)
    rownames(centroids) <- paste(centroids$Treat, centroids$Batch, centroids$TP, sep="_")
    coords_c <- centroids[, -(1:3), drop=FALSE]
    
    # centroid distance matrix
    Dmat_c <- as.matrix(dist(coords_c, method=dist_method))
    ut_c   <- upper.tri(Dmat_c)
    
    treat_c <- centroids$Treat
    batch_c <- centroids$Batch
    tp_c    <- centroids$TP
    
    mask <- ut_c & outer(treat_c,treat_c,"==") & outer(tp_c,tp_c,"==") & !outer(batch_c,batch_c,"==")
    d_within <- Dmat_c[mask]
  }
  
  if (length(d_within) == 0) {
    stop("No valid duplication pairs found under given mode.")
  }
  
  # 3. Permutation test
  T_obs <- mean(d_within)
  T_perm <- numeric(n_iter)
  
  for (b in seq_len(n_iter)) {
    d_fake_all <- c()
    
    if (mode == "sample") {
      # shuffle Treat_ID within each TimePoint
      for (lev in unique(tp)) {
        idx_tp <- which(tp == lev)
        if (length(idx_tp) < 2) next
        tr_shuf <- sample(treat[idx_tp])
        pairs <- combn(seq_along(idx_tp), 2)
        valid <- tr_shuf[pairs[1,]] != tr_shuf[pairs[2,]]
        if (any(valid)) {
          d_fake <- Dmat[idx_tp[pairs[1, valid]], idx_tp[pairs[2, valid]]]
          d_fake_all <- c(d_fake_all, d_fake)
        }
      }
    } else {
      # shuffle Treat_ID within each TimePoint (centroid version)
      for (lev in unique(tp_c)) {
        idx_tp <- which(tp_c == lev)
        if (length(idx_tp) < 2) next
        tr_shuf <- sample(treat_c[idx_tp])
        pairs <- combn(seq_along(idx_tp), 2)
        valid <- tr_shuf[pairs[1,]] != tr_shuf[pairs[2,]]
        if (any(valid)) {
          d_fake <- Dmat_c[idx_tp[pairs[1, valid]], idx_tp[pairs[2, valid]]]
          d_fake_all <- c(d_fake_all, d_fake)
        }
      }
    }
    
    if (length(d_fake_all) > 0) {
      T_perm[b] <- mean(d_fake_all)
    } else {
      T_perm[b] <- NA
    }
  }
  T_perm <- na.omit(T_perm)
  
  # 4. Compute permutation p-value
  pval <- (sum(T_perm <= T_obs) + 1) / (length(T_perm) + 1)
  
  # 5. Return results
  out <- list(
    mode     = mode,
    n_real   = length(d_within),
    T_obs    = T_obs,
    T_perm   = T_perm,
    pval     = pval
  )
  class(out) <- "dup_stability"
  return(out)
}
