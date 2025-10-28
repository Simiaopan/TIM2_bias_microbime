#-------------------------------------------------------------------------------
# Split data based on Duplication, TimePattern, Testproduct, & TimePoint
#-------------------------------------------------------------------------------
make_splits <- function(phy) {
  sd <- sample_data(phy)
  need_cols <- c("Duplication", "TimePattern", "Testproduct", "TimePoint")
  miss <- setdiff(need_cols, colnames(sd))
  if (length(miss) > 0) {
    warning(paste0("Metadata is missing columns: ", paste(miss, collapse = ", "),
                   ". Splitting skipped (returning empty list)."))
    return(list())
  }
  
  pick <- function(cond) sample_names(phy)[which(cond & !is.na(cond))]
  
  lst <- list(
    dt      = prune_samples(pick(sd$Duplication == TRUE), phy),
    dtt1    = prune_samples(pick(sd$Duplication == TRUE & sd$TimePattern == 1), phy),
    dtt1ts  = prune_samples(pick(sd$Duplication == TRUE & sd$TimePattern == 1 & sd$Testproduct == "SIEM"), phy),
    dtt1tns = prune_samples(pick(sd$Duplication == TRUE & sd$TimePattern == 1 & sd$Testproduct != "SIEM"), phy),
    ts      = prune_samples(pick(sd$Testproduct == "SIEM"), phy),
    tns     = prune_samples(pick(sd$Testproduct != "SIEM"), phy),
    t1      = prune_samples(pick(sd$TimePattern == 1), phy)
  )
  
  #-------------------------------------------------------------------------------
  # Add time-point subsets for both dtt1 and t1
  #-------------------------------------------------------------------------------
  for (tpoint in c("T0", "T24", "T48", "T72")) {
    # --- dtt1 by timepoint
    nm_dtt1 <- paste0("dtt1_", tpoint)
    lst[[nm_dtt1]] <- prune_samples(
      pick(sd$Duplication == TRUE & sd$TimePattern == 1 & sd$TimePoint == tpoint), phy
    )
    
    # --- t1 by timepoint
    nm_t1 <- paste0("t1_", tpoint)
    lst[[nm_t1]] <- prune_samples(
      pick(sd$TimePattern == 1 & sd$TimePoint == tpoint), phy
    )
  }
  
  message("✅ Splitting complete: ",
          length(lst), " subsets generated (including time-split for dtt1 and t1).")
  
  return(lst)
}
