filter_phyloseq <- function(data_list, keep_asv) {
  library(phyloseq)
  
  # Step 1 — Apply filtering recursively
  data_filtered <- lapply(data_list, function(sublist) {
    if (inherits(sublist, "phyloseq")) {
      # Case 1: Single phyloseq object
      prune_taxa(taxa_names(sublist) %in% keep_asv, sublist)
    } else if (is.list(sublist)) {
      # Case 2: Nested list of phyloseq objects
      lapply(sublist, function(ps) {
        if (inherits(ps, "phyloseq")) {
          prune_taxa(taxa_names(ps) %in% keep_asv, ps)
        } else {
          ps
        }
      })
    } else {
      sublist
    }
  })
  
  # Step 2 — Print simple check
  cat("---------------------------------------------------\n")
  cat("Check results after filtering:\n")
  
  # Safely check if elements exist
  if ("RAW" %in% names(data_filtered)) {
    raw_counts <- sapply(data_filtered$RAW, ntaxa)
    cat(" RAW:\n")
    print(raw_counts)
  }
  if ("CLR" %in% names(data_filtered)) {
    clr_counts <- sapply(data_filtered$CLR, ntaxa)
    cat(" CLR:\n")
    print(clr_counts)
  }
  
  cat(sprintf("\nTotal ASVs retained: %d\n", length(keep_asv)))
  cat("---------------------------------------------------\n")
  
  return(data_filtered)
}
