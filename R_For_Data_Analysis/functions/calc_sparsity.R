#-------------------------------------------------------------------------------
# calculate sparcity (on phyloseq)
#-------------------------------------------------------------------------------

calc_sparsity <- function(phylo) {
  otu <- as(otu_table(phylo), "matrix")
  if(!taxa_are_rows(phylo)) otu <- t(otu)
  
  # overall
  zero_fraction <- sum(otu == 0) / length(otu)
  
  # per sample
  sample_sparsity <- rowMeans(otu == 0)
  
  # per species
  taxa_sparsity <- colMeans(otu == 0)
  
  return(list(
    zero_fraction = zero_fraction,
    sample_sparsity = sample_sparsity,
    taxa_sparsity = taxa_sparsity
  ))
}