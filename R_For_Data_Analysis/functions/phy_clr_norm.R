#-------------------------------------------------------------------------------
# Normalize ASV count with clr
#-------------------------------------------------------------------------------
phy_clr_norm <- function(phylo, pseudocount = 0.5) {
  require("phyloseq")
  require("matrixStats")
  otu <- as(otu_table(phylo), "matrix")
  if(!taxa_are_rows(phylo)) otu <- t(otu)
  otu <- otu + pseudocount
  logm <- log(otu)
  gm_mean <- colMeans(logm)
  clr <- sweep(logm, 2, gm_mean, FUN = "-")
  otu_table(phylo) <- otu_table(clr, taxa_are_rows = TRUE)
  return(phylo)
}