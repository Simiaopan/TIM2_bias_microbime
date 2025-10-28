#-------------------------------------------------------------------------------
# merge rPCA data with meta
#-------------------------------------------------------------------------------
merge_rpca_meta <- function(dataset_names,
                            indir = "rpca_results",
                            metadir = "exported") {
  library(qiime2R)
  library(dplyr)
  
  results <- list()
  
  for (dataset_name in dataset_names) {
    # import RPCA results
    biplot <- read_qza(file.path(indir, paste0(dataset_name, "_rpca_biplot.qza")))$data
    distance <- read_qza(file.path(indir, paste0(dataset_name, "_rpca_distance.qza")))$data
    if (is.matrix(distance)) distance <- as.dist(distance)
    
    # extract sample coords from Vectors
    sample_coords <- biplot$Vectors

    if (!"SampleID" %in% colnames(sample_coords)) {
      stop("No SampleID column found in Vectors for dataset: ", dataset_name)
    }
    
    # import meta
    meta <- readRDS(file.path(metadir, paste0(dataset_name, "_meta.rds")))
    meta$SampleID <- rownames(meta)
    
    # merge
    merged <- sample_coords %>%
      inner_join(meta, by = "SampleID")
    
    results[[dataset_name]] <- list(
      biplot = merged,
      distance = distance,
      eigvals = biplot$Eigvals,
      proportion = biplot$ProportionExplained,
      species = biplot$Species
    )
  }
  
  return(results)
}
