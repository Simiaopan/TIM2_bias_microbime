#-------------------------------------------------------------------------------
# Extract coords + meta from RPCA results for distance analysis
#-------------------------------------------------------------------------------
extract_rpca_coords_meta <- function(rpca_results, dataset_name) {
  stopifnot(dataset_name %in% names(rpca_results))
  
  biplot_df <- rpca_results[[dataset_name]]$biplot
  
  coords_all <- biplot_df %>%
    tibble::column_to_rownames("SampleID") %>%
    dplyr::select(starts_with("PC"))
  
  meta <- biplot_df %>%
    tibble::column_to_rownames("SampleID") %>%
    dplyr::select(-starts_with("PC"))
  
  return(list(coords_all = coords_all, meta = meta))
}
