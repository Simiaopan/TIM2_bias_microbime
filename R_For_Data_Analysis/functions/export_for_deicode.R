#-------------------------------------------------------------------------------
# Prepare data for rPCA (deicode) - export as TSV + meta
#-------------------------------------------------------------------------------
export_for_deicode <- function(ps_list, outdir) {
  library(phyloseq)
  
  # make output dir
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  for (dataset_name in names(ps_list)) {
    ps_obj <- ps_list[[dataset_name]]
    
    # extract ASV table
    otu <- as(otu_table(ps_obj), "matrix")
    if (!taxa_are_rows(ps_obj)) otu <- t(otu)
    
    # extract meta
    meta <- as(sample_data(ps_obj), "data.frame")
    
    # save meta
    saveRDS(meta, file = file.path(outdir, paste0(dataset_name, "_meta.rds")))
    
    # save OTU as TSV
    tsv_path <- file.path(outdir, paste0(dataset_name, "_asv.tsv"))
    write.table(
      cbind(FeatureID = rownames(otu), otu),
      file = tsv_path,
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    
    # print message for next step (new pipeline: TSV → biom-convert → QZA → RPCA)
    message("Exported ", dataset_name, " → ", tsv_path)
    message("Run in qiime2 env:")
    message("  biom convert -i ", tsv_path,
            " -o ", file.path(outdir, paste0(dataset_name, "_asv.biom")),
            " --to-hdf5 --table-type=\"OTU table\"")
    message("  qiime tools import --input-path ", 
            file.path(outdir, paste0(dataset_name, "_asv.biom")),
            " --type 'FeatureTable[Frequency]' --input-format BIOMV210Format",
            " --output-path ", file.path(outdir, paste0(dataset_name, "_asv.qza")))
    message("  qiime deicode rpca --i-table ", 
            file.path(outdir, paste0(dataset_name, "_asv.qza")),
            " --o-biplot ", file.path(outdir, paste0(dataset_name, "_rpca_biplot.qza")),
            " --o-distance-matrix ", file.path(outdir, paste0(dataset_name, "_rpca_distance.qza")))
    message("------------------------------------------------------------")
  }
}
