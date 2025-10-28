#!/bin/bash
set -e
# Run this script in conda env that has qiime2 + deicode + biom-format

indir="data/rpca_input"
outdir="data/rpca_output_90"

mkdir -p $outdir

for tsv in $indir/*_asv.tsv
do
  name=$(basename $tsv _asv.tsv)
  echo "Processing dataset: $name"

  # TSV → BIOM (HDF5 format)
  biom_path=$outdir/${name}_asv.biom
  biom convert \
    -i $tsv \
    -o $biom_path \
    --to-hdf5 \
    --table-type="OTU table"

  # BIOM → QZA
  qiime tools import \
    --input-path $biom_path \
    --type 'FeatureTable[Frequency]' \
    --input-format BIOMV210Format \
    --output-path $outdir/${name}_asv.qza

  # Run DEICODE RPCA
  qiime deicode rpca \
    --p-n-components 15 \
    --i-table $outdir/${name}_asv.qza \
    --o-biplot $outdir/${name}_rpca_biplot.qza \
    --o-distance-matrix $outdir/${name}_rpca_distance.qza

  echo "✅ Finished $name"
done

echo "All datasets processed!"
