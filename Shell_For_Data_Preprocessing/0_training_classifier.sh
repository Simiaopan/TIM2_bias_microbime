#!/bin/bash
set -e
mkdir -p tools

source /data/tim2_microbiota/software/miniconda/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# 1. Download SILVA 138
wget https://data.qiime2.org/2024.2/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2024.2/common/silva-138-99-tax.qza

# 2. Extract v34
qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-trunc-len 0 \
  --o-reads silva-138-v34-seqs.qza

# 3. Training
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva-138-v34-seqs.qza \
  --i-reference-taxonomy silva-138-99-tax.qza \
  --o-classifier tools/silva-138-v34-classifier.qza

conda deactivate


