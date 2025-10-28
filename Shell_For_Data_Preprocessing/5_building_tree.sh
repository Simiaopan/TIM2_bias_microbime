#!/bin/bash
set -e

##----------parameters setting----------------
projectdir=/data/tim2_microbiota
input=$projectdir/output/dada2
output=$projectdir/output/export/trees
conda_env=/data/tim2_microbiota/software/miniconda/etc/profile.d/conda.sh
n_cores=15

##------- conda activate and mkdir --------------
source "$conda_env"
conda activate qiime2-amplicon-2024.10
mkdir -p $output

qiime alignment mafft \
  --i-sequences $input/rep_seqs.qza \
  --o-alignment $output/aligned-rep-seqs.qza \
  --p-n-threads $n_cores

qiime alignment mask \
  --i-alignment $output/aligned-rep-seqs.qza \
  --o-masked-alignment $output/masked-aligned-rep-seqs.qza

# unrooted-tree
qiime phylogeny fasttree \
  --i-alignment $output/masked-aligned-rep-seqs.qza \
  --o-tree $output/unrooted-tree.qza

# rooted-tree
qiime phylogeny midpoint-root \
  --i-tree $output/unrooted-tree.qza \
  --o-rooted-tree $output/rooted-tree.qza

##----------------------------
conda deactivate
