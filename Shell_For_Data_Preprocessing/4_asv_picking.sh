#!/bin/bash
set -e

##----------parameters setting----------------
projectdir=/data/tim2_microbiota
total_input=$projectdir/flat_fastq
total_output=$projectdir/output
conda_env=/data/tim2_microbiota/software/miniconda/etc/profile.d/conda.sh

tax_classif=$projectdir/tools/silva-138-v34-classifier.qza
n_cores=15

# DADA2
p_trim_left_f=0
p_trim_left_r=0
p_trunc_len_f=250
p_trunc_len_r=215
p_max_ee_f=2
p_max_ee_r=2

# Taxonomy minimum confidence
tax_min_conf=0.7

# out put directory setting
mkdir -p $total_output
mkdir -p $total_output/dada2

out_folder=$total_output/dada2
mkdir -p $out_folder
mkdir -p $out_folder/summary
mkdir -p $out_folder/export

log_file="$projectdir/Log_DADA.log"
echo "🔔 DADA2 pipeline started at $(date)" | tee "$log_file"

##--------- conda activate ---------------------------
source "$conda_env"
conda activate qiime2-amplicon-2024.10

##----------- veiw the first 10 fastq files ---------------------------
echo "🔍 Previewing first 10 FASTQ files (sorted):" | tee -a "$log_file"
find "$total_input" -type f -name "*.fastq.gz" | sort | head -n 10 | tee -a "$log_file"

##-----------Import data using Casava format---------------------------
echo "📦 Importing FASTQ files using Casava format..." | tee -a "$log_file"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$total_input" \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path "$out_folder/demux-paired-end.qza" \
  2>&1 | tee -a "$log_file"

##-----------Demux summary---------------------------
echo "📈 Summarizing demux output..." | tee -a "$log_file"
qiime demux summarize \
  --i-data "$out_folder/demux-paired-end.qza" \
  --o-visualization "$out_folder/summary/reads_summary.qzv" \
  2>&1 | tee -a "$log_file"

##-----------DADA2 denoise-----------------------------
echo "🚀 Running DADA2 denoising (this might take a while)..." | tee -a "$log_file"
temp_dada_log="$out_folder/summary/dada2_time_memory.log"

/usr/bin/time -v \
  qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$out_folder/demux-paired-end.qza" \
  --p-trim-left-f $p_trim_left_f \
  --p-trim-left-r $p_trim_left_r \
  --p-trunc-len-f $p_trunc_len_f \
  --p-trunc-len-r $p_trunc_len_r \
  --p-max-ee-f $p_max_ee_f \
  --p-max-ee-r $p_max_ee_r \
  --p-n-threads $n_cores \
  --p-n-reads-learn 3000000 \
  --p-pooling-method 'pseudo' \
  --o-table "$out_folder/asv_table.qza" \
  --o-representative-sequences "$out_folder/rep_seqs.qza" \
  --o-denoising-stats "$out_folder/summary/denoising_stats.qza" \
  2> "$temp_dada_log"

# veiw the largest RAM usage
max_mem_kb=$(grep "Maximum resident set size" "$temp_dada_log" | awk '{print $6}')
max_mem_mb=$((max_mem_kb / 1024))

if [ "$max_mem_mb" -gt 65536 ]; then
  echo "⚠️ [WARNING] DADA2 used ${max_mem_mb}MB memory, which exceeds the 64GB threshold!" | tee -a "$log_file"
else
  echo "✅ DADA2 memory usage: ${max_mem_mb}MB (within 64GB limit)" | tee -a "$log_file"
fi

# veiw total dada time
elapsed_time=$(grep "Elapsed (wall clock) time" "$temp_dada_log" | awk -F': ' '{print $2}')
echo "⏱  DADA2 runtime: $elapsed_time" | tee -a "$log_file"

# DADA2 statastic visualization
qiime metadata tabulate \
  --m-input-file "$out_folder/summary/denoising_stats.qza" \
  --o-visualization "$out_folder/summary/denoising_stats.qzv" \
  2>&1 | tee -a "$log_file"

##----------- export DADA data ----------------------------
echo "📤 Exporting ASV table and rep_seqs..." | tee -a "$log_file"
qiime tools export \
  --input-path "$out_folder/asv_table.qza" \
  --output-path "$out_folder/export/asv_table"

qiime tools export \
  --input-path "$out_folder/rep_seqs.qza" \
  --output-path "$out_folder/export/rep_seqs"

##----------- Taxonomy -------------------------------
conf_name_add=$(echo "$tax_min_conf" | sed "s|\.||g")
echo "📃 Taxonomy classification..." | tee -a "$log_file"

qiime feature-classifier classify-sklearn \
  --i-classifier "$tax_classif" \
  --i-reads "$out_folder/rep_seqs.qza" \
  --p-n-jobs "$n_cores" \
  --p-confidence "$tax_min_conf" \
  --o-classification "$out_folder/taxonomy_${conf_name_add}.qza" \
  2>&1 | tee -a "$log_file"

qiime metadata tabulate \
  --m-input-file "$out_folder/taxonomy_${conf_name_add}.qza" \
  --o-visualization "$out_folder/summary/taxonomy_${conf_name_add}.qzv" \
  2>&1 | tee -a "$log_file"

qiime tools export \
  --input-path "$out_folder/taxonomy_${conf_name_add}.qza" \
  --output-path "$out_folder/export/taxonomy"

##---------------------------------
echo "✅ Pipeline completed at $(date)" | tee -a "$log_file"
conda deactivate
