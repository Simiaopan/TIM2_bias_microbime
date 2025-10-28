#!/bin/bash

set -e

##-------Set parameters-----------------------
projectdir=/data/tim2_microbiota
fastq_folder=$projectdir/np_fastq
fastqc_output=$projectdir/fastqc_output
conda_env=/data/tim2_microbiota/software/miniconda/etc/profile.d/conda.sh
f_prefix="R1"
r_prefix="R2"

source "$conda_env"
conda activate cut_qc 

##---------------- mkdir ------------------------------------
mkdir -p "$fastqc_output"
mkdir -p "$projectdir/temp"

##------ go into all project folders --------
for proj in "$fastq_folder"/*; do
  projname=$(basename "$proj")
  echo "📕 Processing project: $projname"

  # Storage of fastqc of each project
  mkdir -p "$fastqc_output/$projname"

  #------Get fastq pair count and set sampling amount--------
  files=($proj/*fastq.gz)
  total_files=${#files[@]}
  n_pairs=$((total_files / 2))

  if [ $n_pairs -le 5 ]; then
    n_sample=$n_pairs
  elif [ $n_pairs -le 20 ]; then
    n_sample=10
  elif [ $n_pairs -le 50 ]; then
    n_sample=15
  else
    n_sample=20
  fi

  echo "Total sample pairs: $n_pairs, Sampling: $n_sample per direction (R1 & R2)"

  ##----------------sampling and QC--------------------------
  for i in $f_prefix $r_prefix; do
    mkdir -p "$fastqc_output/$projname/raw_$i"

    # randomly sampling
    fq_rand=$(find "$proj" -type f -name "*$i*fastq.gz" | shuf -n $n_sample)

    # integrate
    echo "$fq_rand" | xargs zcat > "$projectdir/temp/$i.fastq"

    # run fastqc
    fastqc "$projectdir/temp/$i.fastq" \
      -o "$fastqc_output/$projname/raw_$i"
  done

  # clear temp files
  rm "$projectdir/temp/"*.fastq
done

conda deactivate