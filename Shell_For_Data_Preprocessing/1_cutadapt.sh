#!/bin/bash

##------------------- set parameters ---------------------
projectdir=/data/tim2_microbiota
total_input=$projectdir/matched_fastq
inter_output=$projectdir/np_fastq
conda_env=/data/tim2_microbiota/software/miniconda/etc/profile.d/conda.sh
logfile=$projectdir/Logs/Log_Cutadapt.log
errorlog=$projectdir/Logs/Log_Wrong.log
n_cores=15
name_end="001.fastq.gz"
f_primer=CCTACGGGNGGCWGCAG
r_primer=GACTACHVGGGTATCTAATCC
f_prefix=R1
r_prefix=R2

##------------------ conda activate and make logs ----------------------
source "$conda_env"
conda activate cut_qc 

mkdir -p "$inter_output"
mkdir -p "$(dirname "$logfile")"

echo -e "\n [$(date)] Starting cutadapt processing...\n" >> "$logfile"

##---------------- find files and cut adapter --------------------------

for project_path in "$total_input"/*; do
    if [ -d "$project_path" ]; then
        project_name=$(basename "$project_path")
        input_folder="$total_input/$project_name"
        output_folder="$inter_output/$project_name"
        report_folder="$projectdir/QC/cut_report/$project_name"
        mkdir -p "$output_folder"
        mkdir -p "$report_folder"

        echo -e "\n Project: $project_name" | tee -a "$logfile"

        # find all R1 fastq (basename without _R1_...)
        for file in "$input_folder"/*_"$f_prefix"_${name_end}; do
            [ -e "$file" ] || continue  # if no matched files, skip
            base_name=$(basename "$file" | sed "s/_${f_prefix}_.*//")

            out_r1="${output_folder}/np${base_name}_${f_prefix}_${name_end}"
            out_r2="${output_folder}/np${base_name}_${r_prefix}_${name_end}"

            # if treated files already existed, slip
            if [[ -f "$out_r1" && -f "$out_r2" ]]; then
                echo "✅ Sample $base_name already processed. Skipping..." | tee -a "$logfile"
                continue
            fi

            r1_file="${input_folder}/${base_name}_${f_prefix}_${name_end}"
            r2_file="${input_folder}/${base_name}_${r_prefix}_${name_end}"

            echo -e "\n▶️  Processing sample: ${base_name}" | tee -a "$logfile"
            echo "📂 R1: $r1_file" | tee -a "$logfile"
            echo "📂 R2: $r2_file" | tee -a "$logfile"

            # check if the file exist
            if [[ ! -f "$r1_file" || ! -f "$r2_file" ]]; then
                echo "❌ Missing file(s) for $base_name. Skipping..." | tee -a "$logfile"
                echo "[${project_name}] Sample $base_name: missing R1 or R2" >> "$errorlog"
                continue
            fi

            start_time=$(date +%s)

            # run cutadapt and cature error
            if cutadapt \
                --cores=$n_cores \
                --discard-untrimmed \
                -g ^$f_primer \
                -G ^$r_primer \
                -o "$out_r1" \
                -p "$out_r2" \
                "$r1_file" "$r2_file" \
                > "${report_folder}/${base_name}.txt" 2>> "$errorlog"
            then
                end_time=$(date +%s)
                duration=$((end_time - start_time))
                echo "✅ Finished $base_name in ${duration}s" | tee -a "$logfile"
            else
                echo "❌ cutadapt failed for $base_name" | tee -a "$logfile"
                echo "[${project_name}] Sample $base_name: cutadapt failed" >> "$errorlog"
                continue
            fi

        done
    fi
done

echo -e "\n✅ All projects finished at $(date)." | tee -a "$logfile"

conda deactivate
