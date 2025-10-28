#!/bin/bash
set -e

#------ input and output ------------------------
src_dir=/data/tim2_microbiota/np_fastq
dst_dir=/data/tim2_microbiota/flat_fastq

rm -rf "$dst_dir"
mkdir -p "$dst_dir"

echo "Copying FASTQ files (this may take some time depending on size)..."
find "$src_dir" -type f -name "*.fastq.gz" -exec cp {} "$dst_dir/" \;

echo "🔍 Previewing first 10 files (sorted by filename):"
ls "$dst_dir" | sort | head -n 10

echo "✅ Done. All FASTQ files copied to $dst_dir"
