#!/bin/bash
# Usage:
#   ./getUnmappedToRrnaReads.sh -r /path/to/reads -w /path/to/workdir -o /path/to/output
#
# Example:
#   ./getUnmappedToRrnaReads.sh -r ./bbduk_out1 -w ./RNAseq -o ./bowtie2_out_rrna

set -euo pipefail

# ----------- Parse arguments -----------
while getopts "r:w:o:h" opt; do
  case $opt in
    r) reads_dir="$OPTARG" ;;
    w) working_dir="$OPTARG" ;;
    o) bowtie2_folder="$OPTARG" ;;
    h)
      echo "Usage: $0 -r <reads_dir> -w <working_dir> -o <output_dir>"
      exit 0
      ;;
    *)
      echo "Invalid option"; exit 1 ;;
  esac
done

# ----------- Check arguments -----------
if [[ -z "${reads_dir:-}" || -z "${working_dir:-}" || -z "${bowtie2_folder:-}" 
]]; then
  echo "Error: missing arguments."
  echo "Usage: $0 -r <reads_dir> -w <working_dir> -o <output_dir>"
  exit 1
fi

mkdir -p "$bowtie2_folder"
cd "$working_dir"

echo "Extracting unmapped reads from BAM files..."
for f in "${reads_dir}"/*R1_bbduk_001.fastq.gz; do
  file=$(basename "$f")
  strain=$(echo "$file" | cut -f1 -d"_")
  sugar=$(echo "$file" | cut -f2 -d"_")
  prefix="${strain}_${sugar}"

  echo "-------------------------------------------"
  echo "Processing ${prefix}"

  unmapped_bam="${bowtie2_folder}/${prefix}_rrna.unmapped.bam"
  sorted_bam="${bowtie2_folder}/${prefix}_rrna.unmapped.sortedByName.bam"
  
fixmate_bam="${bowtie2_folder}/${prefix}_rrna.unmapped.sortedByName.fixmate.bam"
  out_R1="${bowtie2_folder}/${prefix}_R1.clean.fq"
  out_R2="${bowtie2_folder}/${prefix}_R2.clean.fq"

  if [[ ! -f "$unmapped_bam" ]]; then
    echo "Missing file: $unmapped_bam — skipping"
    continue
  fi

  echo "Sorting by name..."
  samtools sort -n "$unmapped_bam" -o "$sorted_bam"

  echo "Fixing mate information..."
  samtools fixmate "$sorted_bam" "$fixmate_bam"

  echo "Converting BAM → FASTQ..."
  bedtools bamtofastq -i "$fixmate_bam" -fq "$out_R1" -fq2 "$out_R2"

  echo "Done: ${prefix}"
done

echo "All unmapped reads extracted and converted!"

