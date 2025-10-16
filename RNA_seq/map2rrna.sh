#!/bin/bash
# Usage:
#   ./remove_rrna.sh -r /path/to/reads -w /path/to/workdir -o /path/to/output
#
# Example:
#   ./remove_rrna.sh -r ./bbduk_out1 -w ./RNAseq -o ./bowtie2_out_rrna

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

# ----------- Ensure rrna FASTA exist -----------
for strain in k12 nissle; do
  fasta="rrna_${strain}.fasta"
  if [[ ! -f "$fasta" ]]; then
    echo "Missing $fasta. Please prepare it using extract_rrna.py"
    exit 1
  fi
  cp "$fasta" "$bowtie2_folder/"
done

# ----------- Build bowtie2 indices -----------
echo "Building bowtie2 indices..."
bowtie2-build "$bowtie2_folder/rrna_k12.fasta" "$bowtie2_folder/k12_rrna"
bowtie2-build "$bowtie2_folder/rrna_nissle.fasta" 
"$bowtie2_folder/nissle_rrna"

# ----------- Main loop over samples -----------
echo "Processing samples from: $reads_dir"
for f in "${reads_dir}"/*R1_bbduk_001.fastq.gz; do
  file=$(basename "$f")
  strain=$(echo "$file" | cut -f1 -d"_")
  sugar=$(echo "$file" | cut -f2 -d"_")
  pair_R2="${f/R1/R2}"
  prefix="${strain}_${sugar}"

  echo "-------------------------------------------"
  echo "Sample: $prefix"
  echo "R1: $file"
  echo "R2: $(basename "$pair_R2")"

  # Choose rrna index
  if [[ "$strain" == "K12" ]]; then
    index="k12_rrna"
  elif [[ "$strain" == "Nissle" ]]; then
    index="nissle_rrna"
  else
    echo "️ Unknown strain name in $file — skipping"
    continue
  fi

  # Run bowtie2 alignment
  echo " Running bowtie2 for $strain ..."
  bowtie2 -x "${bowtie2_folder}/${index}" \
          -1 "$f" -2 "$pair_R2" \
          -S "${bowtie2_folder}/${prefix}_rrna.sam" \
          --threads 8

  # Convert SAM → BAM
  samtools view -bS "${bowtie2_folder}/${prefix}_rrna.sam" \
    > "${bowtie2_folder}/${prefix}_rrna.bam"

  # Sort BAM
  samtools sort "${bowtie2_folder}/${prefix}_rrna.bam" \
    -o "${bowtie2_folder}/${prefix}_rrna.sorted.bam"
  samtools index "${bowtie2_folder}/${prefix}_rrna.sorted.bam"

  # Extract unmapped reads
  echo " Extracting unmapped reads..."
  samtools view -b -f 4 "${bowtie2_folder}/${prefix}_rrna.sorted.bam" \
    > "${bowtie2_folder}/${prefix}_rrna.unmapped.bam"
  samtools sort -n "${bowtie2_folder}/${prefix}_rrna.unmapped.bam" \
    -o "${bowtie2_folder}/${prefix}_rrna.unmapped.sortedByName.bam"

  # Convert to FASTQ
  bedtools bamtofastq \
    -i "${bowtie2_folder}/${prefix}_rrna.unmapped.sortedByName.bam" \
    -fq "${bowtie2_folder}/${prefix}_R1.clean.fq" \
    -fq2 "${bowtie2_folder}/${prefix}_R2.clean.fq"

  echo "Completed: $prefix"
done

echo "All samples processed successfully!"

