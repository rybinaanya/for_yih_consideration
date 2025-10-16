#!/bin/bash
# Usage:
#   ./run_bowtie2_wo_rrna.sh -r /path/to/reads -w /path/to/workdir -o /path/to/output
#


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
suffix="wo_rrna_f4fixmate"

# ----------- Download reference genomes -----------
echo "Downloading reference genomes..."
declare -A refs
refs=(
  
["k12"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
  
["nissle"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/967/895/GCF_019967895.1_ASM1996789v1/GCF_019967895.1_ASM1996789v1_genomic.fna.gz"
)

for strain in "${!refs[@]}"; do
  url="${refs[$strain]}"
  out_gz="${bowtie2_folder}/$(basename "$url")"
  out_fa="${out_gz%.gz}"

  wget -nc "$url" -P "$bowtie2_folder"
  if [[ -f "$out_gz" && ! -f "$out_fa" ]]; then
    gunzip -k "$out_gz"
  fi
done

# ----------- Build Bowtie2 indices -----------
echo "Building bowtie2 indices..."
bowtie2-build "${bowtie2_folder}/GCF_000005845.2_ASM584v2_genomic.fna" 
"${bowtie2_folder}/ref_k12"
bowtie2-build "${bowtie2_folder}/GCF_019967895.1_ASM1996789v1_genomic.fna" 
"${bowtie2_folder}/ref_nissle"

# ----------- Iterate over read samples -----------
echo "Aligning reads..."
for f in "${reads_dir}"/*unmapped.sortedByName.fixmate.R1.fq; do
  file=$(basename "$f")
  strain=$(echo "$file" | cut -f1 -d"_")
  sugar=$(echo "$file" | cut -f2 -d"_")
  prefix="${strain}_${sugar}"

  echo "-------------------------------------------"
  echo "Sample: $prefix"
  pair_R2="${f/R1/R2}"

  if [[ "$strain" == "K12" ]]; then
    ref="${bowtie2_folder}/ref_k12"
  elif [[ "$strain" == "Nissle" ]]; then
    ref="${bowtie2_folder}/ref_nissle"
  else
    echo "️ Unknown strain in $file — skipping"
    continue
  fi

  echo "Running bowtie2 for $strain ($sugar)..."
  bowtie2 --local -x "$ref" \
          -1 "$f" -2 "$pair_R2" \
          -S "${bowtie2_folder}/${prefix}_${suffix}.sam" \
          --threads 8

  echo "Converting SAM → BAM..."
  samtools view -S "${bowtie2_folder}/${prefix}_${suffix}.sam" -b -o "${bowtie2_folder}/${prefix}_${suffix}.bam"

  echo "Sorting & indexing..."
  samtools sort "${bowtie2_folder}/${prefix}_${suffix}.bam" -o "${bowtie2_folder}/${prefix}_${suffix}.sorted.bam"
  samtools index "${bowtie2_folder}/${prefix}_${suffix}.sorted.bam"

  echo "Done: ${prefix}"
done

echo "All samples aligned successfully!"

