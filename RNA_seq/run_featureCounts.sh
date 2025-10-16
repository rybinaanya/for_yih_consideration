#!/bin/bash
# Usage:
#   ./run_featureCounts.sh -w /path/to/workdir -i /path/to/bowtie2_results -o /path/to/featureCounts_out
#


set -euo pipefail

# ----------- Parse arguments -----------
while getopts "w:i:o:h" opt; do
  case $opt in
    w) working_dir="$OPTARG" ;;
    i) bowtie2_folder="$OPTARG" ;;
    o) featureCounts_folder="$OPTARG" ;;
    h)
      echo "Usage: $0 -w <workdir> -i <bowtie2_input_dir> -o <output_dir>"
      exit 0
      ;;
    *)
      echo "Invalid option"; exit 1 ;;
  esac
done

# ----------- Check arguments -----------
if [[ -z "${working_dir:-}" || -z "${bowtie2_folder:-}" || -z "${featureCounts_folder:-}" ]]; then
  echo "Error: missing arguments."
  echo "Usage: $0 -w <workdir> -i <bowtie2_input_dir> -o <output_dir>"
  exit 1
fi

mkdir -p "$featureCounts_folder"
cd "$working_dir"

# ----------- Download GTF annotations -----------
echo "Downloading GTF annotations..."
declare -A gtfs
gtfs=(
  
["k12"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gtf.gz"
  
["nissle"]="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/967/895/GCF_019967895.1_ASM1996789v1/GCF_019967895.1_ASM1996789v1_genomic.gtf.gz"
)

for strain in "${!gtfs[@]}"; do
  url="${gtfs[$strain]}"
  out_gz="${bowtie2_folder}/$(basename "$url")"
  out_gtf="${out_gz%.gz}"

  wget -nc "$url" -P "$bowtie2_folder"
  if [[ -f "$out_gz" && ! -f "$out_gtf" ]]; then
    gunzip -k "$out_gz"
  fi
done

# ----------- Run featureCounts -----------
echo "Running featureCounts..."

featureCounts \
  -a "${bowtie2_folder}/GCF_000005845.2_ASM584v2_genomic.gtf" \
  -p -t gene -F GTF -f \
  --extraAttributes db_xref,gene \
  -T 5 \
  -o "${featureCounts_folder}/K12_merged.txt" \
  $(ls "${bowtie2_folder}"/K12*sorted.bam)

featureCounts \
  -a "${bowtie2_folder}/GCF_019967895.1_ASM1996789v1_genomic.gtf" \
  -p -t gene -F GTF -f \
  --extraAttributes db_xref,gene \
  -T 5 \
  -o "${featureCounts_folder}/Nissle_merged.txt" \
  $(ls "${bowtie2_folder}"/Nissle*sorted.bam)

echo "featureCounts completed. Results in: ${featureCounts_folder}"

