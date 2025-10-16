#!/bin/bash
#PBS -d .
#PBS -l walltime=100:00:00,mem=20gb
# Usage:
#   ./run_bbduk.sh -r /path/to/reads -a /path/to/adapters.fa -o 
/path/to/output

set -e  # exit if any command fails

# ----------- Parse arguments -----------
while getopts "r:a:o:h" opt; do
  case $opt in
    r) reads_dir="$OPTARG" ;;
    a) adapters_fa="$OPTARG" ;;
    o) out_dir="$OPTARG" ;;
    h)
      echo "Usage: $0 -r <reads_dir> -a <adapters_fa> -o <out_dir>"
      exit 0
      ;;
    *)
      echo "Invalid option"; exit 1 ;;
  esac
done

# ----------- Check arguments -----------
if [[ -z "$reads_dir" || -z "$adapters_fa" || -z "$out_dir" ]]; then
  echo "Error: missing arguments."
  echo "Usage: $0 -r <reads_dir> -a <adapters_fa> -o <out_dir>"
  exit 1
fi

mkdir -p "$out_dir"
cd "$reads_dir"

# ----------- Run BBDuk -----------
for f in *R1_001.fastq.gz; do
  sample=${f%_001.fastq.gz*}
  echo "Processing $sample..."
  
  bbduk.sh -Xmx1g \
    in1="$f" in2="${f/R1/R2}" \
    out1="${out_dir}/${sample}_bbduk_001.fastq.gz" \
    out2="${out_dir}/${sample/R1/R2}_bbduk_001.fastq.gz" \
    ref="$adapters_fa" \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo
done

