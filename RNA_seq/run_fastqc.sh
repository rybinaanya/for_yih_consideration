#!/bin/bash
# Usage:
#   ./run_fastqc.sh -r /path/to/reads -o /path/to/output
#
# Example:
#   ./run_fastqc.sh -r ./bbduk_out1 -o ./fastqc_bbduk_out1

set -e  # exit if any command fails

# ----------- Parse arguments -----------
while getopts "r:o:h" opt; do
  case $opt in
    r) reads_dir="$OPTARG" ;;
    o) out_dir="$OPTARG" ;;
    h)
      echo "Usage: $0 -r <reads_dir> -o <out_dir>"
      exit 0
      ;;
    *)
      echo "Invalid option"; exit 1 ;;
  esac
done

# ----------- Check arguments -----------
if [[ -z "$reads_dir" || -z "$out_dir" ]]; then
  echo "Error: missing arguments."
  echo "Usage: $0 -r <reads_dir> -o <out_dir>"
  exit 1
fi

# ----------- Create output directory -----------
mkdir -p "$out_dir"
cd "$reads_dir"

# ----------- Run FastQC -----------
for f in *.fastq.gz; do
  echo "Running FastQC on $f..."
  fastqc "$f" -o "$out_dir"
done
