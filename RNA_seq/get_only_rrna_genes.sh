#!/bin/bash

#============= Nissle 1917 ============#
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/967/895/GCF_019967895.1_ASM1996789v1/GCF_019967895.1_ASM1996789v1_rna_from_genomic.fna.gz

gunzip GCF_019967895.1_ASM1996789v1_rna_from_genomic.fna.gz


# Run Python script in the bash script
python3 << EOF
from Bio import SeqIO

# Input FASTA file
input_fasta = 'GCF_019967895.1_ASM1996789v1_rna_from_genomic.fna'

# Output FASTA file
output_fasta = 'rrna_nissle.fasta'

# List to hold sequences with 'rrn' in the description
seq_lst = []
for record in SeqIO.parse(input_fasta, 'fasta'):
    if 'rrn' in record.description:
        seq_lst.append(record)

# Write the filtered sequences to output FASTA file
with open(output_fasta, 'w') as fout:
    SeqIO.write(seq_lst, fout, 'fasta')

EOF

echo "Python script has been executed and rrna_k12.fasta has been created."

# ============ K-12 ==============#
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_rna_from_genomic.fna.gz

gunzip GCF_000005845.2_ASM584v2_rna_from_genomic.fna.gz


python3 << EOF
from Bio import SeqIO

# Input FASTA file
input_fasta = 'GCF_000005845.2_ASM584v2_rna_from_genomic.fna'

# Output FASTA file
output_fasta = 'rrna_k12.fasta'

# List to hold sequences with 'rrn' in the description
seq_lst = []
for record in SeqIO.parse(input_fasta, 'fasta'):
    if 'rrn' in record.description:
        seq_lst.append(record)

# Write the filtered sequences to output FASTA file
with open(output_fasta, 'w') as fout:
    SeqIO.write(seq_lst, fout, 'fasta')

EOF

echo "Python script has been executed and rrna_k12.fasta has been created."

