This folder contains basic scripts to process RNA-seq data (trimming, mapping, quantification) 



Quality of reads was assessed using FastQC v0.11.9 (`run_fastqc.sh`) and MultiQC v1.12.dev0. Reads were trimmed to remove adapters using BBduk v35.85 (`run_bbduk.sh`). The genome assemblies of *E. coli* K-12 MG1655 (GenBank accession: GCF_000005845.2) and *E. coli* Nissle 1917 (GenBank accession: GCF_019967895.1) were used as references.Reads were mapped to the reference data using Bowtie2 v2.2.1. Reads were decontaminated of rRNA (`extract_rrna.py`, `map2rrna.sh`). All processing of read alignment data was performed using Samtools v1.11. Read counts were obtained using featureCounts v2.0.1 
