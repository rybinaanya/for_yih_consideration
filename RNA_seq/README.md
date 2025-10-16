This folder contains basic scripts to process RNA-seq data (trimming, mapping, quantification) 


Example of commands:
```{bash}
./run_bbduk.sh \
  -r /mnt/data/RNAseq/fastq/ \
  -a /mnt/data/adapters.fa \
  -o ./results/bbduk_out
```
