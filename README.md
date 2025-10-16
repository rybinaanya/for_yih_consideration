# Sulfoquinovose catabolism in *E. coli* strains: compositional and functional divergence of *yih* gene cassettes

This repository accompanies the manuscript  
*“Sulfoquinovose catabolism in E. coli strains: compositional and functional divergence of yih gene cassettes”*,  
currently under review at *International Journal of Molecular Sciences*.

It includes a minimal set of scripts used to:

- annotate the phylogenetic distribution of *yih* cassette variants across *E. coli* strains (`phylogeny/` folder)
- process RNA-seq data (trimming, mapping, quantification) (`RNA_seq/` folder),
- and perform differential expression analysis (DESeq2, clusterProfiler, etc) (`Diff_expr/` folder).

## Data availability

Raw RNA-seq reads are available from the NCBI Sequence Read Archive under BioProject accession [PRJNA1338229](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1338229).

The phylogenetic tree used in this study was obtained from Seferbekova et al. (2021):  
**High Rates of Genome Rearrangements and Pathogenicity of *Shigella* spp.** *Front. Microbiol.*  
[https://doi.org/10.3389/fmicb.2021.628622](https://doi.org/10.3389/fmicb.2021.628622)

Tree annotations and metadata were generated as part of this repository (see `phylogeny/` folder).
