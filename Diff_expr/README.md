The current folder contains minimal set of scripts used to perform differential expression analysis and functional characterization of differentially expressed genes (DEGs).

R Markdown scripts:
* `run_gseKegg.Rmd`: generates gene expression heatmaps of top-enriched KEGG pathways in *E. coli* K-12 and Nissle 1917 strains (`Figure 6`)
* `plot_qrtPCR_growth.Rmd`: plots results of qRT-PCR (`Figure 4`) and growth curves experiments (`Figure 5`)


Data:
* `growth_curves.tsv` containd OD600 values from bacterial growth experiments and serves as input data for the `plot_qrtPCR_growth.Rmd` script
* `qrt_pcr_2**ddCT.tsv` and `qrt_pcr_extended.txt` provide data on gene expression values obtained by qRT-PCR and serve as input data for the `plot_qrtPCR_growth.Rmd` script
* `K12_merged` and `Nissle_merged` are gene count matrices obtained using the featureCounts software and serve as input data for the `runDESeq2.Rmd` script
* `K12_SQ_vs_gluc_results_padj0.05_DAVID_operonMapper.tsv` and `Nissle_SQ_vs_gluc_results_padj0.05_DAVID_operonMapper_extended.tsv` are dataframes that combine results from DESeq2, operon mapper and DAVID analysis; serve as input data for the `runDESeq2.Rmd` script
* `UpSet_data_Nissle_vs_K12.csv` contains log2FoldChange values and gene IDs for both strains and serves as input data for the `run_gseKegg.Rmd` script
