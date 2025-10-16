The current folder contains minimal set of scripts used to perform differential expression analysis and functional characterization of differentially expressed genes (DEGs).

R Markdown scripts:
* `run_gseKegg.Rmd` to generate gene expression heatmaps of top-enriched KEGG pathways in *E. coli* K-12 and Nissle 1917 strains (`Figure 6`)
* `plot_qrtPCR_growth.Rmd`to plot results of qRT-PCR (`Figure 4`) and growth curves experiments (`Figure 5`)
* `runDESeq2.Rmd` to assess quality of RNA-seq samples (`Figure S1`), perform differential expression analysis, and generate volcano plots (`Figure 3`)
* `plot_KEGG_pathways.Rmd` to generate KEGG pathway maps (`Figure S2- S5`) 


Data:
* `growth_curves.tsv` — contains OD₆₀₀ measurements from bacterial growth experiments; used as input for `plot_qrtPCR_growth.Rmd` (`Figure 5`).
* `qrt_pcr_2**ddCT.tsv` and `qrt_pcr_extended.txt` — contain qRT-PCR–derived gene expression values; used by `plot_qrtPCR_growth.Rmd` to plot `Figure 4`.
* `K12_merged` and `Nissle_merged` — gene count matrices generated with featureCounts; used as input for `runDESeq2.Rmd` (`Figures 3` and `S1`).
* `K12_SQ_vs_gluc_results_padj0.05_DAVID_operonMapper.tsv` and `Nissle_SQ_vs_gluc_results_padj0.05_DAVID_operonMapper_extended.tsv` — data tables combining DESeq2, Operon Mapper, and DAVID outputs; used by `runDESeq2.Rmd` for functional characterization of DEGs.
* `UpSet_data_Nissle_vs_K12.csv` — contains log₂ fold change values and gene IDs for both strains; used as input for `run_gseKegg.Rmd` (`Figure 6`).
* `locus_tags_Nissle2K12_mmseqs2_rbh.tx`t — result of an MMseqs2 reciprocal best-hit search providing gene correspondence between _E. coli_ K-12 MG1655 and _E. coli_ Nissle 1917; used by `plot_KEGG_pathways.Rmd` (`Figures S2–S5`).
