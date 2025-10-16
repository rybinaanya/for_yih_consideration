This folder contains basic scripts to annotate the phylogenetic distribution of *yih* cassette variants across *E. coli* strains.

**Scripts**:
* `run_mmseqs2.py` identifies homologous sequences of _yih_ genes (from _E. coli_ K12) in various _E. coli_ genomes from the Seferbekova et al. (2021) study (ShigellaProject repository) using MMseqs2; the script downloads ShigellaProject genomes, runs MMseqs2 and parses results into final table `yih_genes_mmseq.tsv`
* `get_yih_loci_data.ipynb` analyzes the genomic organization and colocalization of _yih_ locus genes in different _E. coli_ strains; the script takes  `yih_genes_mmseq.tsv` as input and generates output with locus composition and genomic coordinates: `yih_homologs_loci.csv` and `target2info_data_altToGFFs.csv`
* `get_itol_files.ipynb` generates iTOL (Interactive Tree Of Life) annotation files for visualizing _yih_ operon genomic data on phylogenetic tree of _E. coli_ strains from ShigellaProject; the script creates two types of iTOL annotation files: domain annotation file (`Ecoli.domains.txt`) that visualizes genomic organization of _yih_ loci, and binary presence/absence gile (`yih_loci_presence.txt`) that visualizes gene content patterns across strains
* `run_Fisher_test.ipynb` performs statistical analysis to test the association between _yih_ locus variants and bacterial pathogenicity using Fisher's Exact Test.


**Input data**:
*`nt_tree_midpoint.nwk` - Phylogenetic tree of _E. coli_ strains from ShigellaProject (Seferbekova et al. 2021)
*`clusters.txt` - Phylogroup coloring for the tree from ShigellaProject (Seferbekova et al. 2021)
*`labels.txt` - Organism names for tree nodes from ShigellaProject (Seferbekova et al. 2021)
*`total_stats.csv` - Strain statistics and metadata from ShigellaProject (Seferbekova et al. 2021)

**Analysis outputs**:
*`yih_genes_mmseq.tsv` - MMseqs2 homology results (from `run_mmseqs2.py`)
*`yih_homologs_loci.csv` - _yih_ locus composition and genomic coordinates (from `get_yih_loci_data.ipynb`)
*`target2info_data_altToGFFs.csv` - GFF-like annotations for all _yih_ homologs (from `get_yih_loci_data.ipynb`)

**iTOL visualization files**:
`yih_loci_presence.txt` - binary gene presence/absence patterns (from` get_itol_files.ipynb`)
`Ecoli.domains.txt` - genomic organization of _yih_ loci (from `get_itol_files.ipynb`)
