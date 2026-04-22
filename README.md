# Salt-Marshes-Microbiome

This repository contains metadata, bioinformatics and statistical analysis scripts, and processed outputs associated with the study: “Compartment-dependent microbial reassembly reveals early recovery signals in restoring saltmarshes”.

## Overview

This study investigates bacterial (16S rRNA gene) and fungal (ITS) community assembly across saltmarsh restoration treatments. Microbial communities from Natural reference plots were compared with Managed restoration contexts, including BESE, BARE, and Degraded treatments, while also resolving differences between rhizosphere and bulk soil compartments.

## Repository structure

- `metadata/` — sample metadata and environmental variables used in the analyses
- `scripts/` — DADA2 processing scripts, downstream statistical analyses, and HPC job scripts
- `tables/` — processed ASV tables, taxonomy tables, QC summaries, and other statistical outputs


## Sequence data

Raw sequencing reads are deposited in the NCBI Sequence Read Archive (SRA):

- **BioProject:
- **SRA accession(s):

## Bioinformatics pipeline

Amplicon sequences were processed using the DADA2 workflow, including:

1. Quality filtering and trimming
2. Error model learning
3. ASV inference
4. Paired-end merging
5. Chimera removal
6. Taxonomic assignment

Intermediate files generated during DADA2 processing, such as filtered reads, dereplication objects, and learned error models, are not included in this repository because they can be regenerated from the raw data and scripts.

## Statistical analyses

Downstream analyses included:

- Alpha diversity (Chao1, Shannon, Simpson)
- Beta diversity based on Bray–Curtis dissimilarities
- Distance-based redundancy analysis (db-RDA)
- PERMANOVA
- Differential abundance analysis with DESeq2

Analyses were conducted in R using packages including:

`dada2`, `phyloseq`, `vegan`, `DESeq2`, `ggplot2`, `dplyr`, and related dependencies.

## Notes

- This repository contains lightweight, analysis-ready outputs and reproducible scripts.
- Large intermediate files are intentionally excluded.
- Raw sequence data should be obtained from NCBI SRA.
- The analyses in this repository can be reproduced from the raw data, metadata, and scripts provided here.
