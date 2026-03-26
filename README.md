# Salt-Marshes-Microbiome

This repository contains metadata, analysis scripts, and processed outputs for the study:

**"Compartment-dependent microbial reassembly reveals early recovery signals in restoring saltmarshes"**

## Overview

This study investigates bacterial (16S) and fungal (ITS) community assembly across saltmarsh restoration treatments, comparing Natural reference conditions with Managed treatments (BESE, BARE, Degraded), and resolving patterns across soil compartments (rhizosphere vs bulk soil).

## Repository structure

- `metadata/` → sample metadata used in analyses  
- `scripts/` → R scripts for DADA2 processing and downstream analyses  
- `tables/` → processed data and statistical results  
- `figures/` → figures used in the manuscript  
- `docs/` → workflow description  

## Sequence data

Raw sequencing data are available at:

→ NCBI SRA: 

## Bioinformatics pipeline

Sequences were processed using the DADA2 pipeline, including:

1. Quality filtering and trimming  
2. Error model learning  
3. ASV inference  
4. Paired-end merging  
5. Chimera removal  
6. Taxonomic assignment  

Intermediate files (filtered reads, error models, etc.) are not included but can be generated using the scripts provided.

## Statistical analyses

- Alpha diversity (Chao1, Shannon, Simpson)  
- Beta diversity (Bray–Curtis, db-RDA, PERMANOVA)  
- Differential abundance (DESeq2)  

All analyses were conducted in R using packages including:
`phyloseq`, `vegan`, `DESeq2`, `ggplot2`

## Notes

- Large intermediate files are not included in this repository  
- All results can be reproduced from raw data and scripts  
