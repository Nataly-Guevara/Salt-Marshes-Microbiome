# Salt-Marshes-Microbiome

Scripts, metadata, processed sequencing data, and statistical analysis results associated with the saltmarsh restoration microbiome study **“Microbiome composition during early saltmarsh restoration reflects nested environmental, treatment and root-associated structuring.”**

Bacterial and archaeal 16S rRNA gene and fungal ITS communities were studied across three Victorian saltmarshes, comparing the Natural reference, BESE, BARE and Degraded ecological conditions, and bulk soil (BS) with rhizosphere (RZ) where applicable.

## Repository structure

```text
scripts/DADA2/16S/                  Amplicon bioinformatics scripts
scripts/DADA2/ITS/
docs/krona/                         Interactive taxonomic visualisation
metadata/                           Sequencing, soil and vegetation metadata
analyses/soil_environment/          Soil physicochemical analyses
analyses/vegetation/                Vegetation analyses
analyses/alpha_diversity/           Alpha diversity and sequencing depth
analyses/beta_diversity/            Community composition and ordination
analyses/differential_abundance/    ANCOM-BC2 differential abundance
analyses/community_assembly/        Community assembly
analyses/shared_workflows/          Scripts covering multiple analyses
data/processed/                     Analysis-ready datasets
tables/ASV/, tables/taxonomy/       Existing DADA2-related sequence tables
tables/phyloseq/, tables/qc/        Existing phyloseq objects and QC
```

Each analysis folder contains `scripts/`, `inputs/`, and `outputs/` locations; some inputs/outputs are still being transferred from the analysis computers.

## Analysis comparisons

- **H1 (bulk soil):** Natural, BESE, BARE and Degraded.
- **H2 (rhizosphere):** Natural, BESE and BARE.
- **H3 (paired habitats):** Bulk soil and rhizosphere within Natural, BESE and BARE.

Statistical workflows include diversity analyses, Bray–Curtis community composition and reference-similarity tests, and **ANCOM-BC2** differential abundance. Multi-method H1/H2/H3 scripts are in [analyses/shared_workflows/scripts](analyses/shared_workflows/scripts/). DESeq2 is **not** part of the final study.

## Sequencing

Raw reads: NCBI SRA [BioProject PRJNA1446205](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1446205). Reads are not duplicated in this repository.

## Data and reproducibility

The available R scripts and DADA2/Krona assets are committed. **The transfer of the latest local analysis inputs and outputs is still in progress**, and the complete set of scripts has not been rerun from this repository. Some older DADA2-derived data files under `tables/` also require cross-checking against the final analysis-ready datasets. This README will be updated when the local, Gandalf and Gadi files have been synced.
