# Chapter 1 repository refresh — staging checklist (21 October 2026)

This document tracks the migration from legacy Chapter 1 analyses to manuscript-aligned, reviewed scripts, inputs, and outputs from local computer, Gandalf, and NCI Gadi.

**Staging branch:** `chapter1-2026-refresh` (not yet merged into `main`).

## Preserve without alteration
- `scripts/DADA2/16S/`
- `scripts/DADA2/ITS/`
- `metadata/dada2/`
- `docs/krona/index.html`

## Review before replacing
- `scripts/DESeq2/` and `tables/deseq2/`: superseded by the final differential-abundance workflow (AMBCOM, subject to verifying source scripts and outputs).
- `scripts/alpha_diversity/` and `tables/alpha_diversity/`: reconcile against revised manuscript.
- `scripts/beta_diversity/` and `tables/beta_diversity /`: replace obsolete models using verified final scripts and outputs. Note the existing directory name includes a trailing space.
- `scripts/vegetation_cover_color_height/` and environmental/vegetation metadata: check versions.
- `tables/ASV/`, `tables/taxonomy/`, `tables/phyloseq/`: verify scientific provenance and content before altering.

## Proposed structure for final review
```text
README.md
docs/krona/
metadata/dada2/
metadata/environment/
metadata/vegetation/
scripts/DADA2/16S/
scripts/DADA2/ITS/
analyses/soil_environment/{scripts,inputs,outputs}/
analyses/vegetation/{scripts,inputs,outputs}/
analyses/alpha_diversity/{scripts,inputs,outputs}/
analyses/beta_diversity/{scripts,inputs,outputs}/
analyses/differential_abundance/{scripts,inputs,outputs}/
analyses/community_assembly/{scripts,inputs,outputs}/
data/processed/
```
Use a single shared copy of inputs where practical and link to them from analysis READMEs.

## Audit checks
1. Obtain actual source-file inventories and export scripts/outputs from Mac, Gandalf, and Gadi.
2. For each candidate, record source system, source path, size, checksum, analysis group, and manuscript alignment.
3. Flag obsolete DESeq2, duplicate scripts, placeholders, absolute machine-specific paths, and large/raw sequencing files.
4. Stage only validated files, verify script input references, and document software versions.
5. Review a GitHub diff, confirm DADA2 and Krona unchanged, then merge to main.

**Migration state:** Repository surveyed and branch created. Sources from local computer, Gandalf and Gadi have not yet been imported or verified; no old analysis files have been deleted. Never interpret this document as indicating that the new analyses are already published.

Raw reads: NCBI SRA BioProject PRJNA1446205.
