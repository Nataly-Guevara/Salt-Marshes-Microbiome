# Salt-Marshes-Microbiome

Reproducibility repository for **Chapter 1: saltmarsh restoration microbial communities**. This project examines bacterial/archaeal (16S rRNA gene) and fungal (ITS) communities across natural reference, degraded and transplanted saltmarsh conditions, including BESE and BARE treatments, and their bulk-soil and rhizosphere compartments.

> **Refresh in progress.** This `chapter1-2026-refresh` branch establishes the new Chapter 1 repository structure. The analysis folders below are currently placeholders; final source scripts, validated inputs and selected outputs have **not** yet been uploaded. Some older analysis files remain in the previous `scripts/` and `tables/` directories pending verified replacement. Those older files should not be assumed to reflect the final publication.

## Current and planned repository structure

```text
Salt-Marshes-Microbiome/
├── README.md
├── scripts/
│   └── DADA2/
│       ├── 16S/                        # preserved original pipeline
│       └── ITS/                        # preserved original pipeline
├── docs/
│   └── krona/
│       └── index.html                  # preserved visualisation
├── metadata/
│   ├── dada2/                         # preserved
│   ├── soil_parameters/               # existing
│   └── vegetation/                    # existing
├── analyses/
│   ├── soil_environment/
│   ├── vegetation/
│   ├── alpha_diversity/
│   ├── beta_diversity/
│   ├── differential_abundance/
│   └── community_assembly/
└── data/
    └── processed/
```

Each `analyses/` folder has a README describing its purpose. The intended pattern is to provide **scripts**, **input references**, and **selected outputs** once the files have been verified.

| Analysis | Intended coverage | Status |
| --- | --- | --- |
| [Soil environment](analyses/soil_environment/) | Soil physicochemistry and ecological condition | Scaffold |
| [Vegetation](analyses/vegetation/) | Plant cover and related measurements | Scaffold |
| [Alpha diversity](analyses/alpha_diversity/) | 16S and ITS alpha-diversity models and checks | Scaffold |
| [Beta diversity](analyses/beta_diversity/) | Community composition, ordination, PERMANOVA, dispersion | Scaffold |
| [Differential abundance](analyses/differential_abundance/) | Final ANCOM-BC2 workflow, replacing superseded DESeq2 analyses | Scaffold |
| [Community assembly](analyses/community_assembly/) | Manuscript-aligned community assembly analyses | Scaffold |
| [Processed data](data/processed/) | Validated and appropriately sized shared inputs | Scaffold |

## Preserved Chapter 1 components

- **DADA2:** [16S scripts](scripts/DADA2/16S/) and [ITS scripts](scripts/DADA2/ITS/) for amplicon sequence processing and phyloseq preparation.
- **Krona:** [Interactive visualisation](docs/krona/index.html).
- **Metadata:** [Marker metadata](metadata/dada2/), [soil parameters](metadata/soil_parameters/) and [vegetation](metadata/vegetation/).

These existing files were deliberately left unchanged during the structural refresh. Some DADA2 processing intermediates referenced by the scripts are not yet present in the repository and must be supplied or regenerated for full reproducibility; see the migration checklist.

## Sequence reads

Raw sequencing reads are available through NCBI SRA BioProject **[PRJNA1446205](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1446205)** and are not stored here.

## Refresh and provenance

- The source material is being reconciled from the local Chapter 1 analysis folders, Gandalf and NCI Gadi.
- The entire local `legacy Microbial community analysis/` directory contains obsolete analyses not used in the final Chapter 1 publication and is **excluded** from the new release.
- Final workflow selection will follow the published or accepted manuscript, rather than selecting files simply because they are newer.
- Any previous GitHub analyses will be removed from this branch **only after** their replacements have been verified. Git history will preserve older revisions.

See [CHAPTER1_REFRESH_STATUS.md](CHAPTER1_REFRESH_STATUS.md) for the detailed migration audit, missing inputs and cleanup checklist.
