# Chapter 1 refresh — staging audit (21 October 2026)

**Branch:** `chapter1-2026-refresh`. The `main` branch has not been changed.

## Preserve exactly as they are
- `scripts/DADA2/16S/` and `scripts/DADA2/ITS/`
- `metadata/dada2/`
- `docs/krona/index.html`

## Local inventory findings — filenames only
The provided local Chapter 1 CSV lists **734 files (~1.29 GiB)**. It does *not* contain file contents, so no script or result is yet verified.

| Group / source directory | Inventory entries | Next step |
| --- | ---: | --- |
| Current, non-legacy R scripts | 18 | Collect files and inspect actual code, absolute paths and dependencies |
| `Microbial community analysis/outputs/AMBCOM/` | 201 | ANCOM-BC2 Models A/B/C and sensitivity results; identify producing scripts |
| `Microbial community analysis/outputs/Beta/` | 73 | Review final model outputs |
| `Microbial community analysis/outputs/Alpha/` | 63 | Review final model outputs |
| `soil_parameters/` | 74 | Verify soil analysis scripts, input version and output tables |
| `veg_analysis/` | 19 | Verify vegetation analysis scripts and outputs |
| `environmental_association_outputs/` | 23 | Locate producing script and inspect results |
| `legacy Microbial community analysis/` | 108 | **Confirmed obsolete and not used for final Chapter 1 publication. Exclude entire directory from new release; retain locally.** |
| Large `.RData` workspaces | 7 | Exclude from public GitHub (~1.06 GiB total) |

Numbers refer to file inventory locations (some analytical topics overlap) rather than independently verified analyses.

## Confirmed exclusion (user clarification)
The entire local directory `legacy Microbial community analysis/` contains only superseded analyses **not used in the final Chapter 1 publication**. Its 108 inventoried files are fully outside the migration scope: do not collect, review for selection, copy, or upload them. Keep the local originals unchanged. This exclusion does not by itself authorize removal of files already on GitHub; obsolete tracked files will only be removed from the staging branch once the new final workflows are verified.

## Highest-priority validation
1. Replace superseded **DESeq2** analysis with the final **ANCOM-BC2** workflow, after confirming which R script generated the `AMBCOM` outputs. Do not assume that scripts named Model A/B/C necessarily contain complete ANCOM-BC2 code.
2. Replace obsolete GitHub alpha/beta diversity, sensitivity and PERMANOVA work with verified final scripts, inputs and selected publication-aligned outputs. Do **not** source these replacements from the excluded local `legacy Microbial community analysis/` directory.
3. Compare locally listed corrected-metadata phyloseq objects against GitHub objects and the final Chapter 1 manuscript before replacing any shared inputs.
4. Collect missing or more authoritative final script versions from Gandalf and Gadi and reconcile local copies.
5. Document inputs, package versions, execution commands, and output provenance before changing the repository README.

## Suggested final repository layout
```text
README.md
scripts/DADA2/{16S,ITS}/            # preserve
metadata/dada2/                      # preserve
docs/krona/                          # preserve
analyses/soil_environment/{scripts,inputs,outputs}/
analyses/vegetation/{scripts,inputs,outputs}/
analyses/alpha_diversity/{scripts,inputs,outputs}/
analyses/beta_diversity/{scripts,inputs,outputs}/
analyses/differential_abundance/{scripts,inputs,outputs}/
analyses/environment_associations/{scripts,inputs,outputs}/
analyses/community_assembly/{scripts,inputs,outputs}/
analyses/taxonomy/{scripts,inputs,outputs}/
data/processed/
```


## Confirmed pre-existing reproducibility gaps in `main`

- `tables/ASV/16S_ASV_counts.csv`, `tables/ASV/ITS_ASV_counts.csv` and `tables/taxonomy/ITS_taxonomy.csv` each consist only of a CRLF newline (2 bytes), rather than an exported table. Replace placeholders with validated data in the eventual release.
- Preserved 16S phyloseq construction script references `tables/ASV/seqtab_filtered_then_nochim.rds` and `tables/taxonomy/taxa_silva138.1.rds`; those files are not currently in the repository.
- Preserved ITS phyloseq construction script references `tables/ASV/ITS_seqtab_nochim.rds` and `tables/taxonomy/ITS_taxa_assignments.rds`; those files are not currently in the repository.
- The README states phyloseq objects are under `metadata/dada2/`, whereas the currently checked-in objects are under `tables/phyloseq/`. Reconcile documentation with verified final input/output layout.

Do not rewrite the DADA2 scripts as part of this migration unless separately approved. First recover the necessary processing products, or document how to regenerate them from raw reads.


## Folder scaffold progress (21 October 2026)
All six requested `analyses/` directories now exist on the staging branch. Each has its own `README.md` plus empty `scripts/`, `inputs/` and `outputs/` subfolders represented by `.gitkeep` files. `data/processed/README.md` also exists. No final analysis source files or results are yet in those folders.

**Next action:** Collect files from the Mac using the private review tool (`collect_chapter1_review_bundles.py`) and the updated 734-row audit. Collect scripts first (18 locally inventoried candidates), then validated shared inputs and selected outputs. The tool excludes the entire `legacy Microbial community analysis` directory and preserves source-relative paths rather than blindly adopting tentative GitHub target paths. Investigate missing final ANCOM-BC2 and community assembly scripts on Gandalf/Gadi if not represented among the local R files. Do not push an unreviewed ZIP directly to this public repository.

## Publication and cleanup rules
- Inspect files for private material, machine-specific paths, credentials and publication-readiness before uploading to a **public** repository.
- Do not upload manuscript drafts, reviewer replies, photos, recordings, old R workspace images or superseded tests.
- Only remove old DESeq2 scripts/results and other obsolete analyses from this staging branch **after** verified replacement files are ready; review the diff and confirm DADA2/Krona unchanged.
- Keep a recoverable Git history and merge into `main` only after validation.

**Current status:** The local inventory has been classified and the branch audit has been updated. No Mac/Gandalf/Gadi source files have been uploaded to GitHub and no previous analysis files have been deleted. Raw reads: NCBI SRA BioProject PRJNA1446205.
