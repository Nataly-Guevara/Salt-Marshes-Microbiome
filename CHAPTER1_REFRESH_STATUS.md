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

## Privacy check required before any new public upload

- **BLOCKING gate:** no scripts, tables, ZIPs or results from the Mac, Gandalf or Gadi may be published until inspected for usernames, personal directories, absolute computer/HPC paths, contact details, credentials, private keys, hidden/local files and accidental private manuscript content.
- The existing public GitHub repository already contains hard-coded home/HPC machine paths in **11 of 16 existing R scripts** under the old script layout, including references to HPC account identifiers. This finding is from a preliminary pattern search of R scripts; it is not a full security audit and is not evidence of exposed passwords.
- The user requests removal of personal identifiers from any publication-ready scripts. The preserved scientific DADA2 steps should not be changed; configuration paths can later be made relative or environment-driven while the computational steps are retained. No DADA2 code has yet been edited.
- `.gitignore` excludes several common credential, workspace, raw sequencing and legacy-file patterns, but it **cannot prevent hardcoded usernames/passwords inside source files**, nor does it remove previously tracked files.
- The review ZIP from the user's Mac must be supplied privately in chat before its *contents* can be scanned. The local read-only `chapter1_privacy_preflight.py` tool can flag suspect lines but does not edit anything or prove completeness.
- After code inspection, publish only reviewed/sanitized, final-publication-aligned source files, and **verify the GitHub diff** before merging the refresh branch into `main`.
- Removing sensitive text in a future commit will not remove it from Git history; a separate history-remediation procedure may be required if genuinely sensitive secrets are confirmed.

## Local script ZIP privacy review (21 October 2026)
- Received and inspected the 18-script private local review ZIP (21 text members including its manifests).
- Automated preflight found **15 personal working-directory references in four scripts**, all associated with `setwd()` configuration. No password, private-key, token, email-address or HPC-user-ID matches appeared in the uploaded ZIP with the implemented checks.
- Prepared a **private sanitised review ZIP** replacing those 15 `setwd()` blocks with optional `CHAPTER1_WORKDIR` environment-variable configuration. R computational/statistical statements were not deliberately changed. Original private ZIP remains unchanged.
- Re-ran automated preflight on the cleaned ZIP: **zero pattern matches** across 22 inspected text members. This result is not a complete guarantee of privacy or reproducibility.
- **No local R source script has yet been committed to this public repository.** Sanitised review ZIP stays private until classification, human review, script-dependency checks and manuscript alignment are completed.
- Several Model A/B/C scripts contain **multiple analysis types** (alpha diversity, beta diversity, and ANCOM-BC2). Avoid assigning each whole mixed script to a single method folder without documenting the cross-method contents or extracting executable, tested modules.
- The local bundle did not identify a standalone community-assembly/iCAMP source script; find its final workflow on Gandalf or Gadi before filling that analysis folder.

## Publication and cleanup rules
- Inspect files for private material, machine-specific paths, credentials and publication-readiness before uploading to a **public** repository.
- Do not upload manuscript drafts, reviewer replies, photos, recordings, old R workspace images or superseded tests.
- Only remove old DESeq2 scripts/results and other obsolete analyses from this staging branch **after** verified replacement files are ready; review the diff and confirm DADA2/Krona unchanged.
- Keep a recoverable Git history and merge into `main` only after validation.

**Current status:** The local inventory has been classified and the branch audit has been updated. No Mac/Gandalf/Gadi source files have been uploaded to GitHub and no previous analysis files have been deleted. Raw reads: NCBI SRA BioProject PRJNA1446205.
