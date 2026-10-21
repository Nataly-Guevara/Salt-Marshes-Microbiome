# Beta diversity

16S and ITS community composition, ordinations, distance-based comparisons, PERMANOVA, dispersion checks and supplementary analyses from the final Chapter 1 workflow.

## Intended contents

- `scripts/`: checked, publication-aligned source R scripts.
- `inputs/`: manifest of required input files and links to shared processed data.
- `outputs/`: selected analysis results, summaries and figures that support the final manuscript.

## Current status

**Code review only.** Privacy-checked candidate scripts have been added within this analysis folder and/or in [`shared_workflows`](../shared_workflows/), but required inputs, outputs, and manuscript-aligned results remain to be verified. Candidate files will be matched with the local inventory, Gandalf and NCI Gadi, and verified before public release. Do not cite this folder as a complete reproducible analysis until its contents are validated.

See the [Chapter 1 migration checklist](../../CHAPTER1_REFRESH_STATUS.md) for progress and exclusions.

## Hypothesis alignment

- **H1:** Bulk-soil community composition and reference similarity.
- **H2:** Rhizosphere community composition and reference similarity.
- **H3:** Paired bulk-soil–rhizosphere community differentiation.

The H1/H2/H3 multi-method source files are staged under [`analyses/shared_workflows/`](../shared_workflows/). Their ANCOM-BC2, diversity and sensitivity sections have **not** been extracted or rerun.
