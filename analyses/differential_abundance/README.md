# Differential abundance

Final 16S and ITS differential-abundance models, currently identified as ANCOM-BC2 (replacing older DESeq2 scripts/results after provenance validation).

## Intended contents

- `scripts/`: checked, publication-aligned source R scripts.
- `inputs/`: manifest of required input files and links to shared processed data.
- `outputs/`: selected analysis results, summaries and figures that support the final manuscript.

## Current status

**Code review only.** Privacy-checked candidate scripts have been added within this analysis folder and/or in [`shared_workflows`](../shared_workflows/), but required inputs, outputs, and manuscript-aligned results remain to be verified. Candidate files will be matched with the local inventory, Gandalf and NCI Gadi, and verified before public release. Do not cite this folder as a complete reproducible analysis until its contents are validated.

See the [Chapter 1 migration checklist](../../CHAPTER1_REFRESH_STATUS.md) for progress and exclusions.

## Hypothesis alignment

- **H1:** Between-condition contrasts in bulk soil.
- **H2:** Between-condition contrasts in rhizosphere.
- **H3:** Within-condition comparisons between paired bulk soil and rhizosphere.

The public replacement scripts and outputs remain pending validation; old DESeq2 files elsewhere on the branch have not yet been removed.

The H1/H2/H3 multi-method source files are staged under [`analyses/shared_workflows/`](../shared_workflows/). Their ANCOM-BC2, diversity and sensitivity sections have **not** been extracted or rerun.
