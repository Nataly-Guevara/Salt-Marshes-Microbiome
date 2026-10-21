# Differential abundance

Final 16S and ITS differential-abundance models, currently identified as ANCOM-BC2 (replacing older DESeq2 scripts/results after provenance validation).

## Intended contents

- `scripts/`: checked, publication-aligned source R scripts.
- `inputs/`: manifest of required input files and links to shared processed data.
- `outputs/`: selected analysis results, summaries and figures that support the final manuscript.

## Current status

**Scaffold only.** No final analysis scripts, inputs, or outputs have been imported to this directory yet. Candidate files will be matched with the local inventory, Gandalf and NCI Gadi, and verified before public release. Do not cite this folder as a complete reproducible analysis until its contents are validated.

See the [Chapter 1 migration checklist](../../CHAPTER1_REFRESH_STATUS.md) for progress and exclusions.

## Hypothesis alignment

- **H1:** Between-condition contrasts in bulk soil.
- **H2:** Between-condition contrasts in rhizosphere.
- **H3:** Within-condition comparisons between paired bulk soil and rhizosphere.

The public replacement scripts and outputs remain pending validation; old DESeq2 files elsewhere on the branch have not yet been removed.
