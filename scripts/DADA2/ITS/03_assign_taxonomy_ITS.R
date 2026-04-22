# ---------------------------------------------
# TAXONOMIC CLASSIFICATION - 17 April 2025
# ITS dataset using UNITE reference database
# ---------------------------------------------

# Load dada2 package (install if not yet available)
if (!requireNamespace("dada2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("dada2")
}
library(dada2)

# Set working directory to where your filtered sequences and UNITE fasta are located
setwd("~/Library/CloudStorage/OneDrive-DeakinUniversity/Desktop/ALL DESKTOP/ITS_full/filtered")

# Load your cleaned, non-chimeric sequence table
seqtab.nochim <- readRDS("ITS_seqtab_nochim.rds")

# Path to UNITE taxonomy reference
fasta_path <- "sh_general_release_dynamic_19.02.2025.fasta"

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, fasta_path, multithread = TRUE)
print("✅ Taxonomy assignment complete")

# Save as RDS only. Tried with .csv and did not worked
saveRDS(taxa, "ITS_taxa_assignments.rds")
