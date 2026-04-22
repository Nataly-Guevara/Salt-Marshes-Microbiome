#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(phyloseq)
})

# Paths
base_path <- "~/Library/CloudStorage/OneDrive-DeakinUniversity/ALL DESKTOP/ITS_full/filtered"
seqtab_file <- file.path(base_path, "ITS_seqtab_nochim.rds")
taxa_file   <- file.path(base_path, "ITS_taxa_assignments.rds")
meta_file   <- file.path(base_path, "metadata.csv")

out_dir <- file.path(base_path, "Phyloseq")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load inputs
seqtab <- readRDS(seqtab_file)
taxa   <- readRDS(taxa_file)
taxa   <- as.matrix(taxa)

# Standardize taxonomy rank names
expected_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
if (is.null(colnames(taxa))) {
  colnames(taxa) <- expected_ranks[seq_len(ncol(taxa))]
}
if (ncol(taxa) < length(expected_ranks)) {
  taxa <- cbind(taxa, matrix(NA, nrow = nrow(taxa), ncol = length(expected_ranks) - ncol(taxa)))
  colnames(taxa) <- expected_ranks
}
if (ncol(taxa) > length(expected_ranks)) {
  taxa <- taxa[, seq_along(expected_ranks), drop = FALSE]
  colnames(taxa) <- expected_ranks
}

# Match ASVs
common_asvs <- intersect(colnames(seqtab), rownames(taxa))
seqtab <- seqtab[, common_asvs, drop = FALSE]
taxa   <- taxa[common_asvs, , drop = FALSE]

# Build phyloseq
ps <- phyloseq(
  otu_table(seqtab, taxa_are_rows = FALSE),
  tax_table(taxa)
)

# Add metadata
if (file.exists(meta_file)) {
  metadata <- read.csv(meta_file, row.names = 1, check.names = FALSE)
  common_samples <- intersect(sample_names(ps), rownames(metadata))
  if (length(common_samples) > 0) {
    ps <- prune_samples(common_samples, ps)
    metadata <- metadata[common_samples, , drop = FALSE]
    ps <- merge_phyloseq(ps, sample_data(metadata))
  } else {
    warning("No overlapping sample names between metadata and phyloseq object")
  }
} else {
  warning("Metadata file not found: ", meta_file)
}

# Remove likely controls
control_samples <- grep("control|blank|ntc|mock", sample_names(ps), ignore.case = TRUE, value = TRUE)
if (length(control_samples) > 0) {
  ps <- prune_samples(setdiff(sample_names(ps), control_samples), ps)
}

# Keep fungi only
taxmat <- as.data.frame(as.matrix(tax_table(ps)), stringsAsFactors = FALSE)
keep_fungi <- taxmat$Kingdom == "Fungi"
keep_fungi[is.na(keep_fungi)] <- FALSE
ps <- prune_taxa(keep_fungi, ps)

# Remove zero-count taxa/samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps <- prune_samples(sample_sums(ps) > 0, ps)

# Save outputs
saveRDS(ps, file.path(out_dir, "phyloseq_ITS_clean.rds"))
write.csv(as.data.frame(otu_table(ps)), file.path(out_dir, "ASV_table_ITS_clean.csv"))
write.csv(as.data.frame(as.matrix(tax_table(ps))), file.path(out_dir, "taxonomy_ITS_clean.csv"))

cat("Saved cleaned ITS phyloseq object\n")
cat("Samples:", nsamples(ps), "\n")
cat("Taxa:", ntaxa(ps), "\n")