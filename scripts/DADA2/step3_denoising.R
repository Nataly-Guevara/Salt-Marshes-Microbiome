# =========================
# step3_denoising.R (Batch)
# =========================

# Custom library path
.libPaths("/scratch/kl74/ng1317/R/linux_clean")
library(dada2)

# Paths
filt_path <- "/scratch/kl74/ng1317/16S/filtered"

# Load filtered files
filtFs <- sort(list.files(filt_path, pattern = "_R1\\.filtered\\.fastq\\.gz$", full.names = TRUE))
filtRs <- sort(list.files(filt_path, pattern = "_R2\\.filtered\\.fastq\\.gz$", full.names = TRUE))

# Load dereplication objects
derepFs <- readRDS(file.path(filt_path, "derepFs.rds"))
derepRs <- readRDS(file.path(filt_path, "derepRs.rds"))

# Load learned error models (previously generated)
errF <- readRDS(file.path(filt_path, "errF_subsampled.rds"))
errR <- readRDS(file.path(filt_path, "errR_subsampled.rds"))

# Get sample names
samples <- names(derepFs)

# Process in batches of 20
batch_size <- 20
batches <- split(samples, ceiling(seq_along(samples) / batch_size))

# Initialize output lists
dadaFs <- list()
dadaRs <- list()

for (i in seq_along(batches)) {
  cat("\nProcessing batch", i, "of", length(batches), "\n")
  batch <- batches[[i]]

  dadaFs[batch] <- dada(derepFs[batch], err = errF, multithread = TRUE)
  dadaRs[batch] <- dada(derepRs[batch], err = errR, multithread = TRUE)

  # Save intermediate output
  saveRDS(dadaFs, file = file.path(filt_path, paste0("dadaFs_batch", i, ".rds")))
  saveRDS(dadaRs, file = file.path(filt_path, paste0("dadaRs_batch", i, ".rds")))
}

# Optionally save combined outputs if all batches done
saveRDS(dadaFs, file = file.path(filt_path, "dadaFs_all.rds"))
saveRDS(dadaRs, file = file.path(filt_path, "dadaRs_all.rds"))

cat("\n✅ Denoising complete. Output saved to:", filt_path, "\n")

