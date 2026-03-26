# Set custom library path
.libPaths("/scratch/kl74/ng1317/R/linux_clean")
library(dada2)

# Define paths
filt_path <- "/scratch/kl74/ng1317/16S/filtered"
out_path <- "/scratch/kl74/ng1317/16S/denoised_batch3"

# Create output dir if not exists
if (!dir.exists(out_path)) dir.create(out_path)

# Load dereplicated reads and error models
derepFs <- readRDS(file.path(filt_path, "derepFs.rds"))
derepRs <- readRDS(file.path(filt_path, "derepRs.rds"))
errF <- readRDS(file.path(filt_path, "errF_subsampled.rds"))
errR <- readRDS(file.path(filt_path, "errR_subsampled.rds"))

# Select 3rd batch (samples 41–60)
batch_samples <- names(derepFs)[41:60]

# Subset derep objects
derepFs_batch <- derepFs[batch_samples]
derepRs_batch <- derepRs[batch_samples]

# Denoise
dadaFs_batch <- dada(derepFs_batch, err = errF, multithread = 6)
dadaRs_batch <- dada(derepRs_batch, err = errR, multithread = 6)

# Save outputs
saveRDS(dadaFs_batch, file = file.path(out_path, "dadaFs_batch3.rds"))
saveRDS(dadaRs_batch, file = file.path(out_path, "dadaRs_batch3.rds"))

