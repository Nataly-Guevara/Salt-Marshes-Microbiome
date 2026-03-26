# Set custom library path
.libPaths("/scratch/kl74/ng1317/R/linux_clean")
library(dada2)

# Define paths
filt_path <- "/scratch/kl74/ng1317/16S/filtered"
output_path <- "/scratch/kl74/ng1317/16S/denoised_batch1"

# Create output directory if it doesn't exist
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# Load dereplicated reads
derepFs <- readRDS(file.path(filt_path, "derepFs.rds"))
derepRs <- readRDS(file.path(filt_path, "derepRs.rds"))

# Load learned error rates
errF <- readRDS(file.path(filt_path, "errF_subsampled.rds"))
errR <- readRDS(file.path(filt_path, "errR_subsampled.rds"))

# Subset first 20 samples
subset_names <- names(derepFs)[1:20]
derepFs <- derepFs[subset_names]
derepRs <- derepRs[subset_names]

# Run dada denoising
dadaFs <- dada(derepFs, err = errF, multithread = 4)
dadaRs <- dada(derepRs, err = errR, multithread = 4)

# Save results
saveRDS(dadaFs, file = file.path(output_path, "dadaFs.rds"))
saveRDS(dadaRs, file = file.path(output_path, "dadaRs.rds"))

