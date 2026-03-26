#!/usr/bin/env Rscript
.libPaths("/scratch/kl74/ng1317/R/linux_clean")
library(dada2)

den <- "/scratch/kl74/ng1317/16S/denoising_all"

# load the chimera-free table (trimmed to ~25 k ASVs)
seqtab <- readRDS(file.path(den, "seqtab_nochim.rds"))

# SILVA v138.1 training set (gzip stays compressed on disk)
silva_train <- "/g/data/kl74/databases/silva_nr99_v138.1_train_set.fa.gz"

# assign taxonomy (16 threads = 1 thread per CPU you request)
taxa <- assignTaxonomy(seqtab,
                       silva_train,
                       multithread = 16)

saveRDS(taxa, file.path(den, "taxa_silva.rds"))
cat("✅  Taxonomy saved  ➜  taxa_silva.rds\n")

