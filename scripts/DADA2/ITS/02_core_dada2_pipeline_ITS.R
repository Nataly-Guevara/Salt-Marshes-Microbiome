# 8. Learn Errors
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
print("✅ Error learning done")
plotErrors(errF, nominalQ = TRUE)

# 9. Dereplicate
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# 10. Denoise
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
print("✅ Denoising done")

# 11. Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

#how many merged reads per sample you got — we’ll compare that to your input read counts.
merge.stats <- sapply(mergers, function(x) sum(x$abundance))
merge.summary <- summary(merge.stats)
print(merge.summary)


# 12. Build ASV Table
seqtab <- makeSequenceTable(mergers)
write.csv(as.data.frame(dim(seqtab)), "seqtab_dims.csv")
print(table(nchar(getSequences(seqtab))))



# 13. Remove Chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
write.csv(as.data.frame(seqtab.nochim), "ASV_table_nochim.csv")
print("✅ Chimera removal done")

#Inspect your ASV table
## Total number of ASVs
ncol(seqtab.nochim)

# Total number of samples
nrow(seqtab.nochim)

# Number of reads per sample (row sums)
sample_reads <- rowSums(seqtab.nochim)
summary(sample_reads)

# Number of ASVs per sample (non-zero entries)
asvs_per_sample <- rowSums(seqtab.nochim > 0)
summary(asvs_per_sample)

# Histogram of ASVs per sample
hist(asvs_per_sample, main = "ASVs per sample", xlab = "Number of ASVs")

# Optional: ASV abundance distribution
asv_abundance <- colSums(seqtab.nochim)
hist(log10(asv_abundance + 1), main = "Log10 ASV Abundance", xlab = "log10 total reads per ASV")


# 14. Save R objects
saveRDS(seqtab.nochim, "ITS_seqtab_nochim.rds")
saveRDS(sample.names, "sample_names.rds")
