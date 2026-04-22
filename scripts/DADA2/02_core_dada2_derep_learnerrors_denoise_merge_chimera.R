library(dada2)

filt_path <- "~/Amplicon_sequencing/16S_contam_assessment/filtered"
obj_path  <- "~/Amplicon_sequencing/16S_contam_assessment/objects"
res_path  <- "~/Amplicon_sequencing/16S_contam_assessment/results"

dir.create(obj_path, showWarnings = FALSE, recursive = TRUE)
dir.create(res_path, showWarnings = FALSE, recursive = TRUE)

filtFs <- sort(list.files(filt_path, pattern = "_R1\\.filtered\\.fastq\\.gz$", full.names = TRUE))
filtRs <- sort(list.files(filt_path, pattern = "_R2\\.filtered\\.fastq\\.gz$", full.names = TRUE))
sample.names <- sub("_R1\\.filtered\\.fastq\\.gz$", "", basename(filtFs))

cat("Filtered forward files:", length(filtFs), "\n")
cat("Filtered reverse files:", length(filtRs), "\n")
cat("Control present:", any(grepl("CONTROL2", sample.names)), "\n")

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

saveRDS(derepFs, file.path(obj_path, "derepFs.rds"))
saveRDS(derepRs, file.path(obj_path, "derepRs.rds"))

errF <- learnErrors(filtFs, multithread = TRUE, verbose = 1)
errR <- learnErrors(filtRs, multithread = TRUE, verbose = 1)
saveRDS(errF, file.path(obj_path, "errF.rds"))
saveRDS(errR, file.path(obj_path, "errR.rds"))

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
saveRDS(dadaFs, file.path(obj_path, "dadaFs.rds"))
saveRDS(dadaRs, file.path(obj_path, "dadaRs.rds"))

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
saveRDS(mergers, file.path(obj_path, "mergers.rds"))

seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, file.path(obj_path, "seqtab.rds"))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
saveRDS(seqtab.nochim, file.path(obj_path, "seqtab_nochim.rds"))

track <- data.frame(
  sample = sample.names,
  filtered = rowSums(seqtab),
  nonchim = rowSums(seqtab.nochim),
  ASVs_postmerge = rowSums(seqtab > 0),
  ASVs_nonchim = rowSums(seqtab.nochim > 0)
)

write.csv(track, file.path(res_path, "track_core_dada2.csv"), row.names = FALSE)

cat("\nCONTROL summary after core DADA2:\n")
print(track[grep("CONTROL2", track$sample), ])

control_idx <- grep("CONTROL2", rownames(seqtab.nochim), ignore.case = TRUE)
if (length(control_idx) > 0) {
  control_counts <- seqtab.nochim[control_idx, ]
  control_counts <- sort(control_counts[control_counts > 0], decreasing = TRUE)
  control_df <- data.frame(
    ASV = names(control_counts),
    Count = as.numeric(control_counts)
  )
  write.csv(control_df, file.path(res_path, "control2_asvs_nonchim.csv"), row.names = FALSE)
  cat("\nTop ASVs in CONTROL2 after chimera removal:\n")
  print(head(control_df, 20))
} else {
  cat("\nCONTROL2 not found in seqtab.nochim rownames\n")
}
