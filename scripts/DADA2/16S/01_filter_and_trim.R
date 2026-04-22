library(dada2)

path <- "~/Amplicon_sequencing/16S_contam_assessment/raw"
filt_path <- "~/Amplicon_sequencing/16S_contam_assessment/filtered"
dir.create(filt_path, showWarnings = FALSE, recursive = TRUE)

fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz$", full.names = TRUE))
sample.names <- sub("_R1.fastq.gz$", "", basename(fnFs))

cat("Forward files:", length(fnFs), "\n")
cat("Reverse files:", length(fnRs), "\n")
cat("Control present:", any(grepl("CONTROL2", sample.names)), "\n")

filtFs <- file.path(filt_path, paste0(sample.names, "_R1.filtered.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R2.filtered.fastq.gz"))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(240, 200),
  maxN = 0,
  maxEE = c(2, 5),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

out_df <- as.data.frame(out)
out_df$Sample <- rownames(out_df)
out_df <- out_df[, c("Sample", "reads.in", "reads.out")]

write.csv(
  out_df,
  "~/Amplicon_sequencing/16S_contam_assessment/results/filtering_summary_with_control.csv",
  row.names = FALSE
)

cat("Control summary:\n")
print(out_df[grep("CONTROL2", out_df$Sample), ])
