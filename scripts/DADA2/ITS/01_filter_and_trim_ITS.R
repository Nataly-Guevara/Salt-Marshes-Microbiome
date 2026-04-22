# ITS_full_pipeline.R
#16-March-2025
#N.Guevara

# 1. Set working directory
setwd("~/Desktop/ALL DESKTOP/ITS_full")

# 2. Load DADA2
if (!requireNamespace("dada2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("dada2")
}
library(dada2)
packageVersion("dada2")

# 3. Set filtered output directory (creates 'filtered' subfolder inside ITS_full)
filt_path <- file.path(getwd(), "filtered")
if (!dir.exists(filt_path)) dir.create(filt_path)

# 4. List all FASTQ files fnFs: Foward, fnRs: Reverse 
fnFs <- sort(list.files(pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(pattern = "_R2.fastq.gz", full.names = TRUE))

# 5. Extract sample names from filenames
get_sample_name <- function(filename) {
  parts <- strsplit(basename(filename), "_")[[1]]
  paste(parts[1:5], collapse = "_")
}
sample.names <- sapply(fnFs, get_sample_name)

# 6. Set output paths for filtered reads
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))


# 7. Filter and Trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(250, 200),
                     maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     compress = TRUE, multithread = TRUE)
write.csv(out, "filtering_stats.csv")
print("✅ Filtering done")
getwd()

#Assess Filtering Efficiency
#You can calculate the percentage of reads retained with this:
out_df <- read.csv("filtered/output/filtering_stats.csv")  # or the correct path
out_df$percent_retained <- round((out_df$reads.out / out_df$reads.in) * 100, 2)

summary(out_df$percent_retained)  # Quick stats
