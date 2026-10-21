# ============================================================
# H1 ALPHA DIVERSITY — BACTERIA
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(emmeans)
  library(car)
})

# ------------------------------------------------------------
# 1. Load Chapter 1 analysis-ready dataset
# ------------------------------------------------------------

ps_bac <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

cat("\n========================================\n")
cat("FULL BACTERIAL DATASET\n")
cat("========================================\n")
cat("Samples:", nsamples(ps_bac), "\n")
cat("ASVs:", ntaxa(ps_bac), "\n")
cat("Reads:", sum(sample_sums(ps_bac)), "\n")

stopifnot(
  nsamples(ps_bac) == 103,
  ntaxa(ps_bac) == 32530
)

# ------------------------------------------------------------
# 2. H1 — bulk soil
# ------------------------------------------------------------

ps_A <- subset_samples(
  ps_bac,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_A <- prune_taxa(
  taxa_sums(ps_A) > 0,
  ps_A
)

cat("\n========================================\n")
cat("H1 — BACTERIA\n")
cat("========================================\n")
cat("Samples:", nsamples(ps_A), "\n")
cat("ASVs detected:", ntaxa(ps_A), "\n")
cat("Reads:", sum(sample_sums(ps_A)), "\n\n")

print(
  with(
    as.data.frame(sample_data(ps_A)),
    table(Site, Treatment)
  )
)

stopifnot(nsamples(ps_A) == 58)

# ------------------------------------------------------------
# 3. Calculate Chao1 and Shannon
# ------------------------------------------------------------

rich <- as.data.frame(
  estimate_richness(
    ps_A,
    measures = c("Chao1", "Shannon")
  )
)

rich$SampleID <- sub(
  "^X(?=\\d)",
  "",
  rownames(rich),
  perl = TRUE
)

md <- methods::as(
  sample_data(ps_A),
  "data.frame"
)

# Use phyloseq sample names as the unambiguous join key
md$SampleID <- sub(
  "^X(?=\\d)",
  "",
  rownames(md),
  perl = TRUE
)

alpha_A <- merge(
  rich,
  md,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

alpha_A$Site <- factor(
  alpha_A$Site,
  levels = c("ACP", "BASS", "TRUG")
)

alpha_A$Treatment <- factor(
  alpha_A$Treatment,
  levels = c(
    "NATURAL",
    "BESE",
    "BARE",
    "DEGRADED"
  )
)

cat("\nAlpha rows:", nrow(alpha_A), "\n")
cat("Missing Site:", sum(is.na(alpha_A$Site)), "\n")
cat("Missing Treatment:", sum(is.na(alpha_A$Treatment)), "\n")

stopifnot(
  nrow(alpha_A) == 58,
  !anyNA(alpha_A$Site),
  !anyNA(alpha_A$Treatment)
)

# ------------------------------------------------------------
# 4. CHAO1
# ------------------------------------------------------------

cat("\n========================================\n")
cat("CHAO1 — FULL MODEL\n")
cat("========================================\n")

mod_chao <- lm(
  Chao1 ~ Site * Treatment,
  data = alpha_A
)

chao_anova <- car::Anova(
  mod_chao,
  type = 2
)

print(chao_anova)

cat("\nCHAO1 — Treatment within each Site\n")

chao_site <- emmeans::joint_tests(
  mod_chao,
  by = "Site"
)

print(chao_site)

chao_site_df <- as.data.frame(chao_site)

sig_chao_sites <- as.character(
  chao_site_df$Site[
    chao_site_df$p.value < 0.05
  ]
)

if (length(sig_chao_sites) > 0) {
  
  for (s in sig_chao_sites) {
    
    cat(
      "\nCHAO1 Tukey comparisons —",
      s,
      "\n"
    )
    
    emm <- emmeans(
      mod_chao,
      ~ Treatment | Site,
      at = list(Site = s)
    )
    
    print(
      contrast(
        emm,
        method = "pairwise",
        adjust = "tukey"
      )
    )
  }
}

# ------------------------------------------------------------
# 5. SHANNON
# ------------------------------------------------------------

cat("\n========================================\n")
cat("SHANNON — FULL MODEL\n")
cat("========================================\n")

mod_shannon <- lm(
  Shannon ~ Site * Treatment,
  data = alpha_A
)

shannon_anova <- car::Anova(
  mod_shannon,
  type = 2
)

print(shannon_anova)

cat("\nSHANNON — Treatment within each Site\n")

shannon_site <- emmeans::joint_tests(
  mod_shannon,
  by = "Site"
)

print(shannon_site)

shannon_site_df <- as.data.frame(
  shannon_site
)

sig_shannon_sites <- as.character(
  shannon_site_df$Site[
    shannon_site_df$p.value < 0.05
  ]
)

if (length(sig_shannon_sites) > 0) {
  
  for (s in sig_shannon_sites) {
    
    cat(
      "\nSHANNON Tukey comparisons —",
      s,
      "\n"
    )
    
    emm <- emmeans(
      mod_shannon,
      ~ Treatment | Site,
      at = list(Site = s)
    )
    
    print(
      contrast(
        emm,
        method = "pairwise",
        adjust = "tukey"
      )
    )
  }
}

# ------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------

outdir <- "outputs/H1_alpha/bacteria"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  alpha_A,
  file.path(
    outdir,
    "H1_bacteria_alpha_values.csv"
  ),
  row.names = FALSE
)

dataset_summary <- data.frame(
  Dataset = c(
    "Chapter 1 analysis-ready 16S",
    "H1 bulk soil"
  ),
  Samples = c(
    nsamples(ps_bac),
    nsamples(ps_A)
  ),
  ASVs = c(
    ntaxa(ps_bac),
    ntaxa(ps_A)
  ),
  Reads = c(
    sum(sample_sums(ps_bac)),
    sum(sample_sums(ps_A))
  )
)

write.csv(
  dataset_summary,
  file.path(
    outdir,
    "H1_bacteria_dataset_summary.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "sessionInfo_H1_bacteria.txt"
  )
)

cat("\n========================================\n")
cat("DONE\n")
cat("========================================\n")
print(dataset_summary)







# ============================================================
# H1 ALPHA DIVERSITY — FUNGI
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(emmeans)
  library(car)
})

# ------------------------------------------------------------
# 1. Load Chapter 1 analysis-ready ITS dataset
# ------------------------------------------------------------

ps_fun <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)

cat("\n========================================\n")
cat("FULL FUNGAL DATASET\n")
cat("========================================\n")
cat("Samples:", nsamples(ps_fun), "\n")
cat("ASVs:", ntaxa(ps_fun), "\n")
cat("Reads:", sum(sample_sums(ps_fun)), "\n")

stopifnot(
  nsamples(ps_fun) == 87,
  ntaxa(ps_fun) == 11176
)

# ------------------------------------------------------------
# 2. H1 — bulk soil
# ------------------------------------------------------------

ps_A_fun <- subset_samples(
  ps_fun,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_A_fun <- prune_taxa(
  taxa_sums(ps_A_fun) > 0,
  ps_A_fun
)

cat("\n========================================\n")
cat("H1 — FUNGI\n")
cat("========================================\n")
cat("Samples:", nsamples(ps_A_fun), "\n")
cat("ASVs detected:", ntaxa(ps_A_fun), "\n")
cat("Reads:", sum(sample_sums(ps_A_fun)), "\n\n")

print(
  with(
    methods::as(sample_data(ps_A_fun), "data.frame"),
    table(Site, Treatment)
  )
)

# Expected from manuscript df:
# 46 observations -> residual df = 34 in 12 Site × Treatment cells
stopifnot(nsamples(ps_A_fun) == 46)

# ------------------------------------------------------------
# 3. Calculate Chao1 and Shannon
# ------------------------------------------------------------

rich_fun <- as.data.frame(
  estimate_richness(
    ps_A_fun,
    measures = c("Chao1", "Shannon")
  )
)

rich_fun$SampleID <- sub(
  "^X(?=\\d)",
  "",
  rownames(rich_fun),
  perl = TRUE
)

md_fun <- methods::as(
  sample_data(ps_A_fun),
  "data.frame"
)

md_fun$SampleID <- sub(
  "^X(?=\\d)",
  "",
  rownames(md_fun),
  perl = TRUE
)

alpha_fun_A <- merge(
  rich_fun,
  md_fun,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

alpha_fun_A$Site <- factor(
  alpha_fun_A$Site,
  levels = c("ACP", "BASS", "TRUG")
)

alpha_fun_A$Treatment <- factor(
  alpha_fun_A$Treatment,
  levels = c(
    "NATURAL",
    "BESE",
    "BARE",
    "DEGRADED"
  )
)

cat("\nAlpha rows:", nrow(alpha_fun_A), "\n")
cat("Missing Site:", sum(is.na(alpha_fun_A$Site)), "\n")
cat("Missing Treatment:", sum(is.na(alpha_fun_A$Treatment)), "\n")

stopifnot(
  nrow(alpha_fun_A) == 46,
  !anyNA(alpha_fun_A$Site),
  !anyNA(alpha_fun_A$Treatment)
)

# ------------------------------------------------------------
# 4. CHAO1
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FUNGI CHAO1 — FULL MODEL\n")
cat("========================================\n")

mod_chao_fun <- lm(
  Chao1 ~ Site * Treatment,
  data = alpha_fun_A
)

chao_fun_anova <- car::Anova(
  mod_chao_fun,
  type = 2
)

print(chao_fun_anova)

cat("\nFUNGI CHAO1 — Treatment within each Site\n")

chao_fun_site <- emmeans::joint_tests(
  mod_chao_fun,
  by = "Site"
)

print(chao_fun_site)

chao_fun_site_df <- as.data.frame(
  chao_fun_site
)

sig_chao_fun_sites <- as.character(
  chao_fun_site_df$Site[
    chao_fun_site_df$p.value < 0.05
  ]
)

chao_fun_pairs_all <- list()

if (length(sig_chao_fun_sites) > 0) {
  
  for (s in sig_chao_fun_sites) {
    
    cat(
      "\nFUNGI CHAO1 Tukey comparisons —",
      s,
      "\n"
    )
    
    emm <- emmeans(
      mod_chao_fun,
      ~ Treatment | Site,
      at = list(Site = s)
    )
    
    pair_df <- as.data.frame(
      contrast(
        emm,
        method = "pairwise",
        adjust = "tukey"
      )
    )
    
    print(pair_df)
    
    chao_fun_pairs_all[[s]] <- pair_df
  }
}

# ------------------------------------------------------------
# 5. SHANNON
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FUNGI SHANNON — FULL MODEL\n")
cat("========================================\n")

mod_shannon_fun <- lm(
  Shannon ~ Site * Treatment,
  data = alpha_fun_A
)

shannon_fun_anova <- car::Anova(
  mod_shannon_fun,
  type = 2
)

print(shannon_fun_anova)

cat("\nFUNGI SHANNON — Treatment within each Site\n")

shannon_fun_site <- emmeans::joint_tests(
  mod_shannon_fun,
  by = "Site"
)

print(shannon_fun_site)

shannon_fun_site_df <- as.data.frame(
  shannon_fun_site
)

sig_shannon_fun_sites <- as.character(
  shannon_fun_site_df$Site[
    shannon_fun_site_df$p.value < 0.05
  ]
)

shannon_fun_pairs_all <- list()

if (length(sig_shannon_fun_sites) > 0) {
  
  for (s in sig_shannon_fun_sites) {
    
    cat(
      "\nFUNGI SHANNON Tukey comparisons —",
      s,
      "\n"
    )
    
    emm <- emmeans(
      mod_shannon_fun,
      ~ Treatment | Site,
      at = list(Site = s)
    )
    
    pair_df <- as.data.frame(
      contrast(
        emm,
        method = "pairwise",
        adjust = "tukey"
      )
    )
    
    print(pair_df)
    
    shannon_fun_pairs_all[[s]] <- pair_df
  }
}

# ------------------------------------------------------------
# 6. Save statistical outputs
# ------------------------------------------------------------

outdir_fun <- "outputs/H1_alpha/fungi"

dir.create(
  outdir_fun,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  alpha_fun_A,
  file.path(
    outdir_fun,
    "H1_fungi_alpha_values.csv"
  ),
  row.names = FALSE
)

chao_fun_anova_df <- data.frame(
  Term = rownames(as.data.frame(chao_fun_anova)),
  as.data.frame(chao_fun_anova),
  row.names = NULL,
  check.names = FALSE
)

shannon_fun_anova_df <- data.frame(
  Term = rownames(as.data.frame(shannon_fun_anova)),
  as.data.frame(shannon_fun_anova),
  row.names = NULL,
  check.names = FALSE
)

write.csv(
  chao_fun_anova_df,
  file.path(
    outdir_fun,
    "H1_fungi_Chao1_ANOVA.csv"
  ),
  row.names = FALSE
)

write.csv(
  shannon_fun_anova_df,
  file.path(
    outdir_fun,
    "H1_fungi_Shannon_ANOVA.csv"
  ),
  row.names = FALSE
)

write.csv(
  chao_fun_site_df,
  file.path(
    outdir_fun,
    "H1_fungi_Chao1_site_tests.csv"
  ),
  row.names = FALSE
)

write.csv(
  shannon_fun_site_df,
  file.path(
    outdir_fun,
    "H1_fungi_Shannon_site_tests.csv"
  ),
  row.names = FALSE
)

if (length(chao_fun_pairs_all) > 0) {
  chao_pairs_save <- do.call(
    rbind,
    chao_fun_pairs_all
  )
  
  write.csv(
    chao_pairs_save,
    file.path(
      outdir_fun,
      "H1_fungi_Chao1_significant_site_Tukey.csv"
    ),
    row.names = FALSE
  )
}

if (length(shannon_fun_pairs_all) > 0) {
  shannon_pairs_save <- do.call(
    rbind,
    shannon_fun_pairs_all
  )
  
  write.csv(
    shannon_pairs_save,
    file.path(
      outdir_fun,
      "H1_fungi_Shannon_significant_site_Tukey.csv"
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 7. Whole-model R2 + dataset summary
# ------------------------------------------------------------

fungal_model_summary <- data.frame(
  Kingdom = "Fungi",
  Model = "A",
  Metric = c("Chao1", "Shannon"),
  n = c(
    nobs(mod_chao_fun),
    nobs(mod_shannon_fun)
  ),
  R2 = c(
    summary(mod_chao_fun)$r.squared,
    summary(mod_shannon_fun)$r.squared
  ),
  Adjusted_R2 = c(
    summary(mod_chao_fun)$adj.r.squared,
    summary(mod_shannon_fun)$adj.r.squared
  )
)

write.csv(
  fungal_model_summary,
  file.path(
    outdir_fun,
    "H1_fungi_model_R2.csv"
  ),
  row.names = FALSE
)

fungal_dataset_summary <- data.frame(
  Dataset = c(
    "Chapter 1 analysis-ready ITS",
    "H1 bulk soil"
  ),
  Samples = c(
    nsamples(ps_fun),
    nsamples(ps_A_fun)
  ),
  ASVs = c(
    ntaxa(ps_fun),
    ntaxa(ps_A_fun)
  ),
  Reads = c(
    sum(sample_sums(ps_fun)),
    sum(sample_sums(ps_A_fun))
  )
)

write.csv(
  fungal_dataset_summary,
  file.path(
    outdir_fun,
    "H1_fungi_dataset_summary.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    outdir_fun,
    "sessionInfo_H1_fungi.txt"
  )
)

cat("\n========================================\n")
cat("DONE — H1 FUNGI\n")
cat("========================================\n")

print(fungal_dataset_summary)

cat("\nModel R2:\n")
print(fungal_model_summary)


# ============================================================
# H1 — BULK-SOIL CHAO1
# FIX SAMPLE-ID ALIGNMENT + TEST ECOLOGICAL CONTRASTS
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(emmeans)
  library(car)
})

# ============================================================
# 1. LOAD VERIFIED PHYLOSEQ OBJECTS
# ============================================================

ps_bacteria <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)

# ============================================================
# 2. H1 — BULK SOIL
# ============================================================

ps_bac_A <- subset_samples(
  ps_bacteria,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_bac_A <- prune_taxa(
  taxa_sums(ps_bac_A) > 0,
  ps_bac_A
)

ps_fun_A <- subset_samples(
  ps_fungi,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_fun_A <- prune_taxa(
  taxa_sums(ps_fun_A) > 0,
  ps_fun_A
)

# ============================================================
# 3. FUNCTION TO BUILD ALPHA TABLE SAFELY
#    NO left_join()
# ============================================================

make_alpha_table <- function(ps_object) {
  
  # Richness
  rich <- phyloseq::estimate_richness(
    ps_object,
    measures = c(
      "Chao1",
      "Shannon"
    )
  )
  
  # estimate_richness() has converted numeric-leading
  # sample IDs to syntactic R names beginning with X.
  rownames(rich) <- sub(
    "^X(?=\\d)",
    "",
    rownames(rich),
    perl = TRUE
  )
  
  # Metadata directly from phyloseq
  meta <- as(
    sample_data(ps_object),
    "data.frame"
  )
  
  # ----------------------------------------------------------
  # CRITICAL CHECK:
  # both datasets must contain exactly the same samples
  # ----------------------------------------------------------
  
  stopifnot(
    setequal(
      rownames(rich),
      rownames(meta)
    )
  )
  
  # Put metadata in EXACTLY the richness-table order
  meta <- meta[
    rownames(rich),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(rich),
      rownames(meta)
    )
  )
  
  # Bind directly — no join ambiguity
  alpha <- cbind(
    rich,
    meta
  )
  
  alpha$Treatment <- factor(
    alpha$Treatment,
    levels = c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
  )
  
  alpha$Site <- factor(
    alpha$Site,
    levels = c(
      "ACP",
      "BASS",
      "TRUG"
    )
  )
  
  return(alpha)
}

# ============================================================
# 4. BUILD CORRECT TABLES
# ============================================================

alpha_bacteria <- make_alpha_table(
  ps_bac_A
)

alpha_fungi <- make_alpha_table(
  ps_fun_A
)

# ============================================================
# 5. VERIFY — THERE MUST BE NO ARTIFICIAL NA TREATMENT
# ============================================================

cat("\n========================================\n")
cat("CORRECTED BACTERIAL ALPHA TABLE\n")
cat("========================================\n")

cat("Rows:", nrow(alpha_bacteria), "\n")

print(
  table(
    alpha_bacteria$Treatment,
    useNA = "ifany"
  )
)

cat(
  "Missing Treatment:",
  sum(is.na(alpha_bacteria$Treatment)),
  "\n"
)


cat("\n========================================\n")
cat("CORRECTED FUNGAL ALPHA TABLE\n")
cat("========================================\n")

cat("Rows:", nrow(alpha_fungi), "\n")

print(
  table(
    alpha_fungi$Treatment,
    useNA = "ifany"
  )
)

cat(
  "Missing Treatment:",
  sum(is.na(alpha_fungi$Treatment)),
  "\n"
)


# ============================================================
# 6. RE-FIT VERIFIED H1 CHAO1 MODELS
# ============================================================

mod_bac_chao <- lm(
  Chao1 ~ Site * Treatment,
  data = alpha_bacteria
)

mod_fun_chao <- lm(
  Chao1 ~ Site * Treatment,
  data = alpha_fungi
)

cat("\n========================================\n")
cat("BACTERIA — SITE × CONDITION\n")
cat("========================================\n")

print(
  car::Anova(
    mod_bac_chao,
    type = 2
  )
)


cat("\n========================================\n")
cat("FUNGI — SITE × CONDITION\n")
cat("========================================\n")

print(
  car::Anova(
    mod_fun_chao,
    type = 2
  )
)


# ============================================================
# 7. EQUAL-WEIGHT MARGINAL MEANS ACROSS SITES
# ============================================================

emm_bac <- emmeans::emmeans(
  mod_bac_chao,
  ~ Treatment,
  weights = "equal"
)

emm_fun <- emmeans::emmeans(
  mod_fun_chao,
  ~ Treatment,
  weights = "equal"
)

# ============================================================
# 8. ORIGINAL PLANNED CONTRAST
# Natural vs mean(BESE + BARE + Degraded)
# ============================================================

pooled_bac <- emmeans::contrast(
  emm_bac,
  method = list(
    "Natural - mean(BESE, BARE, Degraded)" =
      c(1, -1/3, -1/3, -1/3)
  ),
  adjust = "none"
)

pooled_fun <- emmeans::contrast(
  emm_fun,
  method = list(
    "Natural - mean(BESE, BARE, Degraded)" =
      c(1, -1/3, -1/3, -1/3)
  ),
  adjust = "none"
)

cat("\n========================================\n")
cat("ORIGINAL PLANNED CONTRAST — BACTERIA\n")
cat("========================================\n")
print(pooled_bac)

cat("\n========================================\n")
cat("ORIGINAL PLANNED CONTRAST — FUNGI\n")
cat("========================================\n")
print(pooled_fun)


# ============================================================
# 9. THREE ECOLOGICALLY RELEVANT INDIVIDUAL CONTRASTS
#
# Natural vs BESE
# Natural vs BARE
# Natural vs Degraded
#
# Equal weighting across sites.
# Holm correction across these three planned comparisons.
# ============================================================

natural_contrasts <- list(
  "Natural - BESE" =
    c(1, -1, 0, 0),
  
  "Natural - BARE" =
    c(1, 0, -1, 0),
  
  "Natural - Degraded" =
    c(1, 0, 0, -1)
)

natural_pairs_bac <- emmeans::contrast(
  emm_bac,
  method = natural_contrasts,
  adjust = "holm"
)

natural_pairs_fun <- emmeans::contrast(
  emm_fun,
  method = natural_contrasts,
  adjust = "holm"
)


cat("\n========================================\n")
cat("BACTERIA — NATURAL VS EACH CONDITION\n")
cat("Equal-weight average across Sites\n")
cat("Holm adjusted\n")
cat("========================================\n")

print(
  summary(
    natural_pairs_bac,
    infer = TRUE
  )
)


cat("\n========================================\n")
cat("FUNGI — NATURAL VS EACH CONDITION\n")
cat("Equal-weight average across Sites\n")
cat("Holm adjusted\n")
cat("========================================\n")

print(
  summary(
    natural_pairs_fun,
    infer = TRUE
  )
)


# ============================================================
# 10. SAVE RESULTS
# ============================================================

dir.create(
  "outputs/H1_alpha/contrasts",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  as.data.frame(pooled_bac),
  "outputs/H1_alpha/contrasts/H1_bacteria_Chao1_Natural_vs_nonreference_mean.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(pooled_fun),
  "outputs/H1_alpha/contrasts/H1_fungi_Chao1_Natural_vs_nonreference_mean.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(
    summary(
      natural_pairs_bac,
      infer = TRUE
    )
  ),
  "outputs/H1_alpha/contrasts/H1_bacteria_Chao1_Natural_pairwise.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(
    summary(
      natural_pairs_fun,
      infer = TRUE
    )
  ),
  "outputs/H1_alpha/contrasts/H1_fungi_Chao1_Natural_pairwise.csv",
  row.names = FALSE
)

saveRDS(
  mod_bac_chao,
  "outputs/H1_alpha/contrasts/H1_bacteria_Chao1_model.rds"
)

saveRDS(
  mod_fun_chao,
  "outputs/H1_alpha/contrasts/H1_fungi_Chao1_model.rds"
)


# ============================================================
# H1 — BULK-SOIL BETA DIVERSITY
# PERMANOVA + PERMDISP
#
# Community matrix:
#   ASV relative abundance
#
# Distance:
#   Bray-Curtis
#
# PERMANOVA:
#   Site + Ecological condition
#   Treatment tested AFTER Site
#
# Permutations:
#   9,999
#   restricted within Site
#
# Planned contrasts:
#   Natural vs Degraded
#   Natural vs BESE
#   Natural vs BARE
#   BESE vs BARE
#
# Multiple testing:
#   Holm correction across four planned contrasts
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(permute)
  library(dplyr)
  library(tibble)
})

set.seed(123)


# ============================================================
# 1. LOAD VERIFIED CHAPTER 1 PHYLOSEQ OBJECTS
# ============================================================

ps_bacteria <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)

stopifnot(
  nsamples(ps_bacteria) == 103,
  ntaxa(ps_bacteria) == 32530,
  nsamples(ps_fungi) == 87,
  ntaxa(ps_fungi) == 11176
)


# ============================================================
# 2. H1 — BULK-SOIL SUBSETS
# ============================================================

ps_bac_A <- subset_samples(
  ps_bacteria,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_bac_A <- prune_taxa(
  taxa_sums(ps_bac_A) > 0,
  ps_bac_A
)


ps_fun_A <- subset_samples(
  ps_fungi,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_fun_A <- prune_taxa(
  taxa_sums(ps_fun_A) > 0,
  ps_fun_A
)


stopifnot(
  nsamples(ps_bac_A) == 58,
  ntaxa(ps_bac_A) == 28358,
  nsamples(ps_fun_A) == 46,
  ntaxa(ps_fun_A) == 6495
)


# ============================================================
# 3. PREPARE BETA-DIVERSITY DATA
# ============================================================

prepare_beta <- function(ps_object) {
  
  # ----------------------------------------------------------
  # Relative abundance
  # ----------------------------------------------------------
  
  ps_rel <- transform_sample_counts(
    ps_object,
    function(x) x / sum(x)
  )
  
  # ----------------------------------------------------------
  # Community matrix
  # samples = rows
  # ASVs    = columns
  # ----------------------------------------------------------
  
  comm <- as(
    otu_table(ps_rel),
    "matrix"
  )
  
  if (taxa_are_rows(ps_rel)) {
    comm <- t(comm)
  }
  
  comm <- comm[
    sample_names(ps_rel),
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # Metadata
  # ----------------------------------------------------------
  
  meta <- as(
    sample_data(ps_rel),
    "data.frame"
  )
  
  meta <- meta[
    rownames(comm),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(comm),
      rownames(meta)
    )
  )
  
  meta$Treatment <- factor(
    meta$Treatment,
    levels = c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
  )
  
  meta$Site <- factor(
    meta$Site,
    levels = c(
      "ACP",
      "BASS",
      "TRUG"
    )
  )
  
  # ----------------------------------------------------------
  # Bray-Curtis
  # ----------------------------------------------------------
  
  bray <- vegan::vegdist(
    comm,
    method = "bray"
  )
  
  list(
    ps_rel = ps_rel,
    comm = comm,
    meta = meta,
    bray = bray
  )
}


bac_beta <- prepare_beta(ps_bac_A)
fun_beta <- prepare_beta(ps_fun_A)


# ============================================================
# 4. DATA CHECKS
# ============================================================

cat("\n========================================\n")
cat("H1 — BETA-DIVERSITY INPUT\n")
cat("========================================\n")

cat(
  "\nBACTERIA:",
  nrow(bac_beta$comm),
  "samples\n"
)

print(
  with(
    bac_beta$meta,
    table(
      Site,
      Treatment
    )
  )
)


cat(
  "\nFUNGI:",
  nrow(fun_beta$comm),
  "samples\n"
)

print(
  with(
    fun_beta$meta,
    table(
      Site,
      Treatment
    )
  )
)


# ============================================================
# 5. OMNIBUS PERMANOVA + PERMDISP FUNCTION
# ============================================================

run_omnibus_beta <- function(beta_object, kingdom) {
  
  meta <- beta_object$meta
  bray <- beta_object$bray
  
  # ----------------------------------------------------------
  # Restricted permutations within Site
  # ----------------------------------------------------------
  
  perm_control <- permute::how(
    nperm = 9999,
    blocks = meta$Site
  )
  
  # ----------------------------------------------------------
  # PERMANOVA
  #
  # Sequential terms:
  # Site enters first.
  # Treatment therefore represents ecological-condition
  # variation remaining AFTER accounting for Site.
  # ----------------------------------------------------------
  
  set.seed(123)
  
  permanova_terms <- vegan::adonis2(
    bray ~ Site + Treatment,
    data = meta,
    permutations = perm_control,
    by = "terms"
  )
  
  permanova_df <- as.data.frame(
    permanova_terms
  ) %>%
    tibble::rownames_to_column(
      "Term"
    )
  
  treatment_row <- permanova_df %>%
    filter(
      Term == "Treatment"
    )
  
  site_row <- permanova_df %>%
    filter(
      Term == "Site"
    )
  
  # ----------------------------------------------------------
  # PERMDISP
  #
  # Test whether ecological conditions differ in
  # multivariate dispersion.
  #
  # Permutations also restricted within Site.
  # ----------------------------------------------------------
  
  dispersion <- vegan::betadisper(
    bray,
    group = meta$Treatment,
    type = "median",
    bias.adjust = TRUE
  )
  
  set.seed(123)
  
  dispersion_test <- vegan::permutest(
    dispersion,
    permutations = perm_control
  )
  
  dispersion_df <- as.data.frame(
    dispersion_test$tab
  ) %>%
    tibble::rownames_to_column(
      "Term"
    )
  
  dispersion_row <- dispersion_df %>%
    filter(
      Term == "Groups"
    )
  
  # ----------------------------------------------------------
  # Main result table
  # ----------------------------------------------------------
  
  result <- tibble(
    Kingdom = kingdom,
    
    n = nrow(meta),
    
    PERMANOVA_F =
      treatment_row$F,
    
    PERMANOVA_R2 =
      treatment_row$R2,
    
    PERMANOVA_R2_percent =
      100 * treatment_row$R2,
    
    PERMANOVA_p =
      treatment_row$`Pr(>F)`,
    
    PERMDISP_F =
      dispersion_row$F,
    
    PERMDISP_p =
      dispersion_row$`Pr(>F)`
  )
  
  list(
    result = result,
    permanova = permanova_terms,
    permanova_df = permanova_df,
    dispersion = dispersion,
    dispersion_test = dispersion_test,
    dispersion_df = dispersion_df,
    site_row = site_row
  )
}


# ============================================================
# 6. RUN OMNIBUS TESTS
# ============================================================

omnibus_bac <- run_omnibus_beta(
  bac_beta,
  "Bacteria"
)

omnibus_fun <- run_omnibus_beta(
  fun_beta,
  "Fungi"
)


cat("\n========================================\n")
cat("BACTERIA — PERMANOVA\n")
cat("========================================\n")

print(
  omnibus_bac$permanova
)


cat("\n========================================\n")
cat("BACTERIA — PERMDISP\n")
cat("========================================\n")

print(
  omnibus_bac$dispersion_test
)


cat("\n========================================\n")
cat("FUNGI — PERMANOVA\n")
cat("========================================\n")

print(
  omnibus_fun$permanova
)


cat("\n========================================\n")
cat("FUNGI — PERMDISP\n")
cat("========================================\n")

print(
  omnibus_fun$dispersion_test
)


# ============================================================
# 7. COMBINE OMNIBUS RESULTS
# ============================================================

omnibus_results <- bind_rows(
  omnibus_bac$result,
  omnibus_fun$result
)

cat("\n========================================\n")
cat("H1 — OMNIBUS SUMMARY\n")
cat("========================================\n")

print(
  omnibus_results,
  n = Inf,
  width = Inf
)


# ============================================================
# 8. PLANNED PAIRWISE PERMANOVA + PERMDISP
# ============================================================

planned_contrasts <- list(
  "Natural vs Degraded" =
    c("NATURAL", "DEGRADED"),
  
  "Natural vs BESE" =
    c("NATURAL", "BESE"),
  
  "Natural vs BARE" =
    c("NATURAL", "BARE"),
  
  "BESE vs BARE" =
    c("BESE", "BARE")
)


run_pairwise_beta <- function(
    beta_object,
    kingdom,
    contrast_name,
    groups
) {
  
  meta_full <- beta_object$meta
  comm_full <- beta_object$comm
  
  # ----------------------------------------------------------
  # Subset samples
  # ----------------------------------------------------------
  
  keep <- meta_full$Treatment %in% groups
  
  meta <- droplevels(
    meta_full[
      keep,
      ,
      drop = FALSE
    ]
  )
  
  comm <- comm_full[
    rownames(meta),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(comm),
      rownames(meta)
    )
  )
  
  # Recalculate Bray on exact subset
  bray <- vegan::vegdist(
    comm,
    method = "bray"
  )
  
  # ----------------------------------------------------------
  # Restricted permutations
  # ----------------------------------------------------------
  
  perm_control <- permute::how(
    nperm = 9999,
    blocks = meta$Site
  )
  
  # ----------------------------------------------------------
  # PERMANOVA
  # Treatment after Site
  # ----------------------------------------------------------
  
  set.seed(123)
  
  perm_result <- vegan::adonis2(
    bray ~ Site + Treatment,
    data = meta,
    permutations = perm_control,
    by = "terms"
  )
  
  perm_df <- as.data.frame(
    perm_result
  ) %>%
    tibble::rownames_to_column(
      "Term"
    )
  
  treatment_row <- perm_df %>%
    filter(
      Term == "Treatment"
    )
  
  # ----------------------------------------------------------
  # PERMDISP
  # ----------------------------------------------------------
  
  dispersion <- vegan::betadisper(
    bray,
    group = meta$Treatment,
    type = "median",
    bias.adjust = TRUE
  )
  
  set.seed(123)
  
  disp_result <- vegan::permutest(
    dispersion,
    permutations = perm_control
  )
  
  disp_df <- as.data.frame(
    disp_result$tab
  ) %>%
    tibble::rownames_to_column(
      "Term"
    )
  
  disp_row <- disp_df %>%
    filter(
      Term == "Groups"
    )
  
  tibble(
    Kingdom = kingdom,
    
    Contrast = contrast_name,
    
    Group1 = groups[1],
    
    Group2 = groups[2],
    
    n = nrow(meta),
    
    PERMANOVA_F =
      treatment_row$F,
    
    PERMANOVA_R2 =
      treatment_row$R2,
    
    PERMANOVA_R2_percent =
      100 * treatment_row$R2,
    
    PERMANOVA_p =
      treatment_row$`Pr(>F)`,
    
    PERMDISP_F =
      disp_row$F,
    
    PERMDISP_p =
      disp_row$`Pr(>F)`
  )
}


# ============================================================
# 9. RUN ALL PLANNED CONTRASTS
# ============================================================

pairwise_bac <- bind_rows(
  lapply(
    names(planned_contrasts),
    function(x) {
      run_pairwise_beta(
        bac_beta,
        "Bacteria",
        x,
        planned_contrasts[[x]]
      )
    }
  )
)


pairwise_fun <- bind_rows(
  lapply(
    names(planned_contrasts),
    function(x) {
      run_pairwise_beta(
        fun_beta,
        "Fungi",
        x,
        planned_contrasts[[x]]
      )
    }
  )
)


# ============================================================
# 10. HOLM CORRECTION
#     separately within each kingdom
# ============================================================

pairwise_bac <- pairwise_bac %>%
  mutate(
    PERMANOVA_p_Holm =
      p.adjust(
        PERMANOVA_p,
        method = "holm"
      ),
    
    PERMDISP_p_Holm =
      p.adjust(
        PERMDISP_p,
        method = "holm"
      )
  )


pairwise_fun <- pairwise_fun %>%
  mutate(
    PERMANOVA_p_Holm =
      p.adjust(
        PERMANOVA_p,
        method = "holm"
      ),
    
    PERMDISP_p_Holm =
      p.adjust(
        PERMDISP_p,
        method = "holm"
      )
  )


pairwise_results <- bind_rows(
  pairwise_bac,
  pairwise_fun
)


cat("\n========================================\n")
cat("H1 — PLANNED CONTRASTS\n")
cat("========================================\n")

print(
  pairwise_results,
  n = Inf,
  width = Inf
)


# ============================================================
# 11. FORMATTING FUNCTIONS FOR PUBLICATION TABLE
# ============================================================

format_p <- function(x) {
  
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(
      x < 0.001,
      "<0.001",
      sprintf("%.3f", x)
    )
  )
}


sig_symbol <- function(x) {
  
  case_when(
    is.na(x) ~ "",
    x < 0.001 ~ "***",
    x < 0.01  ~ "**",
    x < 0.05  ~ "*",
    TRUE      ~ ""
  )
}


# ============================================================
# 12. MAIN TABLE — LIKE YOUR EXAMPLE
# ============================================================

main_table_long <- omnibus_results %>%
  mutate(
    Factor =
      "Ecological condition (after accounting for Site)",
    
    F =
      round(
        PERMANOVA_F,
        2
      ),
    
    R2_percent =
      round(
        PERMANOVA_R2_percent,
        1
      ),
    
    p =
      format_p(
        PERMANOVA_p
      ),
    
    Significance =
      sig_symbol(
        PERMANOVA_p
      ),
    
    PERMDISP_F =
      round(
        PERMDISP_F,
        2
      ),
    
    PERMDISP_p_formatted =
      format_p(
        PERMDISP_p
      ),
    
    PERMDISP_significance =
      sig_symbol(
        PERMDISP_p
      )
  ) %>%
  select(
    Kingdom,
    Factor,
    F,
    R2_percent,
    p,
    Significance,
    PERMDISP_F,
    PERMDISP_p_formatted,
    PERMDISP_significance
  )


# ------------------------------------------------------------
# Wide version for manuscript table
# ------------------------------------------------------------

bac_main <- main_table_long %>%
  filter(
    Kingdom == "Bacteria"
  )

fun_main <- main_table_long %>%
  filter(
    Kingdom == "Fungi"
  )


publication_table <- tibble(
  Factor = c(
    "Ecological condition",
    "PERMDISP"
  ),
  
  Bacteria_F = c(
    bac_main$F,
    bac_main$PERMDISP_F
  ),
  
  Bacteria_R2_percent = c(
    bac_main$R2_percent,
    NA
  ),
  
  Bacteria_p = c(
    paste0(
      bac_main$p,
      bac_main$Significance
    ),
    paste0(
      bac_main$PERMDISP_p_formatted,
      bac_main$PERMDISP_significance
    )
  ),
  
  Fungi_F = c(
    fun_main$F,
    fun_main$PERMDISP_F
  ),
  
  Fungi_R2_percent = c(
    fun_main$R2_percent,
    NA
  ),
  
  Fungi_p = c(
    paste0(
      fun_main$p,
      fun_main$Significance
    ),
    paste0(
      fun_main$PERMDISP_p_formatted,
      fun_main$PERMDISP_significance
    )
  )
)


cat("\n========================================\n")
cat("PUBLICATION-STYLE H1 TABLE\n")
cat("========================================\n")

print(
  publication_table,
  n = Inf,
  width = Inf
)


# ============================================================
# 13. PAIRWISE PUBLICATION TABLE
# ============================================================

pairwise_publication <- pairwise_results %>%
  mutate(
    F =
      round(
        PERMANOVA_F,
        2
      ),
    
    R2_percent =
      round(
        PERMANOVA_R2_percent,
        1
      ),
    
    p_Holm_formatted =
      format_p(
        PERMANOVA_p_Holm
      ),
    
    significance =
      sig_symbol(
        PERMANOVA_p_Holm
      ),
    
    PERMDISP_F =
      round(
        PERMDISP_F,
        2
      ),
    
    PERMDISP_p_Holm_formatted =
      format_p(
        PERMDISP_p_Holm
      ),
    
    PERMDISP_significance =
      sig_symbol(
        PERMDISP_p_Holm
      )
  ) %>%
  select(
    Kingdom,
    Contrast,
    n,
    F,
    R2_percent,
    p_Holm_formatted,
    significance,
    PERMDISP_F,
    PERMDISP_p_Holm_formatted,
    PERMDISP_significance
  )


cat("\n========================================\n")
cat("PAIRWISE PUBLICATION TABLE\n")
cat("========================================\n")

print(
  pairwise_publication,
  n = Inf,
  width = Inf
)


# ============================================================
# 14. OUTPUT DIRECTORIES
# ============================================================

dir.create(
  "outputs/H1_beta/PERMANOVA",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "tables/H1_beta",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 15. SAVE RAW RESULTS
# ============================================================

write.csv(
  omnibus_results,
  "outputs/H1_beta/PERMANOVA/H1_omnibus_PERMANOVA_PERMDISP_raw.csv",
  row.names = FALSE
)

write.csv(
  pairwise_results,
  "outputs/H1_beta/PERMANOVA/H1_pairwise_PERMANOVA_PERMDISP_raw.csv",
  row.names = FALSE
)


# ============================================================
# 16. SAVE PUBLICATION TABLES
# ============================================================

write.csv(
  publication_table,
  "tables/H1_beta/H1_PERMANOVA_main_table.csv",
  row.names = FALSE
)

write.csv(
  pairwise_publication,
  "tables/H1_beta/H1_PERMANOVA_pairwise_table.csv",
  row.names = FALSE
)


# ============================================================
# 17. SAVE MODEL OBJECTS
# ============================================================

saveRDS(
  omnibus_bac,
  "outputs/H1_beta/PERMANOVA/H1_bacteria_PERMANOVA_PERMDISP.rds"
)

saveRDS(
  omnibus_fun,
  "outputs/H1_beta/PERMANOVA/H1_fungi_PERMANOVA_PERMDISP.rds"
)


# ============================================================
# 18. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file =
    "outputs/H1_beta/PERMANOVA/sessionInfo_H1_PERMANOVA.txt"
)


cat("\n========================================\n")
cat("H1 PERMANOVA + PERMDISP COMPLETE\n")
cat("========================================\n")


# ============================================================
# H1 — BULK-SOIL REASSEMBLY
# DISTANCE TO NATURAL
#
# Ecological question:
# Are restoration communities closer to the Natural reference
# than Degraded communities?
#
# Distance:
# Bray-Curtis on relative-abundance ASV data
#
# Reference:
# Natural samples from the SAME SITE
#
# Natural samples:
# Leave-one-out mean distance to the other Natural samples
#
# Lower distance = more Natural-like composition
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tibble)
  library(car)
  library(emmeans)
})

set.seed(123)


# ============================================================
# 1. LOAD VERIFIED CHAPTER 1 PHYLOSEQ OBJECTS
# ============================================================

ps_bacteria <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)

stopifnot(
  nsamples(ps_bacteria) == 103,
  ntaxa(ps_bacteria) == 32530,
  nsamples(ps_fungi) == 87,
  ntaxa(ps_fungi) == 11176
)


# ============================================================
# 2. H1 — BULK-SOIL SUBSETS
# ============================================================

ps_bac_A <- subset_samples(
  ps_bacteria,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_bac_A <- prune_taxa(
  taxa_sums(ps_bac_A) > 0,
  ps_bac_A
)


ps_fun_A <- subset_samples(
  ps_fungi,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_fun_A <- prune_taxa(
  taxa_sums(ps_fun_A) > 0,
  ps_fun_A
)


stopifnot(
  nsamples(ps_bac_A) == 58,
  ntaxa(ps_bac_A) == 28358,
  nsamples(ps_fun_A) == 46,
  ntaxa(ps_fun_A) == 6495
)


# ============================================================
# 3. FUNCTION — DISTANCE TO SAME-SITE NATURAL REFERENCE
# ============================================================

calculate_distance_to_natural <- function(ps_object, kingdom) {
  
  # ----------------------------------------------------------
  # Relative abundance
  # ----------------------------------------------------------
  
  ps_rel <- phyloseq::transform_sample_counts(
    ps_object,
    function(x) x / sum(x)
  )
  
  
  
  # ----------------------------------------------------------
  # Community matrix
  # ----------------------------------------------------------
  
  comm <- as(
    otu_table(ps_rel),
    "matrix"
  )
  
  if (taxa_are_rows(ps_rel)) {
    comm <- t(comm)
  }
  
  comm <- comm[
    sample_names(ps_rel),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Metadata
  # ----------------------------------------------------------
  
  meta <- as(
    sample_data(ps_rel),
    "data.frame"
  )
  
  meta <- meta[
    rownames(comm),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(comm),
      rownames(meta)
    )
  )
  
  
  meta$Treatment <- factor(
    meta$Treatment,
    levels = c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
  )
  
  meta$Site <- factor(
    meta$Site,
    levels = c(
      "ACP",
      "BASS",
      "TRUG"
    )
  )
  
  
  # ----------------------------------------------------------
  # Bray-Curtis distance matrix
  # ----------------------------------------------------------
  
  bray <- vegan::vegdist(
    comm,
    method = "bray"
  )
  
  bray_matrix <- as.matrix(
    bray
  )
  
  
  # ----------------------------------------------------------
  # Calculate distance to Natural
  #
  # Non-Natural sample:
  # mean distance to ALL Natural samples from same Site
  #
  # Natural sample:
  # mean distance to OTHER Natural samples from same Site
  # (leave-one-out)
  # ----------------------------------------------------------
  
  sample_ids <- rownames(meta)
  
  distance_values <- numeric(
    length(sample_ids)
  )
  
  n_reference_samples <- integer(
    length(sample_ids)
  )
  
  
  for (i in seq_along(sample_ids)) {
    
    sample_id <- sample_ids[i]
    
    this_site <- as.character(
      meta[sample_id, "Site"]
    )
    
    this_treatment <- as.character(
      meta[sample_id, "Treatment"]
    )
    
    
    natural_ids <- rownames(meta)[
      meta$Site == this_site &
        meta$Treatment == "NATURAL"
    ]
    
    
    # Leave-one-out for Natural samples
    if (this_treatment == "NATURAL") {
      
      natural_ids <- setdiff(
        natural_ids,
        sample_id
      )
      
    }
    
    
    if (length(natural_ids) == 0) {
      
      stop(
        paste(
          "No valid Natural reference samples available for",
          sample_id
        )
      )
      
    }
    
    
    distance_values[i] <- mean(
      bray_matrix[
        sample_id,
        natural_ids,
        drop = TRUE
      ]
    )
    
    
    n_reference_samples[i] <- length(
      natural_ids
    )
  }
  
  
  # ----------------------------------------------------------
  # Assemble output
  # ----------------------------------------------------------
  
  results <- meta %>%
    tibble::rownames_to_column(
      "SampleID"
    ) %>%
    mutate(
      Distance_to_Natural =
        distance_values,
      
      n_Natural_references =
        n_reference_samples,
      
      Kingdom =
        kingdom
    )
  
  
  # ----------------------------------------------------------
  # Checks
  # ----------------------------------------------------------
  
  stopifnot(
    nrow(results) == nsamples(ps_object),
    sum(is.na(results$Distance_to_Natural)) == 0
  )
  
  
  list(
    data = results,
    bray = bray,
    bray_matrix = bray_matrix,
    relative_phyloseq = ps_rel
  )
}


# ============================================================
# 4. CALCULATE DISTANCE TO NATURAL
# ============================================================

dist_bac <- calculate_distance_to_natural(
  ps_bac_A,
  "Bacteria"
)

dist_fun <- calculate_distance_to_natural(
  ps_fun_A,
  "Fungi"
)

# ============================================================
# 5. VERIFY REFERENCE DESIGN
# ============================================================

cat("\n========================================\n")
cat("BACTERIA — NATURAL REFERENCE COUNTS\n")
cat("========================================\n")

print(
  dist_bac$data %>%
    count(
      Site,
      Treatment,
      n_Natural_references
    )
)


cat("\n========================================\n")
cat("FUNGI — NATURAL REFERENCE COUNTS\n")
cat("========================================\n")

print(
  dist_fun$data %>%
    count(
      Site,
      Treatment,
      n_Natural_references
    )
)


# ============================================================
# 6. DESCRIPTIVE SUMMARY
# ============================================================

distance_summary <- bind_rows(
  dist_bac$data,
  dist_fun$data
) %>%
  group_by(
    Kingdom,
    Treatment
  ) %>%
  summarise(
    n = n(),
    mean_distance =
      mean(
        Distance_to_Natural
      ),
    SD =
      sd(
        Distance_to_Natural
      ),
    median_distance =
      median(
        Distance_to_Natural
      ),
    Q1 =
      quantile(
        Distance_to_Natural,
        0.25
      ),
    Q3 =
      quantile(
        Distance_to_Natural,
        0.75
      ),
    .groups = "drop"
  )


cat("\n========================================\n")
cat("DISTANCE TO NATURAL — DESCRIPTIVE SUMMARY\n")
cat("========================================\n")

print(
  distance_summary,
  n = Inf,
  width = Inf
)


# ============================================================
# 7. MODEL FUNCTION
#
# Retain Site × Treatment initially because we need to test
# whether the ecological-condition pattern differs by Site.
# ============================================================

run_distance_model <- function(distance_object, kingdom) {
  
  dat <- distance_object$data
  
  
  mod <- lm(
    Distance_to_Natural ~ Site * Treatment,
    data = dat
  )
  
  
  # ----------------------------------------------------------
  # Type II ANOVA
  # ----------------------------------------------------------
  
  anova_result <- car::Anova(
    mod,
    type = 2
  )
  
  
  cat("\n========================================\n")
  cat(
    kingdom,
    "— TYPE II ANOVA\n"
  )
  cat("========================================\n")
  
  print(
    anova_result
  )
  
  
  # ----------------------------------------------------------
  # Estimated marginal means
  #
  # Equal weighting across sites
  # ----------------------------------------------------------
  
  emm <- emmeans(
    mod,
    ~ Treatment,
    weights = "equal"
  )
  
  
  cat("\n========================================\n")
  cat(
    kingdom,
    "— ESTIMATED MARGINAL MEANS\n"
  )
  cat("========================================\n")
  
  print(
    summary(
      emm,
      infer = TRUE
    )
  )
  
  
  # ==========================================================
  # 8A. PRIMARY RECOVERY CONTRASTS
  #
  # Negative estimate:
  # first group is CLOSER to Natural than second group
  # ==========================================================
  
  recovery_contrasts <- contrast(
    emm,
    
    method = list(
      
      "BESE - Degraded" =
        c(
          0,
          1,
          0,
          -1
        ),
      
      "BARE - Degraded" =
        c(
          0,
          0,
          1,
          -1
        ),
      
      "BESE - BARE" =
        c(
          0,
          1,
          -1,
          0
        )
      
    ),
    
    adjust = "holm"
  )
  
  
  recovery_results <- as.data.frame(
    summary(
      recovery_contrasts,
      infer = TRUE
    )
  )
  
  
  cat("\n========================================\n")
  cat(
    kingdom,
    "— PRIMARY RECOVERY CONTRASTS\n"
  )
  cat("Holm-adjusted p-values\n")
  cat("Negative estimate = first group closer to Natural\n")
  cat("========================================\n")
  
  print(
    recovery_results
  )
  
  
  # ==========================================================
  # 8B. DISTANCE FROM NATURAL BASELINE
  #
  # These tell us whether each ecological condition remains
  # farther from Natural than Natural samples are from one
  # another.
  #
  # NOTE:
  # non-significance does NOT establish equivalence/recovery.
  # ==========================================================
  
  natural_contrasts <- contrast(
    emm,
    
    method = list(
      
      "Natural - BESE" =
        c(
          1,
          -1,
          0,
          0
        ),
      
      "Natural - BARE" =
        c(
          1,
          0,
          -1,
          0
        ),
      
      "Natural - Degraded" =
        c(
          1,
          0,
          0,
          -1
        )
      
    ),
    
    adjust = "holm"
  )
  
  
  natural_results <- as.data.frame(
    summary(
      natural_contrasts,
      infer = TRUE
    )
  )
  
  
  cat("\n========================================\n")
  cat(
    kingdom,
    "— NATURAL BASELINE CONTRASTS\n"
  )
  cat("Holm-adjusted p-values\n")
  cat("========================================\n")
  
  print(
    natural_results
  )
  
  
  # ----------------------------------------------------------
  # Site-specific EMMs — SAVE for audit/supplement
  # ----------------------------------------------------------
  
  site_emm <- emmeans(
    mod,
    ~ Treatment | Site
  )
  
  site_pairs <- pairs(
    site_emm,
    adjust = "tukey"
  )
  
  
  list(
    model = mod,
    anova = anova_result,
    emm = emm,
    recovery_contrasts = recovery_contrasts,
    recovery_results = recovery_results,
    natural_contrasts = natural_contrasts,
    natural_results = natural_results,
    site_emm = site_emm,
    site_pairs = site_pairs
  )
}


# ============================================================
# 9. RUN MODELS
# ============================================================

model_bac <- run_distance_model(
  dist_bac,
  "BACTERIA"
)

model_fun <- run_distance_model(
  dist_fun,
  "FUNGI"
)


# ============================================================
# 10. OUTPUT DIRECTORIES
# ============================================================

dir.create(
  "outputs/H1_beta/distance_to_Natural",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "tables/H1_beta",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 11. SAVE SAMPLE-LEVEL DATA
# ============================================================

write.csv(
  dist_bac$data,
  "outputs/H1_beta/distance_to_Natural/H1_bacteria_distance_to_Natural.csv",
  row.names = FALSE
)

write.csv(
  dist_fun$data,
  "outputs/H1_beta/distance_to_Natural/H1_fungi_distance_to_Natural.csv",
  row.names = FALSE
)


# ============================================================
# 12. SAVE DESCRIPTIVE SUMMARY
# ============================================================

write.csv(
  distance_summary,
  "tables/H1_beta/H1_distance_to_Natural_summary.csv",
  row.names = FALSE
)


# ============================================================
# 13. SAVE PRIMARY RECOVERY CONTRASTS
# ============================================================

bac_recovery <- model_bac$recovery_results %>%
  mutate(
    Kingdom = "Bacteria"
  )

fun_recovery <- model_fun$recovery_results %>%
  mutate(
    Kingdom = "Fungi"
  )

recovery_table <- bind_rows(
  bac_recovery,
  fun_recovery
)


write.csv(
  recovery_table,
  "tables/H1_beta/H1_distance_to_Natural_recovery_contrasts.csv",
  row.names = FALSE
)


# ============================================================
# 14. SAVE NATURAL BASELINE CONTRASTS
# ============================================================

bac_natural <- model_bac$natural_results %>%
  mutate(
    Kingdom = "Bacteria"
  )

fun_natural <- model_fun$natural_results %>%
  mutate(
    Kingdom = "Fungi"
  )

natural_table <- bind_rows(
  bac_natural,
  fun_natural
)


write.csv(
  natural_table,
  "tables/H1_beta/H1_distance_to_Natural_baseline_contrasts.csv",
  row.names = FALSE
)


# ============================================================
# 15. SAVE SITE-SPECIFIC RESULTS FOR SUPPLEMENT / AUDIT
# ============================================================

write.csv(
  as.data.frame(
    model_bac$site_emm
  ),
  "outputs/H1_beta/distance_to_Natural/H1_bacteria_site_EMMs.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(
    model_fun$site_emm
  ),
  "outputs/H1_beta/distance_to_Natural/H1_fungi_site_EMMs.csv",
  row.names = FALSE
)


write.csv(
  as.data.frame(
    model_bac$site_pairs
  ),
  "outputs/H1_beta/distance_to_Natural/H1_bacteria_site_pairwise.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(
    model_fun$site_pairs
  ),
  "outputs/H1_beta/distance_to_Natural/H1_fungi_site_pairwise.csv",
  row.names = FALSE
)


# ============================================================
# 16. SAVE MODELS
# ============================================================

saveRDS(
  model_bac,
  "outputs/H1_beta/distance_to_Natural/H1_bacteria_distance_model.rds"
)

saveRDS(
  model_fun,
  "outputs/H1_beta/distance_to_Natural/H1_fungi_distance_model.rds"
)


# ============================================================
# 17. SAVE BRAY-CURTIS MATRICES
# ============================================================

saveRDS(
  dist_bac$bray,
  "outputs/H1_beta/distance_to_Natural/H1_bacteria_BrayCurtis.rds"
)

saveRDS(
  dist_fun$bray,
  "outputs/H1_beta/distance_to_Natural/H1_fungi_BrayCurtis.rds"
)


# ============================================================
# 18. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file =
    "outputs/H1_beta/distance_to_Natural/sessionInfo_H1_distance_to_Natural.txt"
)


cat("\n========================================\n")
cat("H1 — DISTANCE TO NATURAL COMPLETE\n")
cat("========================================\n")


# ============================================================
# SAVE H1 — DISTANCE TO NATURAL TABLES
# ============================================================

dir.create(
  "outputs/H1_beta/distance_to_Natural",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "tables/H1_beta",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. SAMPLE-LEVEL DATA
# ============================================================

write.csv(
  dist_bac$data,
  "outputs/H1_beta/distance_to_Natural/H1_bacteria_distance_to_Natural.csv",
  row.names = FALSE
)

write.csv(
  dist_fun$data,
  "outputs/H1_beta/distance_to_Natural/H1_fungi_distance_to_Natural.csv",
  row.names = FALSE
)


# ============================================================
# 2. DESCRIPTIVE SUMMARY TABLE
# ============================================================

distance_summary <- dplyr::bind_rows(
  dist_bac$data,
  dist_fun$data
) %>%
  dplyr::group_by(
    Kingdom,
    Treatment
  ) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_distance = mean(Distance_to_Natural),
    SD = sd(Distance_to_Natural),
    median_distance = median(Distance_to_Natural),
    Q1 = quantile(Distance_to_Natural, 0.25),
    Q3 = quantile(Distance_to_Natural, 0.75),
    .groups = "drop"
  )


write.csv(
  distance_summary,
  "tables/H1_beta/H1_distance_to_Natural_descriptive_summary.csv",
  row.names = FALSE
)


# ============================================================
# 3. CHECK
# ============================================================

cat("\n========================================\n")
cat("DISTANCE-TO-NATURAL TABLES SAVED\n")
cat("========================================\n")

print(distance_summary)

cat(
  "\nSaved to:\n",
  "outputs/H1_beta/distance_to_Natural/\n",
  "tables/H1_beta/\n"
)

#ACOMBC 


# ============================================================
# H1 — BACTERIA — GENUS-LEVEL ANCOM-BC2
# FINAL LOCKED ANALYSIS
#
# Research question:
# Which bacterial genera differ in bulk soil among
# Natural, BESE, BARE and Degraded ecological conditions,
# after accounting for Site?
#
# Workflow:
# 1. Verified historical bacterial ASV table
# 2. Bulk-soil samples only
# 3. Raw counts
# 4. Keep taxonomically resolved genera
# 5. Aggregate ASVs -> Genus
# 6. Explicit >=10% prevalence filter
# 7. ANCOM-BC2: Site + Condition
# 8. BH multiple-testing correction
# 9. Pseudo-count sensitivity analysis
# 10. Conservative structural-zero detection
#
# Software currently:
# R 4.6.1
# ANCOMBC 2.14.0
# ============================================================


# ============================================================
# 0. SETUP
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
  library(tibble)
})

set.seed(123)


# ============================================================
# 1. SOFTWARE RECORD
# ============================================================

cat("\n========================================\n")
cat("SOFTWARE\n")
cat("========================================\n")

cat(
  "R version:",
  R.version.string,
  "\n"
)

cat(
  "ANCOMBC version:",
  as.character(packageVersion("ANCOMBC")),
  "\n"
)


# ============================================================
# 2. LOAD VERIFIED HISTORICAL BACTERIAL DATASET
# ============================================================

ps_bacteria <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

stopifnot(
  nsamples(ps_bacteria) == 103,
  ntaxa(ps_bacteria) == 32530
)

cat("\n========================================\n")
cat("FULL VERIFIED BACTERIAL DATASET\n")
cat("========================================\n")

cat(
  "Samples:",
  nsamples(ps_bacteria),
  "\n"
)

cat(
  "ASVs:",
  ntaxa(ps_bacteria),
  "\n"
)

cat(
  "Reads:",
  sum(sample_sums(ps_bacteria)),
  "\n"
)


# ============================================================
# 3. H1 — BULK SOIL ONLY
# ============================================================

ps_A <- subset_samples(
  ps_bacteria,
  SoilPortion == "BS" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
)

ps_A <- prune_taxa(
  taxa_sums(ps_A) > 0,
  ps_A
)

stopifnot(
  nsamples(ps_A) == 58,
  ntaxa(ps_A) == 28358
)

cat("\n========================================\n")
cat("H1 INPUT\n")
cat("========================================\n")

cat(
  "Bulk-soil samples:",
  nsamples(ps_A),
  "\n"
)

cat(
  "ASVs detected in bulk soil:",
  ntaxa(ps_A),
  "\n"
)

cat(
  "Reads:",
  sum(sample_sums(ps_A)),
  "\n"
)


# ============================================================
# 4. EXTRACT RAW ASV COUNTS
#
# Final orientation:
# samples x ASVs
# ============================================================

seqA <- as(
  otu_table(ps_A),
  "matrix"
)

if (taxa_are_rows(ps_A)) {
  seqA <- t(seqA)
}

seqA <- seqA[
  sample_names(ps_A),
  taxa_names(ps_A),
  drop = FALSE
]

stopifnot(
  nrow(seqA) == 58,
  identical(
    rownames(seqA),
    sample_names(ps_A)
  )
)


# ============================================================
# 5. EXTRACT TAXONOMY
# ============================================================

taxA <- as(
  tax_table(ps_A),
  "matrix"
)

taxA <- taxA[
  taxa_names(ps_A),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(seqA),
    rownames(taxA)
  ),
  "Genus" %in% colnames(taxA)
)


# ============================================================
# 6. RETAIN GENUS-ASSIGNED ASVs
#
# Unresolved genus annotations are not assigned a fake genus.
# They are excluded from genus-level DA analysis.
# ============================================================

genus <- taxA[, "Genus"]

assigned <-
  !is.na(genus) &
  trimws(genus) != ""

seq_genus <- seqA[
  ,
  assigned,
  drop = FALSE
]

genus <- trimws(
  genus[assigned]
)

cat("\n========================================\n")
cat("GENUS ASSIGNMENT\n")
cat("========================================\n")

cat(
  "ASVs assigned to a resolved genus:",
  ncol(seq_genus),
  "\n"
)

cat(
  "Unique genus labels:",
  length(unique(genus)),
  "\n"
)

cat(
  "Reads assigned to resolved genera:",
  sum(seq_genus),
  "\n"
)


# ============================================================
# 7. AGGREGATE RAW COUNTS FROM ASV -> GENUS
#
# Result:
# genera x samples
# ============================================================

genus_counts_all <- rowsum(
  t(seq_genus),
  group = genus,
  reorder = FALSE
)

cat(
  "Genus count matrix:",
  nrow(genus_counts_all),
  "genera x",
  ncol(genus_counts_all),
  "samples\n"
)


# ============================================================
# 8. EXPLICIT 10% PREVALENCE FILTER
#
# 58 samples:
# ceiling(0.10 * 58) = 6 samples
#
# This is done BEFORE ANCOM-BC2 so the analytical
# taxon population is transparent and reproducible.
# ============================================================

prevalence_threshold <- ceiling(
  0.10 * ncol(genus_counts_all)
)

genus_prevalence <- rowSums(
  genus_counts_all > 0
)

keep_genus <-
  genus_prevalence >= prevalence_threshold

genus_counts <- genus_counts_all[
  keep_genus,
  ,
  drop = FALSE
]

cat("\n========================================\n")
cat("PREVALENCE FILTER\n")
cat("========================================\n")

cat(
  "Threshold:",
  prevalence_threshold,
  "of",
  ncol(genus_counts_all),
  "samples\n"
)

cat(
  "Resolved genera before prevalence filtering:",
  nrow(genus_counts_all),
  "\n"
)

cat(
  "Genera passing >=10% prevalence:",
  nrow(genus_counts),
  "\n"
)

cat(
  "Genus-assigned reads before filter:",
  sum(genus_counts_all),
  "\n"
)

cat(
  "Reads retained after prevalence filter:",
  sum(genus_counts),
  "\n"
)

cat(
  "Percentage of genus-assigned reads retained:",
  round(
    100 *
      sum(genus_counts) /
      sum(genus_counts_all),
    2
  ),
  "%\n"
)


# ============================================================
# 9. CREATE STABLE GENUS IDs
#
# Avoid potential problems caused by unusual genus names
# inside ANCOM-BC2.
# ============================================================

genus_names <- rownames(genus_counts)

genus_ids <- paste0(
  "GENUS",
  seq_len(nrow(genus_counts))
)

rownames(genus_counts) <- genus_ids

genus_map <- tibble(
  taxon = genus_ids,
  Genus = genus_names,
  Prevalence = as.integer(
    genus_prevalence[keep_genus]
  ),
  Total_reads = as.numeric(
    rowSums(genus_counts)
  )
)


# ============================================================
# 10. METADATA
#
# Natural is reference.
#
# Site is included as an adjustment variable.
# Natural is interpreted as a reference ecological condition,
# not as a randomized fourth restoration treatment.
# ============================================================

meta <- as(
  sample_data(ps_A),
  "data.frame"
)

meta <- meta[
  rownames(seqA),
  ,
  drop = FALSE
]

meta_ancom <- data.frame(
  SampleID = rownames(meta),
  
  Site = factor(
    meta$Site,
    levels = c(
      "ACP",
      "BASS",
      "TRUG"
    )
  ),
  
  Condition = factor(
    meta$Treatment,
    levels = c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
  ),
  
  stringsAsFactors = FALSE
)

rownames(meta_ancom) <-
  meta_ancom$SampleID

genus_counts <- genus_counts[
  ,
  rownames(meta_ancom),
  drop = FALSE
]

stopifnot(
  identical(
    colnames(genus_counts),
    rownames(meta_ancom)
  )
)

cat("\n========================================\n")
cat("SAMPLE DESIGN\n")
cat("========================================\n")

print(
  table(
    meta_ancom$Condition
  )
)

cat("\nSite x ecological condition:\n")

print(
  table(
    meta_ancom$Site,
    meta_ancom$Condition
  )
)


# ============================================================
# 11. OUTPUT DIRECTORY
# ============================================================

outdir <-
  "outputs/H1_ANCOMBC2/bacteria_genus_FINAL"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Save the EXACT analytical inputs
saveRDS(
  genus_counts,
  file.path(
    outdir,
    "H1_Bacteria_Genus_counts_10pct.rds"
  )
)

write.csv(
  genus_map,
  file.path(
    outdir,
    "H1_Bacteria_Genus_map.csv"
  ),
  row.names = FALSE
)

write.csv(
  meta_ancom,
  file.path(
    outdir,
    "H1_Bacteria_metadata.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. RUN ANCOM-BC2
#
# IMPORTANT:
#
# prv_cut = 0
# because the 10% prevalence filter has ALREADY been
# explicitly applied above.
#
# neg_lb = FALSE
# conservative structural-zero detection appropriate
# for this sparse genus table / modest group sizes.
#
# pseudo_sens = TRUE
# retain only findings robust to pseudo-count sensitivity.
# ============================================================

cat("\n========================================\n")
cat("RUNNING FINAL ANCOM-BC2\n")
cat("H1 — BACTERIA — GENUS\n")
cat("========================================\n")

set.seed(123)

fit_A_bac_genus <- ANCOMBC::ancombc2(
  
  data = genus_counts,
  
  taxa_are_rows = TRUE,
  
  meta_data = meta_ancom,
  
  fix_formula = "Site + Condition",
  
  rand_formula = NULL,
  
  p_adj_method = "BH",
  
  pseudo = 0,
  
  pseudo_sens = TRUE,
  
  # Explicit prevalence filtering already performed
  prv_cut = 0,
  
  lib_cut = 0,
  
  s0_perc = 0.05,
  
  group = "Condition",
  
  struc_zero = TRUE,
  
  neg_lb = FALSE,
  
  alpha = 0.05,
  
  n_cl = 1,
  
  verbose = TRUE,
  
  global = TRUE,
  
  pairwise = TRUE,
  
  dunnet = FALSE,
  
  trend = FALSE
)


# ============================================================
# 13. SAVE RAW ANCOM-BC2 OBJECT IMMEDIATELY
# ============================================================

saveRDS(
  fit_A_bac_genus,
  file.path(
    outdir,
    "H1_Bacteria_Genus_ANCOMBC2_full_object.rds"
  )
)


# ============================================================
# 14. PRIMARY QUANTITATIVE RESULTS
# ============================================================

primary <- as.data.frame(
  fit_A_bac_genus$res
)

primary_named <- primary %>%
  left_join(
    genus_map,
    by = "taxon"
  ) %>%
  relocate(
    taxon,
    Genus,
    Prevalence,
    Total_reads
  )

write.csv(
  primary_named,
  file.path(
    outdir,
    "H1_Bacteria_Genus_primary_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. GLOBAL TEST
#
# Robust =
# significant after correction AND
# passed pseudo-count sensitivity
# ============================================================

global <- as.data.frame(
  fit_A_bac_genus$res_global
)

global_named <- global %>%
  left_join(
    genus_map,
    by = "taxon"
  ) %>%
  relocate(
    taxon,
    Genus,
    Prevalence,
    Total_reads
  )

global_robust <- global_named %>%
  filter(
    diff_abn %in% TRUE,
    passed_ss %in% TRUE
  ) %>%
  arrange(
    q_val
  )

cat("\n========================================\n")
cat("ROBUST GLOBAL TEST\n")
cat("========================================\n")

cat(
  "Genera quantitatively modelled:",
  nrow(primary),
  "\n"
)

cat(
  "Robust global genera:",
  nrow(global_robust),
  "\n"
)

write.csv(
  global_named,
  file.path(
    outdir,
    "H1_Bacteria_Genus_global_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  global_robust,
  file.path(
    outdir,
    "H1_Bacteria_Genus_global_ROBUST.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. STRUCTURAL ZEROS
#
# Structural-zero taxa are reported SEPARATELY from
# quantitative differential-abundance results.
# ============================================================

zero_ind <- as.data.frame(
  fit_A_bac_genus$zero_ind
)

zero_named <- zero_ind %>%
  left_join(
    genus_map,
    by = "taxon"
  ) %>%
  relocate(
    taxon,
    Genus,
    Prevalence,
    Total_reads
  )

zero_cols <- grep(
  "^structural_zero",
  colnames(zero_ind),
  value = TRUE
)

has_structural_zero <- apply(
  zero_ind[
    ,
    zero_cols,
    drop = FALSE
  ],
  1,
  any
)

structural_zero_taxa <-
  zero_named[
    has_structural_zero,
    ,
    drop = FALSE
  ]

cat("\n========================================\n")
cat("STRUCTURAL ZEROS\n")
cat("========================================\n")

cat(
  "Prevalence-filtered genera entering ANCOM-BC2:",
  nrow(genus_counts),
  "\n"
)

cat(
  "Genera with >=1 structural zero:",
  nrow(structural_zero_taxa),
  "\n"
)

for (z in zero_cols) {
  
  cat(
    z,
    ":",
    sum(
      zero_ind[[z]],
      na.rm = TRUE
    ),
    "\n"
  )
}

write.csv(
  zero_named,
  file.path(
    outdir,
    "H1_Bacteria_Genus_structural_zero_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  structural_zero_taxa,
  file.path(
    outdir,
    "H1_Bacteria_Genus_structural_zero_taxa.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. PAIRWISE RESULTS
# ============================================================

pair <- as.data.frame(
  fit_A_bac_genus$res_pair
)

pair_named <- pair %>%
  left_join(
    genus_map,
    by = "taxon"
  ) %>%
  relocate(
    taxon,
    Genus,
    Prevalence,
    Total_reads
  )

write.csv(
  pair_named,
  file.path(
    outdir,
    "H1_Bacteria_Genus_pairwise_all.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. FUNCTION TO EXTRACT ROBUST PAIRWISE RESULTS
#
# LFC direction:
# positive = first/non-reference condition is higher
# negative = reference/second condition is higher
#
# NOTE:
# ANCOM-BC2 LFC represents bias-corrected abundance,
# not raw relative abundance.
# ============================================================

extract_robust_contrast <- function(
    pair_table,
    term,
    label_positive,
    label_negative,
    genus_map
) {
  
  lfc_col  <- paste0("lfc_", term)
  q_col    <- paste0("q_", term)
  diff_col <- paste0("diff_", term)
  ss_col   <- paste0("passed_ss_", term)
  
  required <- c(
    "taxon",
    lfc_col,
    q_col,
    diff_col,
    ss_col
  )
  
  if (!all(required %in% colnames(pair_table))) {
    
    stop(
      paste(
        "Missing ANCOM-BC2 pairwise columns for:",
        term
      )
    )
  }
  
  out <- pair_table %>%
    filter(
      .data[[diff_col]] %in% TRUE,
      .data[[ss_col]] %in% TRUE
    ) %>%
    transmute(
      taxon = taxon,
      LFC = .data[[lfc_col]],
      q = .data[[q_col]],
      Higher = ifelse(
        .data[[lfc_col]] > 0,
        label_positive,
        label_negative
      )
    ) %>%
    left_join(
      genus_map,
      by = "taxon"
    ) %>%
    select(
      taxon,
      Genus,
      LFC,
      q,
      Higher,
      Prevalence,
      Total_reads
    ) %>%
    arrange(q)
  
  out
}


# ============================================================
# 19. EXTRACT ALL SIX PAIRWISE CONTRASTS
# ============================================================

BESE_vs_NATURAL <- extract_robust_contrast(
  pair,
  "ConditionBESE",
  "BESE",
  "NATURAL",
  genus_map
)

BARE_vs_NATURAL <- extract_robust_contrast(
  pair,
  "ConditionBARE",
  "BARE",
  "NATURAL",
  genus_map
)

DEGRADED_vs_NATURAL <- extract_robust_contrast(
  pair,
  "ConditionDEGRADED",
  "DEGRADED",
  "NATURAL",
  genus_map
)

BARE_vs_BESE <- extract_robust_contrast(
  pair,
  "ConditionBARE_ConditionBESE",
  "BARE",
  "BESE",
  genus_map
)

DEGRADED_vs_BESE <- extract_robust_contrast(
  pair,
  "ConditionDEGRADED_ConditionBESE",
  "DEGRADED",
  "BESE",
  genus_map
)

DEGRADED_vs_BARE <- extract_robust_contrast(
  pair,
  "ConditionDEGRADED_ConditionBARE",
  "DEGRADED",
  "BARE",
  genus_map
)


# ============================================================
# 20. PAIRWISE SUMMARY
# ============================================================

pairwise_summary <- tibble(
  Contrast = c(
    "BESE - NATURAL",
    "BARE - NATURAL",
    "DEGRADED - NATURAL",
    "BARE - BESE",
    "DEGRADED - BESE",
    "DEGRADED - BARE"
  ),
  
  Robust_genera = c(
    nrow(BESE_vs_NATURAL),
    nrow(BARE_vs_NATURAL),
    nrow(DEGRADED_vs_NATURAL),
    nrow(BARE_vs_BESE),
    nrow(DEGRADED_vs_BESE),
    nrow(DEGRADED_vs_BARE)
  )
)

cat("\n========================================\n")
cat("ROBUST PAIRWISE GENUS RESULTS\n")
cat("========================================\n")

print(
  pairwise_summary,
  n = Inf
)


# ============================================================
# 21. PRINT ROBUST GENERA
# ============================================================

cat("\n========================================\n")
cat("BESE - NATURAL\n")
cat("========================================\n")
print(BESE_vs_NATURAL, n = Inf)

cat("\n========================================\n")
cat("BARE - NATURAL\n")
cat("========================================\n")
print(BARE_vs_NATURAL, n = Inf)

cat("\n========================================\n")
cat("DEGRADED - NATURAL\n")
cat("========================================\n")
print(DEGRADED_vs_NATURAL, n = Inf)

cat("\n========================================\n")
cat("BARE - BESE\n")
cat("========================================\n")
print(BARE_vs_BESE, n = Inf)

cat("\n========================================\n")
cat("DEGRADED - BESE\n")
cat("========================================\n")
print(DEGRADED_vs_BESE, n = Inf)

cat("\n========================================\n")
cat("DEGRADED - BARE\n")
cat("========================================\n")
print(DEGRADED_vs_BARE, n = Inf)


# ============================================================
# 22. SAVE ROBUST PAIRWISE TABLES
# ============================================================

write.csv(
  BESE_vs_NATURAL,
  file.path(
    outdir,
    "H1_Bacteria_Genus_BESE_vs_NATURAL_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  BARE_vs_NATURAL,
  file.path(
    outdir,
    "H1_Bacteria_Genus_BARE_vs_NATURAL_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  DEGRADED_vs_NATURAL,
  file.path(
    outdir,
    "H1_Bacteria_Genus_DEGRADED_vs_NATURAL_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  BARE_vs_BESE,
  file.path(
    outdir,
    "H1_Bacteria_Genus_BARE_vs_BESE_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  DEGRADED_vs_BESE,
  file.path(
    outdir,
    "H1_Bacteria_Genus_DEGRADED_vs_BESE_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  DEGRADED_vs_BARE,
  file.path(
    outdir,
    "H1_Bacteria_Genus_DEGRADED_vs_BARE_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_summary,
  file.path(
    outdir,
    "H1_Bacteria_Genus_pairwise_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. UNIQUE ROBUST PAIRWISE GENERA
# ============================================================

all_pairwise_robust <- bind_rows(
  mutate(
    BESE_vs_NATURAL,
    Contrast = "BESE - NATURAL"
  ),
  mutate(
    BARE_vs_NATURAL,
    Contrast = "BARE - NATURAL"
  ),
  mutate(
    DEGRADED_vs_NATURAL,
    Contrast = "DEGRADED - NATURAL"
  ),
  mutate(
    BARE_vs_BESE,
    Contrast = "BARE - BESE"
  ),
  mutate(
    DEGRADED_vs_BESE,
    Contrast = "DEGRADED - BESE"
  ),
  mutate(
    DEGRADED_vs_BARE,
    Contrast = "DEGRADED - BARE"
  )
)

cat("\nUnique robust pairwise genera:",
    n_distinct(all_pairwise_robust$taxon),
    "\n"
)

write.csv(
  all_pairwise_robust,
  file.path(
    outdir,
    "H1_Bacteria_Genus_all_ROBUST_pairwise_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. SAVE SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "sessionInfo_H1_Bacteria_Genus.txt"
  )
)


# ============================================================
# 25. FINAL AUDIT SUMMARY
# ============================================================

cat("\n\n========================================\n")
cat("H1 BACTERIA GENUS — FINAL SUMMARY\n")
cat("========================================\n")

cat(
  "Samples:",
  nrow(meta_ancom),
  "\n"
)

cat(
  "Resolved genera before prevalence filter:",
  nrow(genus_counts_all),
  "\n"
)

cat(
  "Genera passing >=10% prevalence:",
  nrow(genus_counts),
  "\n"
)

cat(
  "Genera quantitatively modelled:",
  nrow(primary),
  "\n"
)

cat(
  "Genera with >=1 structural zero:",
  nrow(structural_zero_taxa),
  "\n"
)

cat(
  "Robust global genera:",
  nrow(global_robust),
  "\n\n"
)

print(
  pairwise_summary,
  n = Inf
)

cat(
  "\nUnique robust pairwise genera:",
  n_distinct(all_pairwise_robust$taxon),
  "\n"
)

cat(
  "\nResults saved in:\n",
  normalizePath(outdir),
  "\n"
)

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")

# ============================================================
# H1 — COMPLETE ANCOM-BC2
#
# ALREADY COMPLETE:
#   Prokaryotes — Genus
#
# RUN NOW:
#   1. Prokaryotes — ASV
#   2. Fungi       — Genus
#   3. Fungi       — ASV
#
# H1:
#   Bulk soil
#   NATURAL / BESE / BARE / DEGRADED
#
# PRIMARY MODEL:
#   Site + Condition
#
# NO RAREFACTION
# RAW COUNTS
# >=10% prevalence within H1
# ============================================================


suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

set.seed(123)


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

master_outdir <- "outputs/H1_ANCOMBC2"

dir.create(
  master_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. SOFTWARE
# ============================================================

cat("\n========================================\n")
cat("SOFTWARE\n")
cat("========================================\n")

cat("R:", R.version.string, "\n")
cat(
  "ANCOMBC:",
  as.character(packageVersion("ANCOMBC")),
  "\n"
)


# ============================================================
# 3. LOAD HISTORICAL PHYLOSEQ OBJECTS
# ============================================================

ps_prok <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

ps_fun <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)


# ============================================================
# 4. H1 SUBSET
# ============================================================

subset_H1 <- function(ps) {
  
  x <- subset_samples(
    ps,
    SoilPortion == "BS" &
      Treatment %in% c(
        "NATURAL",
        "BESE",
        "BARE",
        "DEGRADED"
      )
  )
  
  x <- prune_taxa(
    taxa_sums(x) > 0,
    x
  )
  
  x
}


ps_A_prok <- subset_H1(ps_prok)
ps_A_fun  <- subset_H1(ps_fun)


cat("\n========================================\n")
cat("H1 INPUT\n")
cat("========================================\n")

cat(
  "Prokaryotes:",
  nsamples(ps_A_prok), "samples;",
  ntaxa(ps_A_prok), "detected ASVs;",
  sum(sample_sums(ps_A_prok)), "reads\n"
)

cat(
  "Fungi:",
  nsamples(ps_A_fun), "samples;",
  ntaxa(ps_A_fun), "detected ASVs;",
  sum(sample_sums(ps_A_fun)), "reads\n"
)


# ============================================================
# 5. RAW COUNTS
# samples x taxa
# ============================================================

get_raw_counts <- function(ps) {
  
  x <- as(
    otu_table(ps),
    "matrix"
  )
  
  if (taxa_are_rows(ps)) {
    x <- t(x)
  }
  
  x[
    sample_names(ps),
    taxa_names(ps),
    drop = FALSE
  ]
}


# ============================================================
# 6. METADATA
# ============================================================

make_metadata <- function(ps) {
  
  meta <- as(
    sample_data(ps),
    "data.frame"
  )
  
  meta <- meta[
    sample_names(ps),
    ,
    drop = FALSE
  ]
  
  out <- data.frame(
    SampleID = rownames(meta),
    
    Site = factor(
      meta$Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Condition = factor(
      meta$Treatment,
      levels = c(
        "NATURAL",
        "BESE",
        "BARE",
        "DEGRADED"
      )
    ),
    
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- out$SampleID
  
  out
}


meta_prok <- make_metadata(ps_A_prok)
meta_fun  <- make_metadata(ps_A_fun)


cat("\n========================================\n")
cat("H1 SAMPLE DESIGN\n")
cat("========================================\n")

cat("\nProkaryotes:\n")
print(
  table(
    meta_prok$Site,
    meta_prok$Condition
  )
)

cat("\nFungi:\n")
print(
  table(
    meta_fun$Site,
    meta_fun$Condition
  )
)


# ============================================================
# 7. PREPARE ASV COUNTS
# ============================================================

prepare_asv <- function(ps, label) {
  
  raw <- get_raw_counts(ps)
  
  prevalence <- colSums(
    raw > 0
  )
  
  threshold <- ceiling(
    0.10 * nrow(raw)
  )
  
  keep <- prevalence >= threshold
  
  raw_filtered <- raw[
    ,
    keep,
    drop = FALSE
  ]
  
  counts <- t(
    raw_filtered
  )
  
  original_asv <- rownames(counts)
  
  stable_ids <- paste0(
    "ASV",
    seq_len(nrow(counts))
  )
  
  rownames(counts) <- stable_ids
  
  
  tax <- as(
    tax_table(ps),
    "matrix"
  )
  
  tax <- tax[
    original_asv,
    ,
    drop = FALSE
  ]
  
  tax_df <- as.data.frame(
    tax,
    stringsAsFactors = FALSE
  )
  
  
  map <- tibble(
    taxon = stable_ids,
    
    Original_ASV =
      original_asv,
    
    Prevalence =
      as.integer(
        prevalence[original_asv]
      ),
    
    Total_reads =
      as.numeric(
        rowSums(counts)
      )
  ) %>%
    bind_cols(
      as_tibble(tax_df)
    )
  
  
  cat("\n========================================\n")
  cat(label, "— ASV PREPROCESSING\n")
  cat("========================================\n")
  
  cat(
    "Samples:",
    nrow(raw),
    "\n"
  )
  
  cat(
    "10% prevalence threshold:",
    threshold,
    "samples\n"
  )
  
  cat(
    "ASVs detected before prevalence filter:",
    ncol(raw),
    "\n"
  )
  
  cat(
    "ASVs entering ANCOM-BC2:",
    nrow(counts),
    "\n"
  )
  
  cat(
    "Reads retained:",
    sum(counts),
    "of",
    sum(raw),
    "(",
    round(
      100 * sum(counts) / sum(raw),
      2
    ),
    "%)\n"
  )
  
  
  list(
    counts = counts,
    map = map,
    threshold = threshold,
    taxa_before_filter = ncol(raw)
  )
}


# ============================================================
# 8. PREPARE GENUS COUNTS
# ============================================================

prepare_genus <- function(ps, label) {
  
  raw <- get_raw_counts(ps)
  
  tax <- as(
    tax_table(ps),
    "matrix"
  )
  
  tax <- tax[
    colnames(raw),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    "Genus" %in% colnames(tax)
  )
  
  genus <- tax[, "Genus"]
  
  assigned <-
    !is.na(genus) &
    trimws(genus) != ""
  
  raw_assigned <- raw[
    ,
    assigned,
    drop = FALSE
  ]
  
  genus_assigned <- trimws(
    genus[assigned]
  )
  
  
  # Aggregate raw ASV counts to genus
  counts_all <- rowsum(
    t(raw_assigned),
    group = genus_assigned,
    reorder = FALSE
  )
  
  
  prevalence <- rowSums(
    counts_all > 0
  )
  
  threshold <- ceiling(
    0.10 * ncol(counts_all)
  )
  
  keep <- prevalence >= threshold
  
  counts <- counts_all[
    keep,
    ,
    drop = FALSE
  ]
  
  
  original_genus <- rownames(counts)
  
  stable_ids <- paste0(
    "GENUS",
    seq_len(nrow(counts))
  )
  
  rownames(counts) <- stable_ids
  
  
  map <- tibble(
    taxon = stable_ids,
    
    Genus =
      original_genus,
    
    Prevalence =
      as.integer(
        prevalence[keep]
      ),
    
    Total_reads =
      as.numeric(
        rowSums(counts)
      )
  )
  
  
  cat("\n========================================\n")
  cat(label, "— GENUS PREPROCESSING\n")
  cat("========================================\n")
  
  cat(
    "Samples:",
    ncol(counts_all),
    "\n"
  )
  
  cat(
    "10% prevalence threshold:",
    threshold,
    "samples\n"
  )
  
  cat(
    "Resolved genera before filter:",
    nrow(counts_all),
    "\n"
  )
  
  cat(
    "Genera entering ANCOM-BC2:",
    nrow(counts),
    "\n"
  )
  
  cat(
    "Reads retained:",
    sum(counts),
    "of",
    sum(counts_all),
    "(",
    round(
      100 * sum(counts) / sum(counts_all),
      2
    ),
    "%)\n"
  )
  
  
  list(
    counts = counts,
    map = map,
    threshold = threshold,
    taxa_before_filter = nrow(counts_all)
  )
}


# ============================================================
# 9. SAFE ROBUST CONTRAST EXTRACTOR
#
# Forces empty tables to retain correct column types.
# ============================================================

extract_robust_contrast <- function(
    pair_table,
    map,
    term,
    contrast_label,
    positive_label,
    negative_label
) {
  
  lfc_col  <- paste0("lfc_", term)
  q_col    <- paste0("q_", term)
  diff_col <- paste0("diff_", term)
  ss_col   <- paste0("passed_ss_", term)
  
  required <- c(
    "taxon",
    lfc_col,
    q_col,
    diff_col,
    ss_col
  )
  
  if (!all(required %in% colnames(pair_table))) {
    
    stop(
      paste(
        "Missing ANCOM-BC2 columns for:",
        term
      )
    )
  }
  
  
  pair_table %>%
    
    filter(
      .data[[diff_col]] %in% TRUE,
      .data[[ss_col]] %in% TRUE
    ) %>%
    
    transmute(
      taxon =
        as.character(taxon),
      
      Contrast =
        as.character(
          contrast_label
        ),
      
      LFC =
        as.numeric(
          .data[[lfc_col]]
        ),
      
      q =
        as.numeric(
          .data[[q_col]]
        ),
      
      Higher =
        as.character(
          ifelse(
            .data[[lfc_col]] > 0,
            positive_label,
            negative_label
          )
        )
    ) %>%
    
    left_join(
      map,
      by = "taxon"
    ) %>%
    
    arrange(q) %>%
    
    as_tibble()
}


# ============================================================
# 10. GENERAL H1 ANCOM-BC2 FUNCTION
# ============================================================

run_H1_ancom <- function(
    input_object,
    meta,
    kingdom_label,
    level_label,
    folder_label
) {
  
  counts <- input_object$counts
  map    <- input_object$map
  
  
  counts <- counts[
    ,
    rownames(meta),
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      colnames(counts),
      rownames(meta)
    )
  )
  
  
  outdir <- file.path(
    master_outdir,
    folder_label
  )
  
  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ----------------------------------------------------------
  # SAVE EXACT INPUT
  # ----------------------------------------------------------
  
  saveRDS(
    counts,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_counts_10pct.rds"
      )
    )
  )
  
  write.csv(
    map,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_taxon_map.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    meta,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_metadata.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # RUN MODEL
  # ----------------------------------------------------------
  
  cat("\n\n========================================\n")
  cat(
    "RUNNING H1 —",
    toupper(kingdom_label),
    "—",
    toupper(level_label),
    "\n"
  )
  cat("========================================\n")
  
  
  set.seed(123)
  
  fit <- ANCOMBC::ancombc2(
    
    data = counts,
    
    taxa_are_rows = TRUE,
    
    meta_data = meta,
    
    fix_formula =
      "Site + Condition",
    
    rand_formula = NULL,
    
    p_adj_method =
      "BH",
    
    pseudo_sens =
      TRUE,
    
    prv_cut =
      0,
    
    lib_cut =
      0,
    
    s0_perc =
      0.05,
    
    group =
      "Condition",
    
    struc_zero =
      TRUE,
    
    neg_lb =
      FALSE,
    
    alpha =
      0.05,
    
    n_cl =
      1,
    
    verbose =
      TRUE,
    
    global =
      TRUE,
    
    pairwise =
      TRUE,
    
    dunnet =
      FALSE,
    
    trend =
      FALSE
  )
  
  
  saveRDS(
    fit,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_ANCOMBC2_full_object.rds"
      )
    )
  )
  
  
  # ==========================================================
  # PRIMARY RESULTS
  # ==========================================================
  
  primary <- as.data.frame(
    fit$res
  )
  
  primary_named <- as_tibble(
    primary
  ) %>%
    left_join(
      map,
      by = "taxon"
    )
  
  
  write.csv(
    primary_named,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_primary_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # GLOBAL RESULTS
  # ==========================================================
  
  global <- as.data.frame(
    fit$res_global
  )
  
  global_named <- as_tibble(
    global
  ) %>%
    left_join(
      map,
      by = "taxon"
    )
  
  
  global_robust <- global_named %>%
    filter(
      diff_abn %in% TRUE,
      passed_ss %in% TRUE
    ) %>%
    arrange(q_val)
  
  
  write.csv(
    global_named,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_global_all.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    global_robust,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_global_ROBUST.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # STRUCTURAL ZEROS
  # ==========================================================
  
  zero_ind <- as.data.frame(
    fit$zero_ind
  )
  
  zero_cols <- grep(
    "^structural_zero",
    colnames(zero_ind),
    value = TRUE
  )
  
  
  if (length(zero_cols) > 0) {
    
    has_structural_zero <- apply(
      zero_ind[
        ,
        zero_cols,
        drop = FALSE
      ],
      1,
      function(x) {
        any(x %in% TRUE)
      }
    )
    
  } else {
    
    has_structural_zero <- rep(
      FALSE,
      nrow(zero_ind)
    )
  }
  
  
  zero_named <- as_tibble(
    zero_ind
  ) %>%
    left_join(
      map,
      by = "taxon"
    )
  
  
  structural_zero_taxa <-
    zero_named[
      has_structural_zero,
      ,
      drop = FALSE
    ]
  
  
  write.csv(
    zero_named,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_structural_zero_all.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    structural_zero_taxa,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_structural_zero_taxa.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # PAIRWISE RESULTS
  # ==========================================================
  
  pair <- as.data.frame(
    fit$res_pair
  )
  
  
  pair_named <- as_tibble(
    pair
  ) %>%
    left_join(
      map,
      by = "taxon"
    )
  
  
  write.csv(
    pair_named,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_pairwise_all.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # SIX H1 CONTRASTS
  # ==========================================================
  
  BESE_NAT <- extract_robust_contrast(
    pair,
    map,
    "ConditionBESE",
    "BESE - NATURAL",
    "BESE",
    "NATURAL"
  )
  
  
  BARE_NAT <- extract_robust_contrast(
    pair,
    map,
    "ConditionBARE",
    "BARE - NATURAL",
    "BARE",
    "NATURAL"
  )
  
  
  DEG_NAT <- extract_robust_contrast(
    pair,
    map,
    "ConditionDEGRADED",
    "DEGRADED - NATURAL",
    "DEGRADED",
    "NATURAL"
  )
  
  
  BARE_BESE <- extract_robust_contrast(
    pair,
    map,
    "ConditionBARE_ConditionBESE",
    "BARE - BESE",
    "BARE",
    "BESE"
  )
  
  
  DEG_BESE <- extract_robust_contrast(
    pair,
    map,
    "ConditionDEGRADED_ConditionBESE",
    "DEGRADED - BESE",
    "DEGRADED",
    "BESE"
  )
  
  
  DEG_BARE <- extract_robust_contrast(
    pair,
    map,
    "ConditionDEGRADED_ConditionBARE",
    "DEGRADED - BARE",
    "DEGRADED",
    "BARE"
  )
  
  
  # ==========================================================
  # SAVE ROBUST CONTRASTS
  # ==========================================================
  
  robust_tables <- list(
    BESE_vs_NATURAL = BESE_NAT,
    BARE_vs_NATURAL = BARE_NAT,
    DEGRADED_vs_NATURAL = DEG_NAT,
    BARE_vs_BESE = BARE_BESE,
    DEGRADED_vs_BESE = DEG_BESE,
    DEGRADED_vs_BARE = DEG_BARE
  )
  
  
  for (nm in names(robust_tables)) {
    
    write.csv(
      robust_tables[[nm]],
      
      file.path(
        outdir,
        paste0(
          "H1_",
          kingdom_label,
          "_",
          level_label,
          "_",
          nm,
          "_ROBUST.csv"
        )
      ),
      
      row.names = FALSE
    )
  }
  
  
  # ==========================================================
  # COMBINED ROBUST RESULTS
  # ==========================================================
  
  all_robust <- bind_rows(
    robust_tables
  )
  
  
  write.csv(
    all_robust,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_all_ROBUST_pairwise_results.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # PAIRWISE SUMMARY
  # ==========================================================
  
  pairwise_summary <- tibble(
    Kingdom = kingdom_label,
    Level = level_label,
    
    Contrast = c(
      "BESE - NATURAL",
      "BARE - NATURAL",
      "DEGRADED - NATURAL",
      "BARE - BESE",
      "DEGRADED - BESE",
      "DEGRADED - BARE"
    ),
    
    Robust_taxa = c(
      nrow(BESE_NAT),
      nrow(BARE_NAT),
      nrow(DEG_NAT),
      nrow(BARE_BESE),
      nrow(DEG_BESE),
      nrow(DEG_BARE)
    )
  )
  
  
  # ==========================================================
  # STRUCTURAL ZERO SUMMARY
  # ==========================================================
  
  structural_zero_summary <- tibble(
    Kingdom = kingdom_label,
    Level = level_label,
    
    Condition_column =
      zero_cols,
    
    Structural_zero_classifications =
      sapply(
        zero_cols,
        function(z) {
          sum(
            zero_ind[[z]],
            na.rm = TRUE
          )
        }
      )
  )
  
  
  # ==========================================================
  # FINAL SUMMARY
  # ==========================================================
  
  summary_final <- tibble(
    Kingdom = kingdom_label,
    
    Level = level_label,
    
    Prevalence_threshold_samples =
      input_object$threshold,
    
    Taxa_before_prevalence_filter =
      input_object$taxa_before_filter,
    
    Taxa_entering_ANCOMBC2 =
      nrow(counts),
    
    Structural_zero_taxa =
      sum(has_structural_zero),
    
    Quantitatively_modelled_taxa =
      nrow(primary),
    
    BH_significant_global =
      sum(
        global$diff_abn %in% TRUE,
        na.rm = TRUE
      ),
    
    Robust_global =
      nrow(global_robust),
    
    Robust_BESE_vs_NATURAL =
      nrow(BESE_NAT),
    
    Robust_BARE_vs_NATURAL =
      nrow(BARE_NAT),
    
    Robust_DEGRADED_vs_NATURAL =
      nrow(DEG_NAT),
    
    Robust_BARE_vs_BESE =
      nrow(BARE_BESE),
    
    Robust_DEGRADED_vs_BESE =
      nrow(DEG_BESE),
    
    Robust_DEGRADED_vs_BARE =
      nrow(DEG_BARE),
    
    Unique_robust_pairwise_taxa =
      n_distinct(
        all_robust$taxon
      )
  )
  
  
  write.csv(
    pairwise_summary,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_pairwise_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    structural_zero_summary,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_structural_zero_summary.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    summary_final,
    file.path(
      outdir,
      paste0(
        "H1_",
        kingdom_label,
        "_",
        level_label,
        "_FINAL_SUMMARY.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  # ==========================================================
  # PRINT
  # ==========================================================
  
  cat("\n========================================\n")
  cat(
    kingdom_label,
    "—",
    level_label,
    "FINAL SUMMARY\n"
  )
  cat("========================================\n")
  
  print(
    summary_final,
    n = Inf,
    width = Inf
  )
  
  
  cat("\nROBUST PAIRWISE COUNTS:\n")
  
  print(
    pairwise_summary,
    n = Inf,
    width = Inf
  )
  
  
  cat("\nSTRUCTURAL ZERO CLASSIFICATIONS:\n")
  
  print(
    structural_zero_summary,
    n = Inf,
    width = Inf
  )
  
  
  cat(
    "\nStructural-zero taxa + quantitatively modelled taxa =",
    sum(has_structural_zero) +
      nrow(primary),
    "\n"
  )
  
  cat(
    "Taxa entering ANCOM-BC2 =",
    nrow(counts),
    "\n"
  )
  
  
  list(
    fit = fit,
    summary = summary_final,
    pairwise_summary = pairwise_summary,
    structural_zero_summary =
      structural_zero_summary,
    global_robust = global_robust,
    robust_tables = robust_tables,
    all_robust = all_robust
  )
}


# ============================================================
# 11. PREPARE THREE REMAINING INPUTS
# ============================================================

prok_asv_input <- prepare_asv(
  ps_A_prok,
  "PROKARYOTES"
)

fun_genus_input <- prepare_genus(
  ps_A_fun,
  "FUNGI"
)

fun_asv_input <- prepare_asv(
  ps_A_fun,
  "FUNGI"
)


# ============================================================
# 12. RUN PROKARYOTIC ASV
# ============================================================

H1_Prokaryotes_ASV <- run_H1_ancom(
  
  input_object =
    prok_asv_input,
  
  meta =
    meta_prok,
  
  kingdom_label =
    "Prokaryotes",
  
  level_label =
    "ASV",
  
  folder_label =
    "prokaryotes_ASV"
)


# ============================================================
# 13. RUN FUNGI GENUS
# ============================================================

H1_Fungi_Genus <- run_H1_ancom(
  
  input_object =
    fun_genus_input,
  
  meta =
    meta_fun,
  
  kingdom_label =
    "Fungi",
  
  level_label =
    "Genus",
  
  folder_label =
    "fungi_genus"
)


# ============================================================
# 14. RUN FUNGI ASV
# ============================================================

H1_Fungi_ASV <- run_H1_ancom(
  
  input_object =
    fun_asv_input,
  
  meta =
    meta_fun,
  
  kingdom_label =
    "Fungi",
  
  level_label =
    "ASV",
  
  folder_label =
    "fungi_ASV"
)


# ============================================================
# 15. MASTER SUMMARY OF NEW ANALYSES
# ============================================================

H1_remaining_summary <- bind_rows(
  
  H1_Prokaryotes_ASV$summary,
  
  H1_Fungi_Genus$summary,
  
  H1_Fungi_ASV$summary
)


cat("\n\n========================================\n")
cat("H1 — REMAINING ANCOM-BC2 MASTER SUMMARY\n")
cat("========================================\n")

print(
  H1_remaining_summary,
  n = Inf,
  width = Inf
)


# ============================================================
# 16. MASTER PAIRWISE SUMMARY
# ============================================================

H1_remaining_pairwise <- bind_rows(
  
  H1_Prokaryotes_ASV$pairwise_summary,
  
  H1_Fungi_Genus$pairwise_summary,
  
  H1_Fungi_ASV$pairwise_summary
)


cat("\n========================================\n")
cat("H1 — ROBUST PAIRWISE SUMMARY\n")
cat("========================================\n")

print(
  H1_remaining_pairwise,
  n = Inf,
  width = Inf
)


# ============================================================
# 17. SAVE MASTER OUTPUTS
# ============================================================

write.csv(
  H1_remaining_summary,
  
  file.path(
    master_outdir,
    "H1_remaining_ANCOMBC2_MASTER_SUMMARY.csv"
  ),
  
  row.names = FALSE
)


write.csv(
  H1_remaining_pairwise,
  
  file.path(
    master_outdir,
    "H1_remaining_ANCOMBC2_MASTER_PAIRWISE_SUMMARY.csv"
  ),
  
  row.names = FALSE
)


saveRDS(
  
  list(
    prokaryotes_ASV =
      H1_Prokaryotes_ASV,
    
    fungi_genus =
      H1_Fungi_Genus,
    
    fungi_ASV =
      H1_Fungi_ASV
  ),
  
  file.path(
    master_outdir,
    "H1_remaining_ANCOMBC2_objects.rds"
  )
)


capture.output(
  sessionInfo(),
  
  file = file.path(
    master_outdir,
    "sessionInfo_H1_remaining_ANCOMBC2.txt"
  )
)


cat("\n========================================\n")
cat("H1 ANCOM-BC2 REMAINING ANALYSES COMPLETE\n")
cat("========================================\n")

# ============================================================
# H1 — ROBUST FUNGAL GENERA
# PRINT IDENTITIES, DIRECTION, LFC AND q
#
# NO MODEL IS RERUN
# ============================================================

library(dplyr)
library(tibble)


fungi_H1_robust <- bind_rows(
  
  H1_Fungi_Genus$robust_tables$BESE_vs_NATURAL,
  
  H1_Fungi_Genus$robust_tables$BARE_vs_NATURAL,
  
  H1_Fungi_Genus$robust_tables$DEGRADED_vs_NATURAL,
  
  H1_Fungi_Genus$robust_tables$BARE_vs_BESE,
  
  H1_Fungi_Genus$robust_tables$DEGRADED_vs_BESE,
  
  H1_Fungi_Genus$robust_tables$DEGRADED_vs_BARE
  
) %>%
  
  select(
    Contrast,
    Genus,
    LFC,
    q,
    Higher,
    Prevalence,
    Total_reads
  ) %>%
  
  arrange(
    Contrast,
    q
  ) %>%
  
  as_tibble()


cat("\n========================================\n")
cat("H1 — ROBUST FUNGAL GENERA\n")
cat("========================================\n")

print(
  fungi_H1_robust,
  n = Inf,
  width = Inf
)


cat("\n========================================\n")
cat("UNIQUE ROBUST FUNGAL GENERA\n")
cat("========================================\n")

print(
  fungi_H1_robust %>%
    distinct(
      Genus
    ),
  n = Inf
)


cat(
  "\nUnique robust fungal genera:",
  n_distinct(fungi_H1_robust$Genus),
  "\n"
)

# ============================================================
# H1 — SUPPLEMENTARY FUNGAL GENUS ANCOM-BC2 TABLE
# NO ANALYSIS IS RERUN
# ============================================================

library(dplyr)

outdir <- "outputs/H1_ANCOMBC2/fungi_genus"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. FULL ROBUST TABLE
# ============================================================

fungi_ANCOM_full <- fungi_H1_robust %>%
  mutate(
    Genus = sub(
      "^g__",
      "",
      Genus
    )
  ) %>%
  select(
    Contrast,
    Genus,
    LFC,
    q,
    Higher,
    Prevalence,
    Total_reads
  )


write.csv(
  fungi_ANCOM_full,
  file.path(
    outdir,
    "Supplementary_Table_H1_Fungi_Genus_ANCOMBC2_ROBUST.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 2. MANUSCRIPT-READY VERSION
#
# More compact:
# Contrast | Genus | LFC | q | Higher
# ============================================================

fungi_ANCOM_manuscript <- fungi_ANCOM_full %>%
  transmute(
    Contrast = Contrast,
    Genus = Genus,
    
    `Log-fold change` =
      round(LFC, 2),
    
    `BH-adjusted p-value (q)` =
      signif(q, 3),
    
    `Higher estimated abundance` =
      case_when(
        Higher == "NATURAL" ~ "Natural",
        Higher == "BESE" ~ "BESE",
        Higher == "BARE" ~ "BARE",
        Higher == "DEGRADED" ~ "Degraded",
        TRUE ~ Higher
      )
  )


write.csv(
  fungi_ANCOM_manuscript,
  file.path(
    outdir,
    "Supplementary_Table_H1_Fungi_Genus_ANCOMBC2_MANUSCRIPT.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 3. PRINT FINAL TABLE
# ============================================================

cat("\n========================================\n")
cat("SUPPLEMENTARY TABLE — H1 FUNGI\n")
cat("========================================\n")

print(
  fungi_ANCOM_manuscript,
  n = Inf,
  width = Inf
)

cat("\nSaved to:\n")
cat(
  normalizePath(outdir),
  "\n"
)
