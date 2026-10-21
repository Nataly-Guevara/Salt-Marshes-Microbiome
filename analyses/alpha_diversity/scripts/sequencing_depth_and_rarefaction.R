# ============================================================
# SEQUENCING-DEPTH / NORMALIZATION DIAGNOSTICS
# STEP 1 — LOAD AND VERIFY FINAL PHYLOSEQ OBJECTS
# ============================================================

library(phyloseq)
library(dplyr)
library(tibble)
library(ggplot2)
library(car)
library(emmeans)

# Create a clean output folder
outdir <- "sequencing_depth_diagnostics"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# Select the FINAL analysis-ready phyloseq objects
# ------------------------------------------------------------

cat("\nSelect the FINAL 16S phyloseq RDS file\n")
file_16S <- file.choose()

cat("\nSelect the FINAL ITS phyloseq RDS file\n")
file_ITS <- file.choose()

ps_16S <- readRDS(file_16S)
ps_ITS <- readRDS(file_ITS)

# ------------------------------------------------------------
# Verify that they really are phyloseq objects
# ------------------------------------------------------------

stopifnot(
  inherits(ps_16S, "phyloseq"),
  inherits(ps_ITS, "phyloseq")
)

# ------------------------------------------------------------
# Basic summary
# ------------------------------------------------------------

summarise_ps <- function(ps, dataset) {
  
  depths <- sample_sums(ps)
  
  tibble(
    Dataset = dataset,
    Samples = nsamples(ps),
    ASVs = ntaxa(ps),
    Min_reads = min(depths),
    Q1_reads = as.numeric(quantile(depths, 0.25)),
    Median_reads = median(depths),
    Mean_reads = mean(depths),
    Q3_reads = as.numeric(quantile(depths, 0.75)),
    Max_reads = max(depths)
  )
}

depth_summary <- bind_rows(
  summarise_ps(ps_16S, "16S"),
  summarise_ps(ps_ITS, "ITS")
)

cat("\n========================================\n")
cat("DATASET SUMMARY\n")
cat("========================================\n")

print(
  depth_summary,
  n = Inf,
  width = Inf
)

# ------------------------------------------------------------
# Metadata columns
# ------------------------------------------------------------

cat("\n========================================\n")
cat("16S METADATA VARIABLES\n")
cat("========================================\n")
print(colnames(as.data.frame(sample_data(ps_16S))))

cat("\n========================================\n")
cat("ITS METADATA VARIABLES\n")
cat("========================================\n")
print(colnames(as.data.frame(sample_data(ps_ITS))))

# Save summary
write.csv(
  depth_summary,
  file.path(outdir, "sequencing_depth_overall_summary.csv"),
  row.names = FALSE
)

# ============================================================
# STEP 3B — SPEARMAN CORRELATIONS:
# ALPHA DIVERSITY vs LIBRARY SIZE
# ============================================================

run_depth_cor <- function(data, metric) {
  
  dat <- data %>%
    filter(
      !is.na(.data[[metric]]),
      !is.na(log10_LibrarySize)
    )
  
  test <- cor.test(
    dat[[metric]],
    dat$log10_LibrarySize,
    method = "spearman",
    exact = FALSE
  )
  
  tibble(
    Metric = metric,
    n = nrow(dat),
    Spearman_rho = unname(test$estimate),
    p_value = test$p.value
  )
}

# Overall
cor_overall <- alpha_depth %>%
  group_by(Dataset) %>%
  group_modify(
    ~ bind_rows(
      run_depth_cor(.x, "Chao1"),
      run_depth_cor(.x, "Shannon")
    )
  ) %>%
  ungroup() %>%
  mutate(
    Analysis = "All samples",
    SoilPortion = "All"
  )

# Within BS and RZ separately
cor_habitat <- alpha_depth %>%
  group_by(Dataset, SoilPortion) %>%
  group_modify(
    ~ bind_rows(
      run_depth_cor(.x, "Chao1"),
      run_depth_cor(.x, "Shannon")
    )
  ) %>%
  ungroup() %>%
  mutate(
    Analysis = "Within soil habitat"
  )

# Combine + BH correction
depth_correlations <- bind_rows(
  cor_overall,
  cor_habitat
) %>%
  select(
    Dataset,
    Analysis,
    SoilPortion,
    Metric,
    n,
    Spearman_rho,
    p_value
  ) %>%
  mutate(
    p_BH = p.adjust(
      p_value,
      method = "BH"
    )
  ) %>%
  arrange(
    Dataset,
    Metric,
    Analysis,
    SoilPortion
  )

cat("\n========================================\n")
cat("ALPHA DIVERSITY vs SEQUENCING DEPTH\n")
cat("========================================\n\n")

print(
  depth_correlations,
  n = Inf,
  width = Inf
)

write.csv(
  depth_correlations,
  file.path(
    outdir,
    "alpha_library_size_correlations.csv"
  ),
  row.names = FALSE
)

# ============================================================
# STEP 4 — DOES SEQUENCING DEPTH DIFFER AMONG
# ECOLOGICAL CONDITIONS?
#
# Separate tests for:
# H1 = Bulk soil
# H2 = Rhizosphere
# 16S and ITS separately
# ============================================================

library(dplyr)
library(car)
library(emmeans)

# ------------------------------------------------------------
# Function
# ------------------------------------------------------------

test_library_size <- function(
    data,
    dataset_name,
    soil_portion,
    treatments
) {
  
  dat <- data %>%
    filter(
      Dataset == dataset_name,
      SoilPortion == soil_portion,
      Treatment %in% treatments
    ) %>%
    mutate(
      Site = factor(Site),
      Treatment = factor(
        Treatment,
        levels = treatments
      )
    )
  
  
  # ----------------------------------------------------------
  # Full Site × Treatment model
  # ----------------------------------------------------------
  
  full_model <- lm(
    log10_LibrarySize ~ Site * Treatment,
    data = dat
  )
  
  full_anova <- car::Anova(
    full_model,
    type = 2
  )
  
  
  # ----------------------------------------------------------
  # Additive model
  # Used to estimate overall treatment differences
  # adjusting for Site
  # ----------------------------------------------------------
  
  additive_model <- lm(
    log10_LibrarySize ~ Site + Treatment,
    data = dat
  )
  
  additive_anova <- car::Anova(
    additive_model,
    type = 2
  )
  
  emm <- emmeans(
    additive_model,
    ~ Treatment
  )
  
  pairs_out <- pairs(
    emm,
    adjust = "tukey"
  )
  
  
  cat("\n\n========================================\n")
  cat(
    dataset_name,
    "|",
    soil_portion,
    "\n"
  )
  cat("========================================\n")
  
  cat("\nSample numbers:\n")
  print(
    table(
      dat$Site,
      dat$Treatment
    )
  )
  
  cat("\nFULL MODEL — Site × Treatment:\n")
  print(full_anova)
  
  cat("\nADDITIVE MODEL — Site + Treatment:\n")
  print(additive_anova)
  
  cat("\nSITE-ADJUSTED LIBRARY SIZE EMMs:\n")
  print(emm)
  
  cat("\nPAIRWISE TREATMENT COMPARISONS:\n")
  print(pairs_out)
  
  
  invisible(
    list(
      data = dat,
      full_model = full_model,
      additive_model = additive_model,
      full_anova = full_anova,
      additive_anova = additive_anova,
      emmeans = emm,
      pairs = pairs_out
    )
  )
}


# ============================================================
# H1 — BULK SOIL
# ============================================================

depth_16S_BS <- test_library_size(
  alpha_depth,
  dataset_name = "16S",
  soil_portion = "BS",
  treatments = c(
    "NATURAL",
    "BESE",
    "BARE",
    "DEGRADED"
  )
)

depth_ITS_BS <- test_library_size(
  alpha_depth,
  dataset_name = "ITS",
  soil_portion = "BS",
  treatments = c(
    "NATURAL",
    "BESE",
    "BARE",
    "DEGRADED"
  )
)


# ============================================================
# H2 — RHIZOSPHERE
# ============================================================

depth_16S_RZ <- test_library_size(
  alpha_depth,
  dataset_name = "16S",
  soil_portion = "RZ",
  treatments = c(
    "NATURAL",
    "BESE",
    "BARE"
  )
)

depth_ITS_RZ <- test_library_size(
  alpha_depth,
  dataset_name = "ITS",
  soil_portion = "RZ",
  treatments = c(
    "NATURAL",
    "BESE",
    "BARE"
  )
)

# ============================================================
# STEP 5 — ITS RHIZOSPHERE CHAO1
# DOES THE CONDITION EFFECT SURVIVE ADJUSTMENT
# FOR SEQUENCING DEPTH?
# ============================================================

library(dplyr)
library(car)
library(emmeans)

# ------------------------------------------------------------
# 1. Prepare ITS rhizosphere data
# ------------------------------------------------------------

its_rz <- alpha_depth %>%
  filter(
    Dataset == "ITS",
    SoilPortion == "RZ",
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE"
    )
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Condition = factor(
      Treatment,
      levels = c(
        "NATURAL",
        "BESE",
        "BARE"
      )
    )
  )

cat("\nSample numbers:\n")
print(
  table(
    its_rz$Site,
    its_rz$Condition
  )
)


# ============================================================
# 2. ORIGINAL MODEL
# Site + Condition
# ============================================================

m_original <- lm(
  Chao1 ~ Site + Condition,
  data = its_rz
)

cat("\n========================================\n")
cat("ORIGINAL MODEL — NO DEPTH ADJUSTMENT\n")
cat("========================================\n")

print(
  car::Anova(
    m_original,
    type = 2
  )
)

emm_original <- emmeans(
  m_original,
  ~ Condition
)

cat("\nOriginal estimated marginal means:\n")
print(emm_original)

cat("\nOriginal pairwise comparisons:\n")
print(
  pairs(
    emm_original,
    adjust = "holm"
  )
)


# ============================================================
# 3. DEPTH-ADJUSTED MODEL
# ============================================================

m_depth <- lm(
  Chao1 ~
    Site +
    Condition +
    log10_LibrarySize,
  data = its_rz
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED MODEL\n")
cat("========================================\n")

print(
  car::Anova(
    m_depth,
    type = 2
  )
)


# ============================================================
# 4. DEPTH-ADJUSTED ESTIMATED MARGINAL MEANS
# ============================================================

emm_depth <- emmeans(
  m_depth,
  ~ Condition
)

cat("\nDepth-adjusted estimated marginal means:\n")
print(emm_depth)

cat("\nDepth-adjusted pairwise comparisons:\n")

depth_pairs <- pairs(
  emm_depth,
  adjust = "holm"
)

print(depth_pairs)


# ============================================================
# 5. MODEL COEFFICIENTS
# ============================================================

cat("\n========================================\n")
cat("DEPTH-ADJUSTED MODEL COEFFICIENTS\n")
cat("========================================\n")

print(
  summary(m_depth)
)


# ============================================================
# 6. COMPARE CONDITION EFFECT DIRECTLY
# ============================================================

original_anova <- as.data.frame(
  car::Anova(
    m_original,
    type = 2
  )
)

depth_anova <- as.data.frame(
  car::Anova(
    m_depth,
    type = 2
  )
)

comparison <- tibble(
  Model = c(
    "Original",
    "Depth-adjusted"
  ),
  
  Condition_F = c(
    original_anova[
      "Condition",
      "F value"
    ],
    depth_anova[
      "Condition",
      "F value"
    ]
  ),
  
  Condition_p = c(
    original_anova[
      "Condition",
      "Pr(>F)"
    ],
    depth_anova[
      "Condition",
      "Pr(>F)"
    ]
  )
)

cat("\n========================================\n")
cat("CONDITION EFFECT — BEFORE vs AFTER DEPTH\n")
cat("========================================\n")

print(
  comparison,
  n = Inf,
  width = Inf
)


# ============================================================
# 7. SAVE
# ============================================================

write.csv(
  comparison,
  file.path(
    outdir,
    "ITS_RZ_Chao1_condition_before_after_depth.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(emm_depth),
  file.path(
    outdir,
    "ITS_RZ_Chao1_depth_adjusted_EMMs.csv"
  ),
  row.names = FALSE
)

write.csv(
  as.data.frame(depth_pairs),
  file.path(
    outdir,
    "ITS_RZ_Chao1_depth_adjusted_pairwise.csv"
  ),
  row.names = FALSE
)


# ============================================================
# STEP 6 — FINAL-SCALE TEST
# ITS RHIZOSPHERE log(Chao1)
# WITH vs WITHOUT SEQUENCING DEPTH
# ============================================================

library(car)
library(emmeans)
library(dplyr)

# ------------------------------------------------------------
# 1. Check Site × Condition interaction on log scale
#    WITHOUT depth adjustment
# ------------------------------------------------------------

m_log_add <- lm(
  log(Chao1) ~ Site + Condition,
  data = its_rz
)

m_log_int <- lm(
  log(Chao1) ~ Site * Condition,
  data = its_rz
)

cat("\n========================================\n")
cat("LOG CHAO1 — ORIGINAL INTERACTION TEST\n")
cat("========================================\n")

print(
  anova(
    m_log_add,
    m_log_int
  )
)


# ------------------------------------------------------------
# 2. Check Site × Condition interaction
#    WITH sequencing depth
# ------------------------------------------------------------

m_log_depth_add <- lm(
  log(Chao1) ~
    Site +
    Condition +
    log10_LibrarySize,
  data = its_rz
)

m_log_depth_int <- lm(
  log(Chao1) ~
    Site * Condition +
    log10_LibrarySize,
  data = its_rz
)

cat("\n========================================\n")
cat("LOG CHAO1 — DEPTH-ADJUSTED INTERACTION TEST\n")
cat("========================================\n")

print(
  anova(
    m_log_depth_add,
    m_log_depth_int
  )
)


# ------------------------------------------------------------
# 3. Original additive log model
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ORIGINAL LOG-CHAO1 MODEL\n")
cat("========================================\n")

print(
  car::Anova(
    m_log_add,
    type = 2
  )
)


# ------------------------------------------------------------
# 4. Depth-adjusted additive log model
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DEPTH-ADJUSTED LOG-CHAO1 MODEL\n")
cat("========================================\n")

print(
  car::Anova(
    m_log_depth_add,
    type = 2
  )
)


# ------------------------------------------------------------
# 5. Original EMMs and pairwise contrasts
# ------------------------------------------------------------

emm_log_original <- emmeans(
  m_log_add,
  ~ Condition
)

pairs_log_original <- pairs(
  emm_log_original,
  adjust = "holm"
)

cat("\n========================================\n")
cat("ORIGINAL LOG EMMs\n")
cat("========================================\n")

print(emm_log_original)

cat("\nBack-transformed original EMMs:\n")

print(
  summary(
    regrid(
      emm_log_original,
      transform = "response"
    )
  )
)

cat("\nOriginal pairwise contrasts:\n")
print(pairs_log_original)


# ------------------------------------------------------------
# 6. Depth-adjusted EMMs and pairwise contrasts
# ------------------------------------------------------------

emm_log_depth <- emmeans(
  m_log_depth_add,
  ~ Condition
)

pairs_log_depth <- pairs(
  emm_log_depth,
  adjust = "holm"
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED LOG EMMs\n")
cat("========================================\n")

print(emm_log_depth)

cat("\nBack-transformed depth-adjusted EMMs:\n")

print(
  summary(
    regrid(
      emm_log_depth,
      transform = "response"
    )
  )
)

cat("\nDepth-adjusted pairwise contrasts:\n")
print(pairs_log_depth)


# ------------------------------------------------------------
# 7. Compact comparison
# ------------------------------------------------------------

a_original <- as.data.frame(
  car::Anova(
    m_log_add,
    type = 2
  )
)

a_depth <- as.data.frame(
  car::Anova(
    m_log_depth_add,
    type = 2
  )
)

log_model_comparison <- tibble(
  Model = c(
    "Original log Chao1",
    "Depth-adjusted log Chao1"
  ),
  Condition_F = c(
    a_original["Condition", "F value"],
    a_depth["Condition", "F value"]
  ),
  Condition_p = c(
    a_original["Condition", "Pr(>F)"],
    a_depth["Condition", "Pr(>F)"]
  ),
  Depth_F = c(
    NA,
    a_depth["log10_LibrarySize", "F value"]
  ),
  Depth_p = c(
    NA,
    a_depth["log10_LibrarySize", "Pr(>F)"]
  )
)

cat("\n========================================\n")
cat("FINAL LOG-MODEL COMPARISON\n")
cat("========================================\n")

print(
  log_model_comparison,
  n = Inf,
  width = Inf
)

write.csv(
  log_model_comparison,
  file.path(
    outdir,
    "ITS_RZ_logChao1_before_after_depth.csv"
  ),
  row.names = FALSE
)

# ============================================================
# STEP 7 — RAREFACTION CURVES
# Reviewer-requested vegan::rarecurve() + abline()
#
# Diagnostic only — NO objects are rarefied or modified
# ============================================================

library(phyloseq)
library(vegan)
library(dplyr)
library(tibble)

# ------------------------------------------------------------
# Helper: phyloseq OTU table -> samples as rows
# ------------------------------------------------------------

get_sample_taxa_matrix <- function(ps) {
  
  otu <- as(
    otu_table(ps),
    "matrix"
  )
  
  if (taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  
  storage.mode(otu) <- "numeric"
  
  otu
}


# ============================================================
# 1. CREATE MATRICES
# ============================================================

otu_16S <- get_sample_taxa_matrix(ps_16S)
otu_ITS <- get_sample_taxa_matrix(ps_ITS)

stopifnot(
  nrow(otu_16S) == nsamples(ps_16S),
  nrow(otu_ITS) == nsamples(ps_ITS)
)


# ============================================================
# 2. COMMON DEPTH = MINIMUM LIBRARY SIZE
#
# This is the highest depth at which ALL samples
# could be retained.
# ============================================================

common_depth_16S <- min(rowSums(otu_16S))
common_depth_ITS <- min(rowSums(otu_ITS))

cat("\n========================================\n")
cat("COMMON DEPTHS RETAINING ALL SAMPLES\n")
cat("========================================\n")

cat(
  "16S:",
  common_depth_16S,
  "reads\n"
)

cat(
  "ITS:",
  common_depth_ITS,
  "reads\n"
)


# ============================================================
# 3. RAREFACTION CURVE FUNCTION
# ============================================================

make_rarecurve_plot <- function(
    otu,
    dataset,
    common_depth,
    filename
) {
  
  max_depth <- max(rowSums(otu))
  
  # Around 150 positions across the largest library
  step_size <- max(
    100,
    floor(max_depth / 150)
  )
  
  png(
    filename = file.path(
      outdir,
      filename
    ),
    width = 1800,
    height = 1400,
    res = 180
  )
  
  vegan::rarecurve(
    otu,
    step = step_size,
    label = FALSE,
    xlab = "Sequencing depth (reads)",
    ylab = "Expected observed ASVs",
    main = paste0(
      dataset,
      " rarefaction curves"
    )
  )
  
  # Reviewer-requested abline:
  # maximum depth retaining every sample
  abline(
    v = common_depth,
    lty = 2,
    lwd = 2
  )
  
  legend(
    "bottomright",
    legend = paste0(
      "All-sample common depth = ",
      format(
        common_depth,
        big.mark = ","
      ),
      " reads"
    ),
    lty = 2,
    lwd = 2,
    bty = "n"
  )
  
  dev.off()
}


# ============================================================
# 4. GENERATE 16S AND ITS CURVES
# ============================================================

make_rarecurve_plot(
  otu = otu_16S,
  dataset = "16S",
  common_depth = common_depth_16S,
  filename = "16S_rarefaction_curves_abline.png"
)

make_rarecurve_plot(
  otu = otu_ITS,
  dataset = "ITS",
  common_depth = common_depth_ITS,
  filename = "ITS_rarefaction_curves_abline.png"
)


# ============================================================
# 5. COMPARE EXPECTED RICHNESS AT COMMON DEPTH
#    WITH OBSERVED RICHNESS
#
# vegan::rarefy() calculates expected number of observed taxa
# at the specified sampling depth WITHOUT actually subsampling
# the phyloseq object.
# ============================================================

make_rarefaction_summary <- function(
    otu,
    dataset,
    common_depth
) {
  
  observed <- vegan::specnumber(
    otu
  )
  
  expected_common <- vegan::rarefy(
    otu,
    sample = common_depth
  )
  
  tibble(
    Dataset = dataset,
    SampleID = rownames(otu),
    LibrarySize = rowSums(otu),
    Observed_ASVs = observed,
    Expected_ASVs_common_depth =
      as.numeric(expected_common),
    Percent_observed_recovered =
      100 *
      as.numeric(expected_common) /
      observed
  )
}


rare_summary_samples <- bind_rows(
  
  make_rarefaction_summary(
    otu_16S,
    "16S",
    common_depth_16S
  ),
  
  make_rarefaction_summary(
    otu_ITS,
    "ITS",
    common_depth_ITS
  )
)


# ============================================================
# 6. DATASET-LEVEL SUMMARY — FIXED
# ============================================================

rare_summary_dataset <- rare_summary_samples %>%
  group_by(Dataset) %>%
  summarise(
    Samples = n(),
    
    Median_observed_ASVs =
      median(Observed_ASVs),
    
    Median_expected_ASVs_common_depth =
      median(Expected_ASVs_common_depth),
    
    Median_percent_recovered =
      median(Percent_observed_recovered),
    
    Q1_percent_recovered =
      as.numeric(
        quantile(
          Percent_observed_recovered,
          0.25
        )
      ),
    
    Q3_percent_recovered =
      as.numeric(
        quantile(
          Percent_observed_recovered,
          0.75
        )
      ),
    
    Min_percent_recovered =
      min(Percent_observed_recovered),
    
    Max_percent_recovered =
      max(Percent_observed_recovered),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    Common_depth = case_when(
      Dataset == "16S" ~ common_depth_16S,
      Dataset == "ITS" ~ common_depth_ITS
    )
  ) %>%
  
  select(
    Dataset,
    Samples,
    Common_depth,
    everything()
  )


cat("\n========================================\n")
cat("RAREFACTION SUMMARY\n")
cat("========================================\n")

print(
  rare_summary_dataset,
  n = Inf,
  width = Inf
)

# ============================================================
# 7. SAVE TABLES
# ============================================================

write.csv(
  rare_summary_samples,
  file.path(
    outdir,
    "rarefaction_expected_vs_observed_sample_level.csv"
  ),
  row.names = FALSE
)

write.csv(
  rare_summary_dataset,
  file.path(
    outdir,
    "rarefaction_expected_vs_observed_summary.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("FILES CREATED\n")
cat("========================================\n")

cat(
  file.path(
    outdir,
    "16S_rarefaction_curves_abline.png"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "ITS_rarefaction_curves_abline.png"
  ),
  "\n"
)



# ============================================================
# STEP 8 — REPEATED RAREFACTION SENSITIVITY
# ITS RHIZOSPHERE CHAO1
#
# 100 independent rarefactions to 8,698 reads
# Primary analysis remains NON-RAREFIED
# ============================================================

library(phyloseq)
library(dplyr)
library(car)
library(emmeans)
library(tibble)

# ------------------------------------------------------------
# 1. ITS rhizosphere phyloseq object
# ------------------------------------------------------------

ps_ITS_RZ <- subset_samples(
  ps_ITS,
  SoilPortion == "RZ" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE"
    )
)

ps_ITS_RZ <- prune_taxa(
  taxa_sums(ps_ITS_RZ) > 0,
  ps_ITS_RZ
)

cat("\nSamples:", nsamples(ps_ITS_RZ), "\n")
cat("Minimum reads:", min(sample_sums(ps_ITS_RZ)), "\n")


# ------------------------------------------------------------
# 2. Rarefaction depth
# Use all-dataset minimum established previously
# ------------------------------------------------------------

rare_depth <- common_depth_ITS

cat(
  "Rarefaction depth:",
  rare_depth,
  "reads\n"
)


# ------------------------------------------------------------
# 3. Function for one rarefaction replicate
# ------------------------------------------------------------

run_one_rarefaction <- function(seed) {
  
  set.seed(seed)
  
  ps_r <- rarefy_even_depth(
    ps_ITS_RZ,
    sample.size = rare_depth,
    rngseed = seed,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
  )
  
  
  # Alpha diversity
  rich <- estimate_richness(
    ps_r,
    measures = "Chao1"
  )
  
  rich$SampleID <- sub(
    "^X(?=\\d)",
    "",
    rownames(rich),
    perl = TRUE
  )
  
  
  # Metadata
  md <- as.data.frame(
    sample_data(ps_r),
    stringsAsFactors = FALSE
  )
  
  md$SampleID <- sub(
    "^X(?=\\d)",
    "",
    rownames(md),
    perl = TRUE
  )
  
  
  dat <- left_join(
    rich,
    md,
    by = "SampleID"
  ) %>%
    mutate(
      Site = factor(
        Site,
        levels = c(
          "ACP",
          "BASS",
          "TRUG"
        )
      ),
      
      Condition = factor(
        Treatment,
        levels = c(
          "NATURAL",
          "BESE",
          "BARE"
        )
      )
    )
  
  
  # Same final-scale model
  mod <- lm(
    log(Chao1) ~ Site + Condition,
    data = dat
  )
  
  a <- as.data.frame(
    car::Anova(
      mod,
      type = 2
    )
  )
  
  
  # EMM contrasts
  emm <- emmeans(
    mod,
    ~ Condition
  )
  
  pw <- as.data.frame(
    pairs(
      emm,
      adjust = "holm"
    )
  )
  
  
  # Extract condition omnibus
  omnibus <- tibble(
    Seed = seed,
    Condition_F =
      a["Condition", "F value"],
    Condition_p =
      a["Condition", "Pr(>F)"]
  )
  
  
  # Extract pairwise results
  pw_out <- pw %>%
    transmute(
      Seed = seed,
      Contrast = contrast,
      Estimate = estimate,
      p_Holm = p.value
    )
  
  
  list(
    omnibus = omnibus,
    pairwise = pw_out
  )
}


# ------------------------------------------------------------
# 4. Run 100 independent rarefactions
# ------------------------------------------------------------

cat("\nRunning 100 rarefaction replicates...\n")

rare_runs <- lapply(
  1:100,
  run_one_rarefaction
)

cat("Completed.\n")


# ------------------------------------------------------------
# 5. Combine results
# ------------------------------------------------------------

rare_omnibus <- bind_rows(
  lapply(
    rare_runs,
    `[[`,
    "omnibus"
  )
)

rare_pairwise <- bind_rows(
  lapply(
    rare_runs,
    `[[`,
    "pairwise"
  )
)


# ------------------------------------------------------------
# 6. Summarise omnibus condition effect
# ------------------------------------------------------------

rare_omnibus_summary <- rare_omnibus %>%
  summarise(
    Replicates = n(),
    
    Median_F =
      median(Condition_F),
    
    Median_p =
      median(Condition_p),
    
    Min_p =
      min(Condition_p),
    
    Max_p =
      max(Condition_p),
    
    Percent_p_below_0_05 =
      100 * mean(
        Condition_p < 0.05
      )
  )


# ------------------------------------------------------------
# 7. Summarise pairwise results
# ------------------------------------------------------------

rare_pairwise_summary <- rare_pairwise %>%
  group_by(Contrast) %>%
  summarise(
    Median_estimate =
      median(Estimate),
    
    Median_p_Holm =
      median(p_Holm),
    
    Min_p_Holm =
      min(p_Holm),
    
    Max_p_Holm =
      max(p_Holm),
    
    Percent_significant =
      100 * mean(
        p_Holm < 0.05
      ),
    
    .groups = "drop"
  )


# ------------------------------------------------------------
# 8. Print
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ITS RZ CHAO1 — REPEATED RAREFACTION\n")
cat("========================================\n")

cat("\nOMNIBUS CONDITION EFFECT:\n")

print(
  rare_omnibus_summary,
  n = Inf,
  width = Inf
)

cat("\nPAIRWISE CONDITION EFFECTS:\n")

print(
  rare_pairwise_summary,
  n = Inf,
  width = Inf
)


# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------

write.csv(
  rare_omnibus,
  file.path(
    outdir,
    "ITS_RZ_Chao1_100_rarefaction_omnibus.csv"
  ),
  row.names = FALSE
)

write.csv(
  rare_pairwise,
  file.path(
    outdir,
    "ITS_RZ_Chao1_100_rarefaction_pairwise.csv"
  ),
  row.names = FALSE
)

write.csv(
  rare_pairwise_summary,
  file.path(
    outdir,
    "ITS_RZ_Chao1_100_rarefaction_summary.csv"
  ),
  row.names = FALSE
)
# ============================================================
# STEP 11 — 16S RHIZOSPHERE CHAO1
# ORIGINAL vs SEQUENCING-DEPTH-ADJUSTED
# ============================================================

library(dplyr)
library(car)
library(emmeans)
library(tibble)

# ------------------------------------------------------------
# 1. Prepare bacterial RZ data
# ------------------------------------------------------------

bac_rz <- alpha_depth %>%
  filter(
    Dataset == "16S",
    SoilPortion == "RZ",
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE"
    )
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Condition = factor(
      Treatment,
      levels = c(
        "NATURAL",
        "BESE",
        "BARE"
      )
    )
  )

cat("\nSample numbers:\n")
print(
  table(
    bac_rz$Site,
    bac_rz$Condition
  )
)


# ============================================================
# 2. ORIGINAL LOG-CHAO1 MODEL
# ============================================================

m_bac_rz_original <- lm(
  log(Chao1) ~
    Site * Condition,
  data = bac_rz
)

cat("\n========================================\n")
cat("ORIGINAL 16S RZ LOG-CHAO1 MODEL\n")
cat("========================================\n")

print(
  car::Anova(
    m_bac_rz_original,
    type = 2
  )
)

cat("\nJOINT TESTS:\n")

print(
  emmeans::joint_tests(
    m_bac_rz_original
  )
)


# ============================================================
# 3. DEPTH-ADJUSTED LOG-CHAO1 MODEL
# ============================================================

m_bac_rz_depth <- lm(
  log(Chao1) ~
    Site * Condition +
    log10_LibrarySize,
  data = bac_rz
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED 16S RZ LOG-CHAO1 MODEL\n")
cat("========================================\n")

print(
  car::Anova(
    m_bac_rz_depth,
    type = 2
  )
)

cat("\nJOINT TESTS:\n")

print(
  emmeans::joint_tests(
    m_bac_rz_depth
  )
)


# ============================================================
# 4. DIRECTLY COMPARE INTERACTION BEFORE vs AFTER DEPTH
# ============================================================

a_orig <- as.data.frame(
  car::Anova(
    m_bac_rz_original,
    type = 2
  )
)

a_depth <- as.data.frame(
  car::Anova(
    m_bac_rz_depth,
    type = 2
  )
)

bac_rz_depth_comparison <- tibble(
  
  Model = c(
    "Original",
    "Depth-adjusted"
  ),
  
  Interaction_F = c(
    a_orig[
      "Site:Condition",
      "F value"
    ],
    
    a_depth[
      "Site:Condition",
      "F value"
    ]
  ),
  
  Interaction_p = c(
    a_orig[
      "Site:Condition",
      "Pr(>F)"
    ],
    
    a_depth[
      "Site:Condition",
      "Pr(>F)"
    ]
  ),
  
  Depth_F = c(
    NA,
    a_depth[
      "log10_LibrarySize",
      "F value"
    ]
  ),
  
  Depth_p = c(
    NA,
    a_depth[
      "log10_LibrarySize",
      "Pr(>F)"
    ]
  )
)

cat("\n========================================\n")
cat("16S RZ CHAO1 — BEFORE vs AFTER DEPTH\n")
cat("========================================\n")

print(
  bac_rz_depth_comparison,
  n = Inf,
  width = Inf
)


# ============================================================
# STEP 12 — 16S BULK-SOIL CHAO1
# ORIGINAL vs SEQUENCING-DEPTH-ADJUSTED
# ============================================================

library(dplyr)
library(car)
library(emmeans)
library(tibble)

# ------------------------------------------------------------
# 1. Prepare bacterial bulk-soil data
# ------------------------------------------------------------

bac_bs <- alpha_depth %>%
  filter(
    Dataset == "16S",
    SoilPortion == "BS",
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    )
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "NATURAL",
        "BESE",
        "BARE",
        "DEGRADED"
      )
    )
  )

cat("\nSample numbers:\n")
print(
  table(
    bac_bs$Site,
    bac_bs$Treatment
  )
)


# ============================================================
# 2. ORIGINAL MODEL
# ============================================================

m_bac_bs_original <- lm(
  Chao1 ~ Site * Treatment,
  data = bac_bs
)

cat("\n========================================\n")
cat("ORIGINAL 16S BS CHAO1 MODEL\n")
cat("========================================\n")

print(
  car::Anova(
    m_bac_bs_original,
    type = 2
  )
)

cat("\nSITE-SPECIFIC TREATMENT TESTS:\n")

print(
  emmeans::joint_tests(
    m_bac_bs_original,
    by = "Site"
  )
)


# ============================================================
# 3. DEPTH-ADJUSTED MODEL
# ============================================================

m_bac_bs_depth <- lm(
  Chao1 ~
    Site * Treatment +
    log10_LibrarySize,
  data = bac_bs
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED 16S BS CHAO1 MODEL\n")
cat("========================================\n")

print(
  car::Anova(
    m_bac_bs_depth,
    type = 2
  )
)

cat("\nSITE-SPECIFIC TREATMENT TESTS:\n")

print(
  emmeans::joint_tests(
    m_bac_bs_depth,
    by = "Site"
  )
)


# ============================================================
# 4. SITE-SPECIFIC PAIRWISE COMPARISONS AFTER DEPTH ADJUSTMENT
# ============================================================

emm_bac_bs_depth <- emmeans(
  m_bac_bs_depth,
  ~ Treatment | Site
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED SITE-SPECIFIC PAIRWISE TESTS\n")
cat("========================================\n")

print(
  pairs(
    emm_bac_bs_depth,
    adjust = "tukey"
  )
)


# ============================================================
# 5. BEFORE vs AFTER DEPTH
# ============================================================

a_orig_bs <- as.data.frame(
  car::Anova(
    m_bac_bs_original,
    type = 2
  )
)

a_depth_bs <- as.data.frame(
  car::Anova(
    m_bac_bs_depth,
    type = 2
  )
)

bac_bs_depth_comparison <- tibble(
  
  Model = c(
    "Original",
    "Depth-adjusted"
  ),
  
  Interaction_F = c(
    a_orig_bs[
      "Site:Treatment",
      "F value"
    ],
    
    a_depth_bs[
      "Site:Treatment",
      "F value"
    ]
  ),
  
  Interaction_p = c(
    a_orig_bs[
      "Site:Treatment",
      "Pr(>F)"
    ],
    
    a_depth_bs[
      "Site:Treatment",
      "Pr(>F)"
    ]
  ),
  
  Depth_F = c(
    NA,
    a_depth_bs[
      "log10_LibrarySize",
      "F value"
    ]
  ),
  
  Depth_p = c(
    NA,
    a_depth_bs[
      "log10_LibrarySize",
      "Pr(>F)"
    ]
  )
)

cat("\n========================================\n")
cat("16S BS CHAO1 — BEFORE vs AFTER DEPTH\n")
cat("========================================\n")

print(
  bac_bs_depth_comparison,
  n = Inf,
  width = Inf
)


# ============================================================
# STEP 14 — PAIRED 16S CHAO1
# ORIGINAL vs SEQUENCING-DEPTH-ADJUSTED
#
# Same verified H3 paired structure
# sqrt(Chao1)
# ============================================================

library(dplyr)
library(lme4)
library(lmerTest)
library(tibble)

# ------------------------------------------------------------
# 1. Load exact paired bacterial dataset
# ------------------------------------------------------------

prok_pair_file <-
  "outputs/H3_alpha/tables/H3_prokaryotes_paired_alpha_data.csv"

if (exists("alpha_prok_c")) {
  
  prok_pairs_depth <- as.data.frame(alpha_prok_c)
  
  cat("Using alpha_prok_c already loaded in R.\n")
  
} else if (file.exists(prok_pair_file)) {
  
  prok_pairs_depth <- read.csv(
    prok_pair_file,
    stringsAsFactors = FALSE
  )
  
  cat(
    "Loaded paired bacterial dataset from:\n",
    prok_pair_file,
    "\n"
  )
  
} else {
  
  stop(
    paste0(
      "Paired bacterial H3 dataset not found:\n",
      prok_pair_file
    )
  )
}


# ------------------------------------------------------------
# 2. Clean IDs
# ------------------------------------------------------------

prok_pairs_depth$SampleID <- sub(
  "^X(?=\\d)",
  "",
  prok_pairs_depth$SampleID,
  perl = TRUE
)


# ------------------------------------------------------------
# 3. Attach 16S sequencing depth
# ------------------------------------------------------------

depth_16S_lookup <- alpha_depth %>%
  filter(
    Dataset == "16S"
  ) %>%
  transmute(
    SampleID = sub(
      "^X(?=\\d)",
      "",
      SampleID,
      perl = TRUE
    ),
    LibrarySize,
    log10_LibrarySize
  )

prok_pairs_depth <- prok_pairs_depth %>%
  select(
    -any_of(
      c(
        "LibrarySize",
        "log10_LibrarySize"
      )
    )
  ) %>%
  left_join(
    depth_16S_lookup,
    by = "SampleID"
  ) %>%
  mutate(
    Site = factor(Site),
    Ecological_condition = factor(
      Ecological_condition,
      levels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    ),
    Soil_habitat = factor(
      Soil_habitat,
      levels = c(
        "Bulk soil",
        "Rhizosphere"
      )
    ),
    PlotID = factor(PlotID)
  )


cat("\n========================================\n")
cat("PAIRED 16S DATA CHECK\n")
cat("========================================\n")

cat(
  "Rows:",
  nrow(prok_pairs_depth),
  "\n"
)

cat(
  "Paired plots:",
  n_distinct(prok_pairs_depth$PlotID),
  "\n"
)

cat(
  "Missing library sizes:",
  sum(is.na(prok_pairs_depth$LibrarySize)),
  "\n"
)

print(
  with(
    prok_pairs_depth,
    table(
      Ecological_condition,
      Soil_habitat
    )
  )
)


# ============================================================
# 4. ORIGINAL MODEL — EXACT H3 STRUCTURE
# ============================================================

m_prok_C_orig_add <- lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition +
    Soil_habitat +
    (1 | PlotID),
  data = prok_pairs_depth,
  REML = FALSE
)

m_prok_C_orig_int <- lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition * Soil_habitat +
    (1 | PlotID),
  data = prok_pairs_depth,
  REML = FALSE
)

cat("\n========================================\n")
cat("ORIGINAL PAIRED 16S CHAO1 INTERACTION\n")
cat("========================================\n")

orig_lrt <- anova(
  m_prok_C_orig_add,
  m_prok_C_orig_int
)

print(orig_lrt)


# ============================================================
# 5. DEPTH-ADJUSTED MODEL
# ============================================================

m_prok_C_depth_add <- lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition +
    Soil_habitat +
    log10_LibrarySize +
    (1 | PlotID),
  data = prok_pairs_depth,
  REML = FALSE
)

m_prok_C_depth_int <- lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition * Soil_habitat +
    log10_LibrarySize +
    (1 | PlotID),
  data = prok_pairs_depth,
  REML = FALSE
)

cat("\n========================================\n")
cat("DEPTH-ADJUSTED PAIRED 16S CHAO1 INTERACTION\n")
cat("========================================\n")

depth_lrt <- anova(
  m_prok_C_depth_add,
  m_prok_C_depth_int
)

print(depth_lrt)


# ============================================================
# 6. DEPTH EFFECT IN FULL ADJUSTED MODEL
# ============================================================

m_prok_C_depth_REML <- lmerTest::lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition * Soil_habitat +
    log10_LibrarySize +
    (1 | PlotID),
  data = prok_pairs_depth,
  REML = TRUE
)

cat("\n========================================\n")
cat("FULL DEPTH-ADJUSTED MODEL\n")
cat("========================================\n")

print(
  anova(
    m_prok_C_depth_REML,
    type = 3,
    ddf = "Satterthwaite"
  )
)


# ============================================================
# 7. BEFORE vs AFTER
# ============================================================

prok_C_depth_comparison <- tibble(
  
  Model = c(
    "Original paired model",
    "Depth-adjusted paired model"
  ),
  
  Chisq_interaction = c(
    orig_lrt$Chisq[2],
    depth_lrt$Chisq[2]
  ),
  
  df = c(
    orig_lrt$`Chi Df`[2],
    depth_lrt$`Chi Df`[2]
  ),
  
  Interaction_p = c(
    orig_lrt$`Pr(>Chisq)`[2],
    depth_lrt$`Pr(>Chisq)`[2]
  )
)

cat("\n========================================\n")
cat("PAIRED 16S CHAO1 — BEFORE vs AFTER DEPTH\n")
cat("========================================\n")

print(
  prok_C_depth_comparison,
  n = Inf,
  width = Inf
)


