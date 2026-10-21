# ============================================================
# H3 — PAIRED ALPHA DIVERSITY
# Root filtering: does the BS -> RZ change depend on condition?
#
# Chao1:
#   log(RZ / BS)
#   0 = no change
#   >0 = higher richness in RZ
#
# Shannon:
#   RZ - BS
#   0 = no change
#   >0 = higher diversity in RZ
#
# Model:
#   paired change ~ Site + Condition
#
# Natural, BESE and BARE only
# Degraded excluded because it has no RZ
# ============================================================


# ============================================================
# H3 — FUNGI CHAO1 REPRODUCTION TEST
#
# Question:
# Does the within-plot RZ/BS richness response differ among
# Natural, BESE and BARE?
#
# Expected historical result:
# F(2,29) = 8.53, p = 0.0012
# ============================================================
# H3 — FUNGI
# LEGACY PAIRING LOGIC, USING CORRECTED ps_fun_C
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(lmerTest)
  library(lme4)
  library(emmeans)
})

# ------------------------------------------------------------
# 1. Metadata directly from corrected fungal phyloseq object
# ------------------------------------------------------------

metadata_fun_C <- data.frame(
  sample_data(ps_fun_C),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata_fun_C$SampleID <- rownames(metadata_fun_C)

# Make sure SampleID itself does NOT have artificial X
metadata_fun_C$SampleID <- sub(
  "^X(?=\\d)",
  "",
  metadata_fun_C$SampleID,
  perl = TRUE
)


# ------------------------------------------------------------
# 2. Reproduce legacy alpha-diversity construction
# ------------------------------------------------------------

alpha_fungi_c <- phyloseq::estimate_richness(
  ps_fun_C,
  measures = c("Chao1", "Shannon")
) %>%
  tibble::rownames_to_column("SampleID") %>%
  
  # THIS WAS THE IMPORTANT LEGACY STEP
  dplyr::mutate(
    SampleID = sub(
      "^X(?=\\d)",
      "",
      SampleID,
      perl = TRUE
    )
  ) %>%
  
  # Join to corrected fungal metadata
  dplyr::left_join(
    metadata_fun_C,
    by = "SampleID"
  ) %>%
  
  dplyr::filter(
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE"
    ),
    SoilPortion %in% c(
      "BS",
      "RZ"
    )
  ) %>%
  
  dplyr::mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Ecological_condition = factor(
      Treatment,
      levels = c(
        "NATURAL",
        "BESE",
        "BARE"
      ),
      labels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    ),
    
    Soil_habitat = factor(
      SoilPortion,
      levels = c(
        "BS",
        "RZ"
      ),
      labels = c(
        "Bulk soil",
        "Rhizosphere"
      )
    ),
    
    PlotID = interaction(
      Site,
      Treatment,
      Block,
      drop = TRUE
    )
  ) %>%
  
  # Legacy definition of a complete paired plot
  dplyr::group_by(PlotID) %>%
  dplyr::filter(
    dplyr::n() == 2,
    dplyr::n_distinct(Soil_habitat) == 2
  ) %>%
  dplyr::ungroup()


# ============================================================
# 3. VERIFY THAT LEGACY PAIRING NOW WORKS
# ============================================================

cat("\n========================================\n")
cat("H3 FUNGI — LEGACY PAIR CHECK\n")
cat("========================================\n")

cat(
  "Rows:",
  nrow(alpha_fungi_c),
  "\n"
)

cat(
  "Paired plots:",
  dplyr::n_distinct(alpha_fungi_c$PlotID),
  "\n"
)

cat(
  "Unique samples:",
  dplyr::n_distinct(alpha_fungi_c$SampleID),
  "\n"
)

cat("\nPairs by ecological condition:\n")

print(
  alpha_fungi_c %>%
    dplyr::distinct(
      PlotID,
      Ecological_condition
    ) %>%
    dplyr::count(
      Ecological_condition,
      name = "Paired_plots"
    )
)

cat("\nSamples by condition and habitat:\n")

print(
  alpha_fungi_c %>%
    dplyr::count(
      Ecological_condition,
      Soil_habitat
    )
)


# ============================================================
# 4. CHECK THE ACTUAL PAIRED IDS
# ============================================================

cat("\n========================================\n")
cat("PAIRED PLOTS\n")
cat("========================================\n")

pair_check <- alpha_fungi_c %>%
  dplyr::select(
    PlotID,
    Site,
    Ecological_condition,
    Block,
    Soil_habitat,
    SampleID
  ) %>%
  dplyr::arrange(
    Site,
    Ecological_condition,
    Block,
    Soil_habitat
  )

print(
  pair_check,
  n = Inf
)


# ============================================================
# 5. LEGACY H3 legacy variant 1 — FUNGI CHAO1
# ============================================================

model_chao_fun_C_add_ML <- lmerTest::lmer(
  Chao1 ~
    Site +
    Ecological_condition +
    Soil_habitat +
    (1 | PlotID),
  data = alpha_fungi_c,
  REML = FALSE
)

model_chao_fun_C_int_ML <- lmerTest::lmer(
  Chao1 ~
    Site +
    Ecological_condition * Soil_habitat +
    (1 | PlotID),
  data = alpha_fungi_c,
  REML = FALSE
)


cat("\n========================================\n")
cat("CHAO1 — ADDITIVE vs INTERACTION MODEL\n")
cat("========================================\n")

print(
  anova(
    model_chao_fun_C_add_ML,
    model_chao_fun_C_int_ML
  )
)


# ============================================================
# 6. CONDITION-SPECIFIC BS vs RZ EFFECTS
# Use interaction model because this is what is being tested
# ============================================================

emm_chao_fun_C <- emmeans::emmeans(
  model_chao_fun_C_int_ML,
  ~ Soil_habitat | Ecological_condition
)

chao_simple_fun_C <- emmeans::contrast(
  emm_chao_fun_C,
  method = "revpairwise",
  adjust = "holm"
)

cat("\n========================================\n")
cat("CHAO1 — RZ vs BS WITHIN CONDITION\n")
cat("========================================\n")

print(
  summary(
    chao_simple_fun_C,
    infer = c(TRUE, TRUE)
  )
)


anova(
  model_chao_fun_C_add_ML,
  model_chao_fun_C_int_ML
)


# ============================================================
# H3 — FUNGI — SHANNON
# Legacy paired BS-RZ analysis
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(car)
  library(emmeans)
})

# ------------------------------------------------------------
# 1. Convert paired long data to one row per plot
# ------------------------------------------------------------

shannon_fun_C_diff <- alpha_fungi_c %>%
  dplyr::select(
    PlotID,
    Site,
    Ecological_condition,
    Soil_habitat,
    Shannon
  ) %>%
  tidyr::pivot_wider(
    names_from = Soil_habitat,
    values_from = Shannon
  ) %>%
  dplyr::mutate(
    Delta_Shannon =
      Rhizosphere - `Bulk soil`
  )


# ------------------------------------------------------------
# 2. Verify paired dataset
# ------------------------------------------------------------

cat("\n========================================\n")
cat("FUNGI SHANNON — PAIR CHECK\n")
cat("========================================\n")

cat(
  "Paired plots:",
  nrow(shannon_fun_C_diff),
  "\n"
)

print(
  with(
    shannon_fun_C_diff,
    table(Ecological_condition)
  )
)


# ------------------------------------------------------------
# 3. Model:
# Does RZ-BS Shannon difference vary among conditions?
# ------------------------------------------------------------

model_shannon_fun_C_diff <- lm(
  Delta_Shannon ~
    Site +
    Ecological_condition,
  data = shannon_fun_C_diff
)


cat("\n========================================\n")
cat("FUNGI SHANNON — CONDITION EFFECT\n")
cat("========================================\n")

print(
  car::Anova(
    model_shannon_fun_C_diff,
    type = 2
  )
)


# ------------------------------------------------------------
# 4. Diagnostics
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DIAGNOSTICS\n")
cat("========================================\n")

cat(
  "Residual df:",
  df.residual(model_shannon_fun_C_diff),
  "\n"
)

cat(
  "Shapiro-Wilk p:",
  shapiro.test(
    residuals(model_shannon_fun_C_diff)
  )$p.value,
  "\n"
)


# ------------------------------------------------------------
# 5. Condition-specific estimated RZ-BS differences
# ------------------------------------------------------------

emm_shannon_fun_C <- emmeans::emmeans(
  model_shannon_fun_C_diff,
  ~ Ecological_condition
)

cat("\n========================================\n")
cat("RZ - BS SHANNON BY CONDITION\n")
cat("========================================\n")

print(
  summary(
    emm_shannon_fun_C,
    infer = c(TRUE, TRUE)
  )
)


# ------------------------------------------------------------
# 6. Test each condition against zero
# Holm adjustment
# ------------------------------------------------------------

cat("\n========================================\n")
cat("WITHIN-CONDITION TESTS AGAINST ZERO\n")
cat("========================================\n")

print(
  summary(
    test(
      emm_shannon_fun_C,
      null = 0,
      adjust = "holm"
    )
  )
)


# ------------------------------------------------------------
# 7. Pairwise differences among conditions
# ------------------------------------------------------------

cat("\n========================================\n")
cat("PAIRWISE CONDITION COMPARISONS\n")
cat("========================================\n")

print(
  summary(
    pairs(
      emm_shannon_fun_C,
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE)
  )
)

# ============================================================
# H3 — FUNGI — CHAO1
# Paired log-ratio analysis
#
# Response:
# log(RZ Chao1) - log(BS Chao1)
# = log(RZ / BS)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(car)
  library(emmeans)
})


# ------------------------------------------------------------
# 1. One row per paired plot
# ------------------------------------------------------------

fungi_C_chao_logratio <- alpha_fungi_c %>%
  dplyr::select(
    PlotID,
    Site,
    Ecological_condition,
    Soil_habitat,
    Chao1
  ) %>%
  tidyr::pivot_wider(
    names_from = Soil_habitat,
    values_from = Chao1
  ) %>%
  dplyr::mutate(
    Log_ratio_Chao1 =
      log(Rhizosphere) -
      log(`Bulk soil`),
    
    Site = factor(Site),
    
    Ecological_condition = factor(
      Ecological_condition,
      levels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    )
  )


# ------------------------------------------------------------
# 2. Verify
# ------------------------------------------------------------

cat("\n========================================\n")
cat("CHAO1 LOG-RATIO — PAIR CHECK\n")
cat("========================================\n")

cat(
  "Paired plots:",
  nrow(fungi_C_chao_logratio),
  "\n"
)

print(
  with(
    fungi_C_chao_logratio,
    table(Ecological_condition)
  )
)


# ------------------------------------------------------------
# 3. Model
# ------------------------------------------------------------

model_chao_fun_C_logratio <- lm(
  Log_ratio_Chao1 ~
    Site +
    Ecological_condition,
  data = fungi_C_chao_logratio
)


# ------------------------------------------------------------
# 4. Does paired BS→RZ response vary among conditions?
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ECOLOGICAL-CONDITION EFFECT\n")
cat("========================================\n")

print(
  car::Anova(
    model_chao_fun_C_logratio,
    type = 2
  )
)


# ------------------------------------------------------------
# 5. Estimated log RZ/BS response by condition
# ------------------------------------------------------------

emm_chao_fun_C_logratio <- emmeans::emmeans(
  model_chao_fun_C_logratio,
  ~ Ecological_condition
)

cat("\n========================================\n")
cat("LOG RZ/BS RESPONSE BY CONDITION\n")
cat("========================================\n")

print(
  summary(
    emm_chao_fun_C_logratio,
    infer = c(TRUE, TRUE)
  )
)


# ------------------------------------------------------------
# 6. Convert estimates to RZ / BS fold changes
# ------------------------------------------------------------

fold_chao_fun_C <- as.data.frame(
  summary(
    emm_chao_fun_C_logratio,
    infer = c(TRUE, TRUE)
  )
) %>%
  dplyr::mutate(
    RZ_BS_ratio = exp(emmean),
    Ratio_lower = exp(lower.CL),
    Ratio_upper = exp(upper.CL)
  )

cat("\n========================================\n")
cat("RZ / BS CHAO1 FOLD CHANGE\n")
cat("========================================\n")

print(
  fold_chao_fun_C[
    ,
    c(
      "Ecological_condition",
      "RZ_BS_ratio",
      "Ratio_lower",
      "Ratio_upper"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Planned contrast:
# Natural response versus equal-weight mean of BESE + BARE
# ------------------------------------------------------------

C2_chao_fun_final <- emmeans::contrast(
  emm_chao_fun_C_logratio,
  method = list(
    "Natural vs mean(BESE + BARE)" =
      c(
        1,
        -0.5,
        -0.5
      )
  ),
  adjust = "none"
)

cat("\n========================================\n")
cat("PLANNED CONTRAST\n")
cat("========================================\n")

print(
  summary(
    C2_chao_fun_final,
    infer = c(TRUE, TRUE)
  )
)


# ------------------------------------------------------------
# 8. Diagnostics
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DIAGNOSTICS\n")
cat("========================================\n")

cat(
  "Residual df:",
  df.residual(model_chao_fun_C_logratio),
  "\n"
)

cat(
  "Shapiro-Wilk p:",
  shapiro.test(
    residuals(model_chao_fun_C_logratio)
  )$p.value,
  "\n"
)

# ============================================================
# H3 — PROKARYOTES — CHAO1
# Paired bulk-soil and rhizosphere samples
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

# ============================================================
# 1. Metadata from corrected prokaryotic phyloseq object
# ============================================================

metadata_prok_C <- data.frame(
  sample_data(ps_prok_C),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata_prok_C$SampleID <- rownames(metadata_prok_C)

metadata_prok_C$SampleID <- sub(
  "^X(?=\\d)",
  "",
  metadata_prok_C$SampleID,
  perl = TRUE
)


# ============================================================
# 2. Reconstruct legacy paired alpha-diversity dataset
# ============================================================

alpha_prok_c <- phyloseq::estimate_richness(
  ps_prok_C,
  measures = c("Chao1", "Shannon")
) %>%
  tibble::rownames_to_column("SampleID") %>%
  
  dplyr::mutate(
    SampleID = sub(
      "^X(?=\\d)",
      "",
      SampleID,
      perl = TRUE
    )
  ) %>%
  
  dplyr::left_join(
    metadata_prok_C,
    by = "SampleID"
  ) %>%
  
  dplyr::filter(
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE"
    ),
    SoilPortion %in% c(
      "BS",
      "RZ"
    )
  ) %>%
  
  dplyr::mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Ecological_condition = factor(
      Treatment,
      levels = c(
        "NATURAL",
        "BESE",
        "BARE"
      ),
      labels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    ),
    
    Soil_habitat = factor(
      SoilPortion,
      levels = c(
        "BS",
        "RZ"
      ),
      labels = c(
        "Bulk soil",
        "Rhizosphere"
      )
    ),
    
    PlotID = interaction(
      Site,
      Treatment,
      Block,
      drop = TRUE
    )
  ) %>%
  
  dplyr::group_by(PlotID) %>%
  dplyr::filter(
    dplyr::n() == 2,
    dplyr::n_distinct(Soil_habitat) == 2
  ) %>%
  dplyr::ungroup()


# ============================================================
# 3. VERIFY PAIRING
# ============================================================

cat("\n========================================\n")
cat("PROKARYOTIC H3 — PAIR CHECK\n")
cat("========================================\n")

cat(
  "Rows:",
  nrow(alpha_prok_c),
  "\n"
)

cat(
  "Paired plots:",
  dplyr::n_distinct(alpha_prok_c$PlotID),
  "\n"
)

cat("\nPairs by condition:\n")

print(
  alpha_prok_c %>%
    dplyr::distinct(
      PlotID,
      Ecological_condition
    ) %>%
    dplyr::count(
      Ecological_condition,
      name = "Paired_plots"
    )
)


# ============================================================
# 4. CHAO1 — RAW, LOG AND SQRT DIAGNOSTICS
# ============================================================

model_chao_prok_C_raw <- lmerTest::lmer(
  Chao1 ~
    Site +
    Ecological_condition * Soil_habitat +
    (1 | PlotID),
  data = alpha_prok_c,
  REML = TRUE
)

model_chao_prok_C_log <- lmerTest::lmer(
  log(Chao1) ~
    Site +
    Ecological_condition * Soil_habitat +
    (1 | PlotID),
  data = alpha_prok_c,
  REML = TRUE
)

model_chao_prok_C_sqrt <- lmerTest::lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition * Soil_habitat +
    (1 | PlotID),
  data = alpha_prok_c,
  REML = TRUE
)


cat("\n========================================\n")
cat("CHAO1 TRANSFORMATION DIAGNOSTICS\n")
cat("========================================\n")

diagnostics_prok_chao <- data.frame(
  Model = c(
    "Raw",
    "Log",
    "Square-root"
  ),
  
  Singular = c(
    lme4::isSingular(model_chao_prok_C_raw),
    lme4::isSingular(model_chao_prok_C_log),
    lme4::isSingular(model_chao_prok_C_sqrt)
  ),
  
  Shapiro_W = c(
    shapiro.test(
      residuals(model_chao_prok_C_raw)
    )$statistic,
    
    shapiro.test(
      residuals(model_chao_prok_C_log)
    )$statistic,
    
    shapiro.test(
      residuals(model_chao_prok_C_sqrt)
    )$statistic
  ),
  
  Shapiro_p = c(
    shapiro.test(
      residuals(model_chao_prok_C_raw)
    )$p.value,
    
    shapiro.test(
      residuals(model_chao_prok_C_log)
    )$p.value,
    
    shapiro.test(
      residuals(model_chao_prok_C_sqrt)
    )$p.value
  )
)

print(
  diagnostics_prok_chao,
  row.names = FALSE
)


# ============================================================
# 5. PRIMARY TEST
# Does BS-RZ differentiation vary among conditions?
#
# Compare additive and interaction models by ML
# ============================================================

model_chao_prok_C_sqrt_add_ML <- lme4::lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition +
    Soil_habitat +
    (1 | PlotID),
  data = alpha_prok_c,
  REML = FALSE
)

model_chao_prok_C_sqrt_int_ML <- lme4::lmer(
  sqrt(Chao1) ~
    Site +
    Ecological_condition * Soil_habitat +
    (1 | PlotID),
  data = alpha_prok_c,
  REML = FALSE
)


cat("\n========================================\n")
cat("CHAO1 — CONDITION × SOIL-HABITAT TEST\n")
cat("========================================\n")

print(
  anova(
    model_chao_prok_C_sqrt_add_ML,
    model_chao_prok_C_sqrt_int_ML
  )
)

# ============================================================
# H3 — PROKARYOTES — CHAO1
# Condition-specific BS→RZ effects + planned contrast
# ============================================================

suppressPackageStartupMessages({
  library(emmeans)
})

# ------------------------------------------------------------
# 1. Estimated BS and RZ richness within each condition
# Full square-root interaction model
# ------------------------------------------------------------

emm_chao_prok_C <- emmeans::emmeans(
  model_chao_prok_C_sqrt,
  ~ Soil_habitat | Ecological_condition,
  lmer.df = "satterthwaite"
)


# ------------------------------------------------------------
# 2. RZ - BS effect within each condition
# ------------------------------------------------------------

chao_prok_C_habitat_effects <- emmeans::contrast(
  emm_chao_prok_C,
  method = "revpairwise"
)

cat("\n========================================\n")
cat("PROKARYOTIC CHAO1 — RZ vs BS\n")
cat("========================================\n")

print(
  summary(
    chao_prok_C_habitat_effects,
    infer = c(TRUE, TRUE)
  )
)


# ------------------------------------------------------------
# 3. Holm-adjusted within-condition tests
# ------------------------------------------------------------

cat("\n========================================\n")
cat("RZ vs BS — HOLM-ADJUSTED\n")
cat("========================================\n")

print(
  summary(
    update(
      chao_prok_C_habitat_effects,
      by = NULL
    ),
    infer = c(TRUE, TRUE),
    adjust = "holm"
  )
)


# ------------------------------------------------------------
# 4. Compare strength of BS→RZ response among conditions
# ------------------------------------------------------------

cat("\n========================================\n")
cat("PAIRWISE DIFFERENCES IN BS→RZ RESPONSE\n")
cat("========================================\n")

print(
  summary(
    pairs(
      update(
        chao_prok_C_habitat_effects,
        by = NULL
      ),
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE)
  )
)


# ------------------------------------------------------------
# 5. Planned contrast:
# Natural versus equal-weight mean of BESE + BARE
# ------------------------------------------------------------

chao_prok_C_planned <- emmeans::contrast(
  chao_prok_C_habitat_effects,
  method = list(
    "Natural vs mean(BESE + BARE)" =
      c(
        1,
        -0.5,
        -0.5
      )
  ),
  by = NULL,
  adjust = "none"
)

cat("\n========================================\n")
cat("PLANNED CONTRAST\n")
cat("========================================\n")

print(
  summary(
    chao_prok_C_planned,
    infer = c(TRUE, TRUE)
  )
)

# ============================================================
# H3 — PROKARYOTES — SHANNON
# Paired RZ - BS analysis
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(car)
  library(emmeans)
})

# ------------------------------------------------------------
# 1. Convert paired long data to one row per plot
# ------------------------------------------------------------

shannon_prok_C_diff <- alpha_prok_c %>%
  dplyr::select(
    PlotID,
    Site,
    Ecological_condition,
    Soil_habitat,
    SampleID,
    Shannon
  ) %>%
  tidyr::pivot_wider(
    names_from = Soil_habitat,
    values_from = c(
      SampleID,
      Shannon
    )
  ) %>%
  dplyr::mutate(
    Delta_Shannon =
      Shannon_Rhizosphere -
      `Shannon_Bulk soil`
  )


# ============================================================
# 2. VERIFY PAIRED DATA
# ============================================================

cat("\n========================================\n")
cat("PROKARYOTIC SHANNON — PAIR CHECK\n")
cat("========================================\n")

cat(
  "Paired plots:",
  nrow(shannon_prok_C_diff),
  "\n"
)

print(
  with(
    shannon_prok_C_diff,
    table(Ecological_condition)
  )
)


# ============================================================
# 3. STANDARD LINEAR MODEL
# ============================================================

model_shannon_prok_C_diff <- lm(
  Delta_Shannon ~
    Site +
    Ecological_condition,
  data = shannon_prok_C_diff
)


cat("\n========================================\n")
cat("TYPE-II ANOVA\n")
cat("========================================\n")

print(
  car::Anova(
    model_shannon_prok_C_diff,
    type = 2
  )
)


# ============================================================
# 4. DIAGNOSTICS
# ============================================================

cat("\n========================================\n")
cat("DIAGNOSTICS\n")
cat("========================================\n")

cat(
  "Residual df:",
  df.residual(model_shannon_prok_C_diff),
  "\n"
)

cat(
  "Shapiro-Wilk p:",
  shapiro.test(
    residuals(model_shannon_prok_C_diff)
  )$p.value,
  "\n"
)


# ============================================================
# 5. SITE-RESTRICTED RESIDUAL PERMUTATION TEST
#
# H0: ecological condition does not explain the
# paired RZ-BS Shannon response after accounting for Site
# ============================================================

model_shannon_prok_C_reduced <- lm(
  Delta_Shannon ~ Site,
  data = shannon_prok_C_diff
)

model_shannon_prok_C_full <- lm(
  Delta_Shannon ~
    Site +
    Ecological_condition,
  data = shannon_prok_C_diff
)

F_observed <- anova(
  model_shannon_prok_C_reduced,
  model_shannon_prok_C_full
)$F[2]

set.seed(123)

F_permuted <- replicate(
  9999,
  {
    
    permuted_residuals <- ave(
      residuals(
        model_shannon_prok_C_reduced
      ),
      shannon_prok_C_diff$Site,
      FUN = sample
    )
    
    permuted_response <-
      fitted(
        model_shannon_prok_C_reduced
      ) +
      permuted_residuals
    
    reduced_permuted <- lm(
      permuted_response ~ Site,
      data = shannon_prok_C_diff
    )
    
    full_permuted <- lm(
      permuted_response ~
        Site +
        Ecological_condition,
      data = shannon_prok_C_diff
    )
    
    anova(
      reduced_permuted,
      full_permuted
    )$F[2]
  }
)

permutation_p <-
  (
    1 +
      sum(
        F_permuted >= F_observed
      )
  ) / 10000


cat("\n========================================\n")
cat("PERMUTATION TEST\n")
cat("========================================\n")

print(
  data.frame(
    Effect =
      "Ecological condition on RZ-BS Shannon difference",
    F = F_observed,
    df1 = 2,
    df2 = df.residual(
      model_shannon_prok_C_full
    ),
    permutation_p = permutation_p
  ),
  row.names = FALSE
)


# ============================================================
# 6. ESTIMATED RZ-BS RESPONSE BY CONDITION
# ============================================================

emm_shannon_prok_C <- emmeans::emmeans(
  model_shannon_prok_C_diff,
  ~ Ecological_condition
)

cat("\n========================================\n")
cat("RZ - BS SHANNON BY CONDITION\n")
cat("========================================\n")

print(
  summary(
    emm_shannon_prok_C,
    infer = c(TRUE, TRUE)
  )
)


# ============================================================
# 7. WITHIN-CONDITION TESTS AGAINST ZERO
# Holm correction
# ============================================================

cat("\n========================================\n")
cat("WITHIN-CONDITION TESTS\n")
cat("========================================\n")

print(
  summary(
    test(
      emm_shannon_prok_C,
      null = 0,
      adjust = "holm"
    )
  )
)


# ============================================================
# 8. PAIRWISE CONDITION DIFFERENCES
# ============================================================

cat("\n========================================\n")
cat("PAIRWISE CONDITION COMPARISONS\n")
cat("========================================\n")

print(
  summary(
    pairs(
      emm_shannon_prok_C,
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE)
  )
)


# ============================================================
# 9. PLANNED CONTRAST:
# Natural vs mean(BESE + BARE)
# ============================================================

planned_shannon_prok_C <- emmeans::contrast(
  emm_shannon_prok_C,
  method = list(
    "Natural vs mean(BESE + BARE)" =
      c(
        1,
        -0.5,
        -0.5
      )
  ),
  adjust = "none"
)

cat("\n========================================\n")
cat("PLANNED CONTRAST\n")
cat("========================================\n")

print(
  summary(
    planned_shannon_prok_C,
    infer = c(TRUE, TRUE)
  )
)

###SAVING
# ============================================================
# SAVE H3 — ALPHA DIVERSITY
# Prokaryotes + Fungi
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
})

# ------------------------------------------------------------
# 1. Output folders
# ------------------------------------------------------------

dir.create(
  "outputs/H3_alpha",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_alpha/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_alpha/models",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. SAVE CORRECTED WORKING PHYLOSEQ OBJECTS
# ============================================================

saveRDS(
  ps_prok_C,
  "outputs/H3_alpha/models/ps_prok_H3_corrected.rds"
)

saveRDS(
  ps_fun_C,
  "outputs/H3_alpha/models/ps_fungi_H3_corrected.rds"
)


# ============================================================
# 3. SAVE PAIRED ALPHA DATA
# ============================================================

write.csv(
  alpha_prok_c,
  "outputs/H3_alpha/tables/H3_prokaryotes_paired_alpha_data.csv",
  row.names = FALSE
)

write.csv(
  alpha_fungi_c,
  "outputs/H3_alpha/tables/H3_fungi_paired_alpha_data.csv",
  row.names = FALSE
)

write.csv(
  shannon_prok_C_diff,
  "outputs/H3_alpha/tables/H3_prokaryotes_Shannon_paired_differences.csv",
  row.names = FALSE
)

write.csv(
  shannon_fun_C_diff,
  "outputs/H3_alpha/tables/H3_fungi_Shannon_paired_differences.csv",
  row.names = FALSE
)

write.csv(
  fungi_C_chao_logratio,
  "outputs/H3_alpha/tables/H3_fungi_Chao1_logratio_data.csv",
  row.names = FALSE
)


# ============================================================
# 4. FUNGI — CHAO1 RESULTS
# ============================================================

fungi_chao_simple_table <- as.data.frame(
  summary(
    chao_simple_fun_C,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  fungi_chao_simple_table,
  "outputs/H3_alpha/tables/H3_fungi_Chao1_RZ_vs_BS.csv",
  row.names = FALSE
)


fungi_chao_logratio_table <- as.data.frame(
  summary(
    emm_chao_fun_C_logratio,
    infer = c(TRUE, TRUE)
  )
) %>%
  mutate(
    RZ_BS_ratio = exp(emmean),
    Ratio_lower = exp(lower.CL),
    Ratio_upper = exp(upper.CL)
  )

write.csv(
  fungi_chao_logratio_table,
  "outputs/H3_alpha/tables/H3_fungi_Chao1_RZ_BS_fold_change.csv",
  row.names = FALSE
)


fungi_chao_planned_table <- as.data.frame(
  summary(
    C2_chao_fun_final,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  fungi_chao_planned_table,
  "outputs/H3_alpha/tables/H3_fungi_Chao1_planned_contrast.csv",
  row.names = FALSE
)


fungi_chao_interaction_table <- data.frame(
  Test = "Ecological condition x soil habitat",
  Chisq = 22.314,
  df = 2,
  p_value = 1.427e-05
)

write.csv(
  fungi_chao_interaction_table,
  "outputs/H3_alpha/tables/H3_fungi_Chao1_interaction_test.csv",
  row.names = FALSE
)


# ============================================================
# 5. FUNGI — SHANNON RESULTS
# ============================================================

fungi_shannon_anova_table <- as.data.frame(
  car::Anova(
    model_shannon_fun_C_diff,
    type = 2
  )
)

fungi_shannon_anova_table$Effect <-
  rownames(fungi_shannon_anova_table)

rownames(fungi_shannon_anova_table) <- NULL

write.csv(
  fungi_shannon_anova_table,
  "outputs/H3_alpha/tables/H3_fungi_Shannon_ANOVA.csv",
  row.names = FALSE
)


fungi_shannon_emm_table <- as.data.frame(
  summary(
    emm_shannon_fun_C,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  fungi_shannon_emm_table,
  "outputs/H3_alpha/tables/H3_fungi_Shannon_RZ_minus_BS.csv",
  row.names = FALSE
)


fungi_shannon_within_table <- as.data.frame(
  summary(
    test(
      emm_shannon_fun_C,
      null = 0,
      adjust = "holm"
    )
  )
)

write.csv(
  fungi_shannon_within_table,
  "outputs/H3_alpha/tables/H3_fungi_Shannon_within_condition_Holm.csv",
  row.names = FALSE
)


# ============================================================
# 6. PROKARYOTES — CHAO1 RESULTS
# ============================================================

write.csv(
  diagnostics_prok_chao,
  "outputs/H3_alpha/tables/H3_prokaryotes_Chao1_diagnostics.csv",
  row.names = FALSE
)


prok_chao_interaction_table <- data.frame(
  Test = "Ecological condition x soil habitat",
  Chisq = 5.2911,
  df = 2,
  p_value = 0.07097
)

write.csv(
  prok_chao_interaction_table,
  "outputs/H3_alpha/tables/H3_prokaryotes_Chao1_interaction_test.csv",
  row.names = FALSE
)


prok_chao_within_table <- as.data.frame(
  summary(
    update(
      chao_prok_C_habitat_effects,
      by = NULL
    ),
    infer = c(TRUE, TRUE),
    adjust = "holm"
  )
)

write.csv(
  prok_chao_within_table,
  "outputs/H3_alpha/tables/H3_prokaryotes_Chao1_RZ_vs_BS_Holm.csv",
  row.names = FALSE
)


prok_chao_planned_table <- as.data.frame(
  summary(
    chao_prok_C_planned,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  prok_chao_planned_table,
  "outputs/H3_alpha/tables/H3_prokaryotes_Chao1_planned_contrast.csv",
  row.names = FALSE
)


# ============================================================
# 7. PROKARYOTES — SHANNON RESULTS
# ============================================================

prok_shannon_anova_table <- as.data.frame(
  car::Anova(
    model_shannon_prok_C_diff,
    type = 2
  )
)

prok_shannon_anova_table$Effect <-
  rownames(prok_shannon_anova_table)

rownames(prok_shannon_anova_table) <- NULL

write.csv(
  prok_shannon_anova_table,
  "outputs/H3_alpha/tables/H3_prokaryotes_Shannon_ANOVA.csv",
  row.names = FALSE
)


prok_shannon_perm_table <- data.frame(
  Effect =
    "Ecological condition on paired RZ-BS Shannon difference",
  F = F_observed,
  df1 = 2,
  df2 = df.residual(
    model_shannon_prok_C_full
  ),
  permutation_p = permutation_p
)

write.csv(
  prok_shannon_perm_table,
  "outputs/H3_alpha/tables/H3_prokaryotes_Shannon_permutation_test.csv",
  row.names = FALSE
)


prok_shannon_emm_table <- as.data.frame(
  summary(
    emm_shannon_prok_C,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  prok_shannon_emm_table,
  "outputs/H3_alpha/tables/H3_prokaryotes_Shannon_RZ_minus_BS.csv",
  row.names = FALSE
)


prok_shannon_planned_table <- as.data.frame(
  summary(
    planned_shannon_prok_C,
    infer = c(TRUE, TRUE)
  )
)

write.csv(
  prok_shannon_planned_table,
  "outputs/H3_alpha/tables/H3_prokaryotes_Shannon_planned_contrast.csv",
  row.names = FALSE
)


# ============================================================
# 8. SAVE ALL MODEL OBJECTS TO ONE RDS FILE
# ============================================================

H3_alpha_models <- list(
  
  # Corrected phyloseq
  ps_prok_C = ps_prok_C,
  ps_fun_C = ps_fun_C,
  
  # Paired data
  alpha_prok_c = alpha_prok_c,
  alpha_fungi_c = alpha_fungi_c,
  
  # Prokaryotic Chao1
  model_chao_prok_C_raw =
    model_chao_prok_C_raw,
  
  model_chao_prok_C_log =
    model_chao_prok_C_log,
  
  model_chao_prok_C_sqrt =
    model_chao_prok_C_sqrt,
  
  model_chao_prok_C_sqrt_add_ML =
    model_chao_prok_C_sqrt_add_ML,
  
  model_chao_prok_C_sqrt_int_ML =
    model_chao_prok_C_sqrt_int_ML,
  
  chao_prok_C_habitat_effects =
    chao_prok_C_habitat_effects,
  
  chao_prok_C_planned =
    chao_prok_C_planned,
  
  # Prokaryotic Shannon
  shannon_prok_C_diff =
    shannon_prok_C_diff,
  
  model_shannon_prok_C_diff =
    model_shannon_prok_C_diff,
  
  permutation_p =
    permutation_p,
  
  planned_shannon_prok_C =
    planned_shannon_prok_C,
  
  # Fungal Chao1
  model_chao_fun_C_add_ML =
    model_chao_fun_C_add_ML,
  
  model_chao_fun_C_int_ML =
    model_chao_fun_C_int_ML,
  
  chao_simple_fun_C =
    chao_simple_fun_C,
  
  fungi_C_chao_logratio =
    fungi_C_chao_logratio,
  
  model_chao_fun_C_logratio =
    model_chao_fun_C_logratio,
  
  C2_chao_fun_final =
    C2_chao_fun_final,
  
  # Fungal Shannon
  shannon_fun_C_diff =
    shannon_fun_C_diff,
  
  model_shannon_fun_C_diff =
    model_shannon_fun_C_diff,
  
  emm_shannon_fun_C =
    emm_shannon_fun_C
)

saveRDS(
  H3_alpha_models,
  "outputs/H3_alpha/models/H3_alpha_all_models_and_results.rds"
)


# ============================================================
# 9. CONFIRM FILES
# ============================================================

cat("\n========================================\n")
cat("H3 ALPHA RESULTS SAVED\n")
cat("========================================\n")

print(
  list.files(
    "outputs/H3_alpha",
    recursive = TRUE
  )
)

#PCoA
# ============================================================
# H3 — PROKARYOTIC BETA DIVERSITY
# PCoA of paired bulk-soil and rhizosphere communities
#
# Complete pairs only:
# Natural = 14
# BESE    = 15
# BARE    = 15
# Total   = 44 paired plots / 88 samples
#
# Bray-Curtis on relative abundance
# Lingoes correction
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. CHECK OBJECTS
# ------------------------------------------------------------

stopifnot(
  exists("ps_prok_C"),
  exists("alpha_prok_c")
)

# ------------------------------------------------------------
# 2. METADATA FOR VERIFIED COMPLETE PAIRS
# ------------------------------------------------------------

pair_meta <- alpha_prok_c %>%
  dplyr::select(
    SampleID,
    Site,
    Ecological_condition,
    Soil_habitat,
    PlotID
  ) %>%
  dplyr::distinct()

cat("\n========================================\n")
cat("H3 — PROKARYOTIC PAIRED DATA\n")
cat("========================================\n")

cat("Samples:", nrow(pair_meta), "\n")
cat("Paired plots:", dplyr::n_distinct(pair_meta$PlotID), "\n\n")

print(
  pair_meta %>%
    distinct(PlotID, Site, Ecological_condition) %>%
    count(Site, Ecological_condition)
)

# Safety checks
stopifnot(
  nrow(pair_meta) == 88,
  dplyr::n_distinct(pair_meta$PlotID) == 44
)

# ------------------------------------------------------------
# 3. KEEP ONLY VERIFIED COMPLETE-PAIR SAMPLES
# ------------------------------------------------------------

missing_ids <- setdiff(
  pair_meta$SampleID,
  sample_names(ps_prok_C)
)

if (length(missing_ids) > 0) {
  stop(
    "These paired SampleIDs are missing from ps_prok_C: ",
    paste(missing_ids, collapse = ", ")
  )
}

ps_C_beta <- prune_samples(
  pair_meta$SampleID,
  ps_prok_C
)

# remove taxa absent from this subset
ps_C_beta <- prune_taxa(
  taxa_sums(ps_C_beta) > 0,
  ps_C_beta
)

cat("\nTaxa retained:", ntaxa(ps_C_beta), "\n")
cat("Samples retained:", nsamples(ps_C_beta), "\n")

# ------------------------------------------------------------
# 4. RELATIVE ABUNDANCE
# ------------------------------------------------------------

ps_C_rel <- transform_sample_counts(
  ps_C_beta,
  function(x) x / sum(x)
)

# ------------------------------------------------------------
# 5. BRAY-CURTIS DISTANCE
# ------------------------------------------------------------

bray_C_prok <- phyloseq::distance(
  ps_C_rel,
  method = "bray"
)

# ------------------------------------------------------------
# 6. PCoA — LINGOES CORRECTION
# ------------------------------------------------------------

pcoa_C_prok <- vegan::wcmdscale(
  bray_C_prok,
  eig = TRUE,
  add = "lingoes"
)

scores <- as.data.frame(
  pcoa_C_prok$points[, 1:2]
)

colnames(scores) <- c(
  "PCoA1",
  "PCoA2"
)

scores$SampleID <- rownames(scores)

# Variance explained
positive_eig <- pcoa_C_prok$eig[
  pcoa_C_prok$eig > 0
]

axis1_pct <- round(
  100 * positive_eig[1] / sum(positive_eig),
  1
)

axis2_pct <- round(
  100 * positive_eig[2] / sum(positive_eig),
  1
)

cat("\nPCoA1:", axis1_pct, "%\n")
cat("PCoA2:", axis2_pct, "%\n")
cat(
  "PCoA1 + PCoA2:",
  round(axis1_pct + axis2_pct, 1),
  "%\n"
)

# ------------------------------------------------------------
# 7. JOIN VERIFIED PAIR METADATA
# ------------------------------------------------------------

pcoa_dat <- scores %>%
  left_join(
    pair_meta,
    by = "SampleID"
  ) %>%
  mutate(
    Site = factor(
      Site,
      levels = c("ACP", "BASS", "TRUG")
    ),
    
    Ecological_condition = factor(
      Ecological_condition,
      levels = c("Natural", "BESE", "BARE")
    ),
    
    Soil_habitat = factor(
      Soil_habitat,
      levels = c("Bulk soil", "Rhizosphere")
    )
  )

stopifnot(
  !anyNA(pcoa_dat$PlotID),
  nrow(pcoa_dat) == 88
)

# ------------------------------------------------------------
# 8. COLOURS + SHAPES
# ------------------------------------------------------------

condition_cols <- c(
  "Natural" = "#9BC1BC",
  "BESE"    = "#FF5A5F",
  "BARE"    = "#FFA630"
)

habitat_shapes <- c(
  "Bulk soil"   = 21,   # circle
  "Rhizosphere" = 24    # triangle
)

# ------------------------------------------------------------
# 9. PCoA FIGURE
# ------------------------------------------------------------

p_pcoa_C_prok <- ggplot(
  pcoa_dat,
  aes(
    x = PCoA1,
    y = PCoA2
  )
) +
  
  # connect each BS-RZ pair
  geom_line(
    aes(
      group = PlotID
    ),
    colour = "grey70",
    linewidth = 0.45,
    alpha = 0.55
  ) +
  
  # samples
  geom_point(
    aes(
      fill = Ecological_condition,
      shape = Soil_habitat
    ),
    colour = "black",
    size = 3.1,
    stroke = 0.7,
    alpha = 0.90
  ) +
  
  # look at consistency across sites
  facet_wrap(
    ~ Site,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = condition_cols,
    name = "Ecological condition"
  ) +
  
  scale_shape_manual(
    values = habitat_shapes,
    name = "Soil habitat"
  ) +
  
  labs(
    x = paste0(
      "PCoA1 (",
      axis1_pct,
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      axis2_pct,
      "%)"
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text = element_text(
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 11
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      face = "bold"
    ),
    
    panel.spacing = unit(
      1.1,
      "lines"
    )
  )

print(p_pcoa_C_prok)

# ------------------------------------------------------------
# 10. SAVE
# ------------------------------------------------------------

dir.create(
  "outputs/H3_beta/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_beta/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "outputs/H3_beta/figures/H3_prokaryotic_PCoA_paired.png",
  p_pcoa_C_prok,
  width = 9.5,
  height = 4.5,
  dpi = 600
)

ggsave(
  "outputs/H3_beta/figures/H3_prokaryotic_PCoA_paired.pdf",
  p_pcoa_C_prok,
  width = 9.5,
  height = 4.5
)

write.csv(
  pcoa_dat,
  "outputs/H3_beta/tables/H3_prokaryotic_PCoA_scores.csv",
  row.names = FALSE
)

saveRDS(
  list(
    phyloseq_relative = ps_C_rel,
    bray = bray_C_prok,
    pcoa = pcoa_C_prok,
    scores = pcoa_dat
  ),
  "outputs/H3_beta/H3_prokaryotic_PCoA_objects.rds"
)

# ============================================================
# H3 — PROKARYOTIC BETA DIVERSITY
# PAIRED PERMANOVA
#
# Bray-Curtis, relative abundance
# 9999 permutations restricted within PlotID
#
# Tests:
# 1. Overall BS vs RZ effect
# 2. Ecological condition × soil habitat interaction
# 3. BS vs RZ within Natural, BESE and BARE
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(permute)
})

set.seed(123)

# ------------------------------------------------------------
# 1. CHECK REQUIRED OBJECTS
# ------------------------------------------------------------

stopifnot(
  exists("bray_C_prok"),
  exists("pcoa_dat")
)

# ------------------------------------------------------------
# 2. ALIGN METADATA EXACTLY TO BRAY-CURTIS MATRIX
# ------------------------------------------------------------

sample_order <- attr(bray_C_prok, "Labels")

meta_C <- pcoa_dat %>%
  dplyr::select(
    SampleID,
    Site,
    Ecological_condition,
    Soil_habitat,
    PlotID
  ) %>%
  dplyr::distinct() %>%
  dplyr::filter(SampleID %in% sample_order)

meta_C <- meta_C[
  match(sample_order, meta_C$SampleID),
]

rownames(meta_C) <- meta_C$SampleID

# Factors
meta_C$Site <- factor(
  meta_C$Site,
  levels = c("ACP", "BASS", "TRUG")
)

meta_C$Ecological_condition <- factor(
  meta_C$Ecological_condition,
  levels = c("Natural", "BESE", "BARE")
)

meta_C$Soil_habitat <- factor(
  meta_C$Soil_habitat,
  levels = c("Bulk soil", "Rhizosphere")
)

meta_C$PlotID <- factor(meta_C$PlotID)

# Safety checks
stopifnot(
  identical(rownames(meta_C), sample_order),
  nrow(meta_C) == 88,
  dplyr::n_distinct(meta_C$PlotID) == 44
)

cat("\n========================================\n")
cat("H3 — PAIRED PERMANOVA\n")
cat("========================================\n")
cat("Samples:", nrow(meta_C), "\n")
cat("Paired plots:", dplyr::n_distinct(meta_C$PlotID), "\n\n")

print(
  meta_C %>%
    distinct(PlotID, Site, Ecological_condition) %>%
    count(Site, Ecological_condition)
)


# ============================================================
# 3. RESTRICT PERMUTATIONS WITHIN PLOT
# ============================================================

perm_control <- permute::how(
  blocks = meta_C$PlotID,
  nperm = 9999
)

perm_matrix <- permute::shuffleSet(
  n = nrow(meta_C),
  nset = 9999,
  control = perm_control
)


# ============================================================
# 4. TEST 1 — OVERALL BS vs RZ EFFECT
#
# Site and ecological condition are included first as
# between-plot covariates.
#
# IMPORTANT:
# Because permutations occur only within PlotID, inferential
# interpretation is for Soil_habitat, NOT Site or Condition.
# ============================================================

perm_overall <- vegan::adonis2(
  bray_C_prok ~
    Site +
    Ecological_condition +
    Soil_habitat,
  data = meta_C,
  permutations = perm_matrix,
  by = "terms"
)

cat("\n\n========================================\n")
cat("1. OVERALL PAIRED BS vs RZ EFFECT\n")
cat("========================================\n")
print(perm_overall)


# ============================================================
# 5. TEST 2 — CONDITION × SOIL-HABITAT INTERACTION
#
# The interaction asks:
# Does the BS → RZ compositional shift differ among
# Natural, BESE and BARE?
# ============================================================

perm_interaction <- vegan::adonis2(
  bray_C_prok ~
    Site +
    Ecological_condition +
    Soil_habitat +
    Ecological_condition:Soil_habitat,
  data = meta_C,
  permutations = perm_matrix,
  by = "terms"
)

cat("\n\n========================================\n")
cat("2. CONDITION × SOIL-HABITAT INTERACTION\n")
cat("========================================\n")
print(perm_interaction)


# ============================================================
# 6. TEST 3 — BS vs RZ WITHIN EACH CONDITION
# ============================================================

bray_matrix <- as.matrix(bray_C_prok)

conditions <- c(
  "Natural",
  "BESE",
  "BARE"
)

within_results <- list()

for (cond in conditions) {
  
  # ------------------------------------------
  # subset metadata
  # ------------------------------------------
  
  meta_sub <- meta_C %>%
    dplyr::filter(
      Ecological_condition == cond
    ) %>%
    droplevels()
  
  ids <- rownames(meta_sub)
  
  # ------------------------------------------
  # subset Bray-Curtis matrix
  # ------------------------------------------
  
  bray_sub <- as.dist(
    bray_matrix[ids, ids]
  )
  
  # ------------------------------------------
  # paired permutations within PlotID
  # ------------------------------------------
  
  perm_sub_control <- permute::how(
    blocks = meta_sub$PlotID,
    nperm = 9999
  )
  
  perm_sub <- permute::shuffleSet(
    n = nrow(meta_sub),
    nset = 9999,
    control = perm_sub_control
  )
  
  # ------------------------------------------
  # PERMANOVA
  # Site first, then Soil_habitat
  # ------------------------------------------
  
  fit_sub <- vegan::adonis2(
    bray_sub ~
      Site +
      Soil_habitat,
    data = meta_sub,
    permutations = perm_sub,
    by = "terms"
  )
  
  within_results[[cond]] <- fit_sub
  
  cat("\n\n----------------------------------------\n")
  cat(cond, ": BULK SOIL vs RHIZOSPHERE\n")
  cat("----------------------------------------\n")
  
  print(fit_sub)
}


# ============================================================
# 7. EXTRACT WITHIN-CONDITION RESULTS
# ============================================================

within_table <- dplyr::bind_rows(
  lapply(
    conditions,
    function(cond) {
      
      tab <- as.data.frame(
        within_results[[cond]]
      )
      
      habitat_row <- tab[
        rownames(tab) == "Soil_habitat",
        ,
        drop = FALSE
      ]
      
      data.frame(
        Ecological_condition = cond,
        F = habitat_row$F,
        R2 = habitat_row$R2,
        p_raw = habitat_row$`Pr(>F)`
      )
    }
  )
)

# Holm correction across the three planned within-condition tests
within_table$p_Holm <- p.adjust(
  within_table$p_raw,
  method = "holm"
)

within_table <- within_table %>%
  mutate(
    R2_percent = 100 * R2
  )

cat("\n\n========================================\n")
cat("3. WITHIN-CONDITION SUMMARY\n")
cat("========================================\n")

print(
  within_table,
  row.names = FALSE
)


# ============================================================
# 8. EXTRACT PRIMARY RESULTS
# ============================================================

overall_df <- as.data.frame(perm_overall)

overall_habitat <- data.frame(
  Test = "Overall BS vs RZ",
  F = overall_df["Soil_habitat", "F"],
  R2 = overall_df["Soil_habitat", "R2"],
  p = overall_df["Soil_habitat", "Pr(>F)"]
)

interaction_df <- as.data.frame(perm_interaction)

interaction_result <- data.frame(
  Test = "Ecological condition × Soil habitat",
  F = interaction_df[
    "Ecological_condition:Soil_habitat",
    "F"
  ],
  R2 = interaction_df[
    "Ecological_condition:Soil_habitat",
    "R2"
  ],
  p = interaction_df[
    "Ecological_condition:Soil_habitat",
    "Pr(>F)"
  ]
)

primary_results <- bind_rows(
  overall_habitat,
  interaction_result
)

primary_results$R2_percent <- 100 * primary_results$R2

cat("\n\n========================================\n")
cat("PRIMARY H3 PERMANOVA RESULTS\n")
cat("========================================\n")

print(
  primary_results,
  row.names = FALSE
)


# ============================================================
# 9. SAVE
# ============================================================

dir.create(
  "outputs/H3_beta/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_beta/models",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  primary_results,
  "outputs/H3_beta/tables/H3_prok_PERMANOVA_primary.csv",
  row.names = FALSE
)

write.csv(
  within_table,
  "outputs/H3_beta/tables/H3_prok_PERMANOVA_within_condition.csv",
  row.names = FALSE
)

saveRDS(
  list(
    overall = perm_overall,
    interaction = perm_interaction,
    within_condition = within_results,
    primary_results = primary_results,
    within_table = within_table,
    metadata = meta_C
  ),
  "outputs/H3_beta/models/H3_prok_paired_PERMANOVA.rds"
)

cat("\nSaved to outputs/H3_beta/\n")



# ============================================================
# H3 — PROKARYOTIC BETA DIVERSITY
# PERMDISP DIAGNOSTICS
#
# Bray-Curtis distances from H3
#
# Questions:
# 1. Do the six Condition × Habitat groups differ in dispersion?
# 2. Within Natural, BESE and BARE, do BS and RZ differ
#    in multivariate dispersion?
#
# NOTE:
# PERMDISP is a diagnostic for the PERMANOVA results.
# It does not reproduce the paired PlotID design of PERMANOVA.
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
})

set.seed(123)

# ------------------------------------------------------------
# 1. CHECK OBJECTS
# ------------------------------------------------------------

stopifnot(
  exists("bray_C_prok"),
  exists("meta_C")
)

# Make sure metadata follows Bray-Curtis sample order
sample_order <- attr(bray_C_prok, "Labels")

stopifnot(
  identical(
    rownames(meta_C),
    sample_order
  )
)

# ============================================================
# 2. SIX-GROUP PERMDISP
# Condition × Soil habitat
# ============================================================

meta_C <- meta_C %>%
  mutate(
    Condition_habitat = interaction(
      Ecological_condition,
      Soil_habitat,
      sep = " | ",
      drop = TRUE
    )
  )

bd_all <- vegan::betadisper(
  bray_C_prok,
  group = meta_C$Condition_habitat,
  type = "median",
  bias.adjust = TRUE
)

perm_bd_all <- vegan::permutest(
  bd_all,
  permutations = 9999
)

cat("\n========================================\n")
cat("H3 — GLOBAL PERMDISP\n")
cat("Condition × Soil habitat groups\n")
cat("========================================\n")

print(perm_bd_all)

cat("\nMean distance to group centroid:\n")

dispersion_means <- data.frame(
  Condition_habitat = names(
    tapply(
      bd_all$distances,
      bd_all$group,
      mean
    )
  ),
  Mean_distance = as.numeric(
    tapply(
      bd_all$distances,
      bd_all$group,
      mean
    )
  )
)

print(
  dispersion_means,
  row.names = FALSE
)


# ============================================================
# 3. WITHIN-CONDITION PERMDISP
# BS vs RZ separately in Natural / BESE / BARE
# ============================================================

conditions <- c(
  "Natural",
  "BESE",
  "BARE"
)

bray_matrix <- as.matrix(
  bray_C_prok
)

permdisp_results <- list()
permdisp_summary <- list()

for (cond in conditions) {
  
  # ----------------------------------------------------------
  # Subset metadata
  # ----------------------------------------------------------
  
  meta_sub <- meta_C %>%
    filter(
      Ecological_condition == cond
    ) %>%
    droplevels()
  
  ids <- rownames(meta_sub)
  
  # ----------------------------------------------------------
  # Subset Bray-Curtis matrix
  # ----------------------------------------------------------
  
  bray_sub <- as.dist(
    bray_matrix[
      ids,
      ids
    ]
  )
  
  # ----------------------------------------------------------
  # PERMDISP
  # ----------------------------------------------------------
  
  bd_sub <- vegan::betadisper(
    bray_sub,
    group = meta_sub$Soil_habitat,
    type = "median",
    bias.adjust = TRUE
  )
  
  perm_sub <- vegan::permutest(
    bd_sub,
    permutations = 9999
  )
  
  permdisp_results[[cond]] <- list(
    betadisper = bd_sub,
    permutest = perm_sub
  )
  
  # ----------------------------------------------------------
  # Extract result
  # ----------------------------------------------------------
  
  perm_tab <- as.data.frame(
    perm_sub$tab
  )
  
  p_value <- perm_tab[
    1,
    "Pr(>F)"
  ]
  
  F_value <- perm_tab[
    1,
    "F"
  ]
  
  mean_disp <- tapply(
    bd_sub$distances,
    bd_sub$group,
    mean
  )
  
  permdisp_summary[[cond]] <- data.frame(
    Ecological_condition = cond,
    F = F_value,
    p_raw = p_value,
    BS_mean_dispersion =
      unname(mean_disp["Bulk soil"]),
    RZ_mean_dispersion =
      unname(mean_disp["Rhizosphere"])
  )
  
  cat("\n\n----------------------------------------\n")
  cat(cond, ": BS vs RZ PERMDISP\n")
  cat("----------------------------------------\n")
  
  print(perm_sub)
  
  cat("\nMean distance to centroid:\n")
  print(mean_disp)
}


# ============================================================
# 4. COMBINE RESULTS
# ============================================================

permdisp_table <- bind_rows(
  permdisp_summary
)

# Holm correction across the three within-condition tests
permdisp_table$p_Holm <- p.adjust(
  permdisp_table$p_raw,
  method = "holm"
)

permdisp_table <- permdisp_table %>%
  mutate(
    Difference_RZ_minus_BS =
      RZ_mean_dispersion -
      BS_mean_dispersion
  )

cat("\n\n========================================\n")
cat("H3 — WITHIN-CONDITION PERMDISP SUMMARY\n")
cat("========================================\n")

print(
  permdisp_table,
  row.names = FALSE
)


# ============================================================
# 5. SAVE
# ============================================================

dir.create(
  "outputs/H3_beta/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_beta/models",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  dispersion_means,
  "outputs/H3_beta/tables/H3_prok_PERMDISP_six_groups.csv",
  row.names = FALSE
)

write.csv(
  permdisp_table,
  "outputs/H3_beta/tables/H3_prok_PERMDISP_within_condition.csv",
  row.names = FALSE
)

saveRDS(
  list(
    global_betadisper = bd_all,
    global_permutest = perm_bd_all,
    within_condition = permdisp_results,
    summary = permdisp_table
  ),
  "outputs/H3_beta/models/H3_prok_PERMDISP.rds"
)

cat("\nSaved to outputs/H3_beta/\n")


# ============================================================
# H3 — PROKARYOTIC BETA DIVERSITY
# MAGNITUDE OF WITHIN-PLOT BS → RZ DIFFERENTIATION
#
# Response:
# Bray-Curtis distance between paired BS and RZ samples
#
# Question:
# Does the magnitude of BS-RZ differentiation differ among
# Natural, BESE and BARE after accounting for Site?
#
# n = 44 paired plots
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
  library(ggplot2)
})

set.seed(123)

# ------------------------------------------------------------
# 1. CHECK REQUIRED OBJECTS
# ------------------------------------------------------------

stopifnot(
  exists("bray_C_prok"),
  exists("meta_C")
)

bray_matrix <- as.matrix(bray_C_prok)


# ============================================================
# 2. CALCULATE ONE BS-RZ BRAY DISTANCE PER PLOT
# ============================================================

plot_ids <- unique(meta_C$PlotID)

paired_bray_list <- lapply(
  plot_ids,
  function(pid) {
    
    tmp <- meta_C %>%
      filter(PlotID == pid)
    
    # Safety check: exactly one BS and one RZ
    if (
      nrow(tmp) != 2 ||
      !all(c("Bulk soil", "Rhizosphere") %in% tmp$Soil_habitat)
    ) {
      stop(
        "Plot ", pid,
        " does not contain exactly one BS and one RZ sample."
      )
    }
    
    bs_id <- rownames(tmp)[
      tmp$Soil_habitat == "Bulk soil"
    ]
    
    rz_id <- rownames(tmp)[
      tmp$Soil_habitat == "Rhizosphere"
    ]
    
    data.frame(
      PlotID = as.character(pid),
      Site = as.character(tmp$Site[1]),
      Ecological_condition =
        as.character(tmp$Ecological_condition[1]),
      
      Bray_BS_RZ =
        bray_matrix[bs_id, rz_id]
    )
  }
)

paired_bray <- bind_rows(
  paired_bray_list
)


# ============================================================
# 3. FACTORS
# ============================================================

paired_bray <- paired_bray %>%
  mutate(
    Site = factor(
      Site,
      levels = c("ACP", "BASS", "TRUG")
    ),
    
    Ecological_condition = factor(
      Ecological_condition,
      levels = c("Natural", "BESE", "BARE")
    )
  )


# ============================================================
# 4. VERIFY DESIGN
# ============================================================

cat("\n========================================\n")
cat("H3 — PAIRED BS-RZ BRAY DISTANCE\n")
cat("========================================\n")

cat(
  "Total paired plots:",
  nrow(paired_bray),
  "\n\n"
)

print(
  paired_bray %>%
    count(
      Site,
      Ecological_condition
    )
)

cat("\nMean raw Bray distance by condition:\n")

print(
  paired_bray %>%
    group_by(Ecological_condition) %>%
    summarise(
      n = n(),
      Mean = mean(Bray_BS_RZ),
      SD = sd(Bray_BS_RZ),
      Median = median(Bray_BS_RZ),
      .groups = "drop"
    )
)

stopifnot(
  nrow(paired_bray) == 44
)


# ============================================================
# 5. ADDITIVE MODEL
# Site-adjusted ecological-condition effect
# ============================================================

mod_bray_magnitude <- lm(
  Bray_BS_RZ ~
    Site +
    Ecological_condition,
  data = paired_bray
)

cat("\n\n========================================\n")
cat("LINEAR MODEL\n")
cat("Bray BS-RZ ~ Site + Ecological condition\n")
cat("========================================\n")

print(
  summary(mod_bray_magnitude)
)


# ============================================================
# 6. PARTIAL F TESTS
# drop1 gives each term adjusted for the other
# ============================================================

anova_bray <- drop1(
  mod_bray_magnitude,
  test = "F"
)

cat("\n\n========================================\n")
cat("ADJUSTED TERM TESTS\n")
cat("========================================\n")

print(anova_bray)


# ============================================================
# 7. RESIDUAL DIAGNOSTIC
# ============================================================

shapiro_bray <- shapiro.test(
  residuals(mod_bray_magnitude)
)

cat("\n\n========================================\n")
cat("RESIDUAL NORMALITY\n")
cat("========================================\n")

print(shapiro_bray)


# ============================================================
# 8. SITE-RESTRICTED RESIDUAL PERMUTATION TEST
#
# Freedman-Lane style:
# reduced model = Site only
# residuals permuted within Site
#
# Gives robust test of Ecological_condition after Site
# ============================================================

mod_reduced <- lm(
  Bray_BS_RZ ~ Site,
  data = paired_bray
)

obs_F <- drop1(
  mod_bray_magnitude,
  test = "F"
)["Ecological_condition", "F value"]

reduced_fitted <- fitted(
  mod_reduced
)

reduced_resid <- residuals(
  mod_reduced
)

nperm <- 9999

perm_F <- numeric(
  nperm
)

for (i in seq_len(nperm)) {
  
  perm_resid <- reduced_resid
  
  # Permute only within Site
  for (s in levels(paired_bray$Site)) {
    
    idx <- which(
      paired_bray$Site == s
    )
    
    perm_resid[idx] <- sample(
      reduced_resid[idx],
      length(idx),
      replace = FALSE
    )
  }
  
  y_perm <-
    reduced_fitted +
    perm_resid
  
  perm_dat <- paired_bray
  perm_dat$Bray_perm <- y_perm
  
  fit_perm <- lm(
    Bray_perm ~
      Site +
      Ecological_condition,
    data = perm_dat
  )
  
  perm_F[i] <- drop1(
    fit_perm,
    test = "F"
  )[
    "Ecological_condition",
    "F value"
  ]
}

p_perm <- (
  1 +
    sum(perm_F >= obs_F)
) / (
  nperm + 1
)

cat("\n\n========================================\n")
cat("SITE-RESTRICTED PERMUTATION TEST\n")
cat("========================================\n")

cat(
  "Observed F =",
  round(obs_F, 4),
  "\n"
)

cat(
  "Permutation p =",
  format(
    p_perm,
    digits = 4
  ),
  "\n"
)


# ============================================================
# 9. SITE-ADJUSTED ESTIMATED MEANS
# ============================================================

emm_bray <- emmeans(
  mod_bray_magnitude,
  ~ Ecological_condition
)

emm_bray_df <- as.data.frame(
  summary(
    emm_bray,
    infer = c(TRUE, FALSE)
  )
)

cat("\n\n========================================\n")
cat("SITE-ADJUSTED ESTIMATED MEANS\n")
cat("========================================\n")

print(
  emm_bray_df,
  row.names = FALSE
)


# ============================================================
# 10. PAIRWISE CONDITION CONTRASTS
# Holm adjustment
# ============================================================

pairwise_bray <- as.data.frame(
  summary(
    pairs(
      emm_bray,
      adjust = "holm"
    ),
    infer = c(TRUE, TRUE)
  )
)

cat("\n\n========================================\n")
cat("PAIRWISE CONDITION COMPARISONS\n")
cat("Holm adjusted\n")
cat("========================================\n")

print(
  pairwise_bray,
  row.names = FALSE
)


# ============================================================
# 11. PLANNED CONTRAST
# Natural vs mean(BESE + BARE)
# ============================================================

planned_bray <- contrast(
  emm_bray,
  method = list(
    "Natural vs mean(BESE + BARE)" =
      c(
        Natural = 1,
        BESE = -0.5,
        BARE = -0.5
      )
  )
)

planned_bray_df <- as.data.frame(
  summary(
    planned_bray,
    infer = c(TRUE, TRUE)
  )
)

cat("\n\n========================================\n")
cat("PLANNED CONTRAST\n")
cat("Natural vs mean(BESE + BARE)\n")
cat("========================================\n")

print(
  planned_bray_df,
  row.names = FALSE
)


# ============================================================
# 12. SAVE RESULTS
# ============================================================

dir.create(
  "outputs/H3_beta/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_beta/models",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  paired_bray,
  "outputs/H3_beta/tables/H3_prok_paired_BS_RZ_Bray_distances.csv",
  row.names = FALSE
)

write.csv(
  emm_bray_df,
  "outputs/H3_beta/tables/H3_prok_Bray_magnitude_EMMs.csv",
  row.names = FALSE
)

write.csv(
  pairwise_bray,
  "outputs/H3_beta/tables/H3_prok_Bray_magnitude_pairwise.csv",
  row.names = FALSE
)

write.csv(
  planned_bray_df,
  "outputs/H3_beta/tables/H3_prok_Bray_magnitude_planned_contrast.csv",
  row.names = FALSE
)

saveRDS(
  list(
    data = paired_bray,
    model = mod_bray_magnitude,
    anova = anova_bray,
    shapiro = shapiro_bray,
    observed_F = obs_F,
    permutation_p = p_perm,
    emmeans = emm_bray,
    pairwise = pairwise_bray,
    planned = planned_bray
  ),
  "outputs/H3_beta/models/H3_prok_Bray_magnitude_analysis.rds"
)

cat("\n========================================\n")
cat("RESULTS SAVED\n")
cat("========================================\n")


# ============================================================
# H3 — PROKARYOTIC ANCOM-BC2
# STEP 1: INPUT + TAXONOMY PREFLIGHT
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(ANCOMBC)
})

# ------------------------------------------------------------
# 1. CHECK OBJECTS
# ------------------------------------------------------------

stopifnot(
  exists("ps_C_beta"),
  exists("meta_C")
)

cat("\n========================================\n")
cat("H3 ANCOM-BC2 — PREFLIGHT\n")
cat("========================================\n")

cat(
  "ANCOMBC version:",
  as.character(packageVersion("ANCOMBC")),
  "\n"
)

cat(
  "Samples:",
  nsamples(ps_C_beta),
  "\n"
)

cat(
  "ASVs:",
  ntaxa(ps_C_beta),
  "\n"
)

cat(
  "Total reads:",
  sum(sample_sums(ps_C_beta)),
  "\n"
)


# ============================================================
# 2. VERIFY RAW COUNTS
# ============================================================

otu_C <- as(
  otu_table(ps_C_beta),
  "matrix"
)

if (!taxa_are_rows(ps_C_beta)) {
  otu_C <- t(otu_C)
}

cat("\nRaw-count check:\n")
cat(
  "All counts non-negative:",
  all(otu_C >= 0),
  "\n"
)

cat(
  "All counts integer-like:",
  all(abs(otu_C - round(otu_C)) < 1e-8),
  "\n"
)


# ============================================================
# 3. ADD VERIFIED H3 METADATA
# ============================================================

meta_ancom <- meta_C[
  sample_names(ps_C_beta),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(meta_ancom),
    sample_names(ps_C_beta)
  )
)

sample_data(ps_C_beta)$Site <-
  meta_ancom$Site

sample_data(ps_C_beta)$Ecological_condition <-
  meta_ancom$Ecological_condition

sample_data(ps_C_beta)$Soil_habitat <-
  meta_ancom$Soil_habitat

sample_data(ps_C_beta)$PlotID <-
  meta_ancom$PlotID


# ============================================================
# 4. VERIFY PAIRED DESIGN
# ============================================================

design_check <- meta_ancom %>%
  as.data.frame() %>%
  mutate(
    PlotID = as.character(PlotID)
  ) %>%
  group_by(
    Site,
    Ecological_condition,
    PlotID
  ) %>%
  summarise(
    n_samples = n(),
    n_habitats = n_distinct(Soil_habitat),
    .groups = "drop"
  )

cat("\n========================================\n")
cat("PAIRING CHECK\n")
cat("========================================\n")

print(
  design_check %>%
    count(
      Site,
      Ecological_condition,
      n_samples,
      n_habitats
    )
)

cat(
  "\nComplete paired plots:",
  sum(
    design_check$n_samples == 2 &
      design_check$n_habitats == 2
  ),
  "\n"
)


# ============================================================
# 5. SAMPLE COUNTS BY CONDITION × HABITAT
# ============================================================

cat("\n========================================\n")
cat("SAMPLE COUNTS\n")
cat("========================================\n")

print(
  meta_ancom %>%
    as.data.frame() %>%
    count(
      Site,
      Ecological_condition,
      Soil_habitat
    )
)


# ============================================================
# 6. TAXONOMY AUDIT
# ============================================================

tax_C <- as.data.frame(
  tax_table(ps_C_beta),
  stringsAsFactors = FALSE
)

cat("\n========================================\n")
cat("TAXONOMY RANKS\n")
cat("========================================\n")

print(
  colnames(tax_C)
)

if ("Genus" %in% colnames(tax_C)) {
  
  genus_clean <- trimws(
    as.character(tax_C$Genus)
  )
  
  missing_genus <- is.na(genus_clean) |
    genus_clean == "" |
    genus_clean %in% c(
      "NA",
      "Unclassified",
      "unclassified"
    )
  
  cat(
    "\nASVs with genus assignment:",
    sum(!missing_genus),
    "\n"
  )
  
  cat(
    "ASVs without genus assignment:",
    sum(missing_genus),
    "\n"
  )
  
  cat(
    "Unique named genera:",
    length(
      unique(
        genus_clean[!missing_genus]
      )
    ),
    "\n"
  )
  
} else {
  
  stop(
    "No Genus column found in taxonomy table."
  )
}


# ============================================================
# 7. H3 COUNTS BY ECOLOGICAL CONDITION
# ============================================================

cat("\n========================================\n")
cat("H3 CONDITION TOTALS\n")
cat("========================================\n")

print(
  meta_ancom %>%
    as.data.frame() %>%
    distinct(
      PlotID,
      Site,
      Ecological_condition
    ) %>%
    count(
      Ecological_condition,
      name = "Paired_plots"
    )
)

cat("\nPreflight complete.\n")

# ============================================================
# H3 — PROKARYOTIC ANCOM-BC2
# GENUS LEVEL — PAIRED BS vs RZ
#
# Separate models within:
#   Natural
#   BESE
#   BARE
#
# Fixed effects:
#   Site + Soil_habitat
#
# Random effect:
#   (1 | PlotID)
#
# Interpretation:
#   Positive LFC = higher estimated abundance in Rhizosphere
#   Negative LFC = higher estimated abundance in Bulk soil
#
# Manual prevalence filter:
#   detected in >=10% of samples within each condition
#
# Structural zeros kept separate from quantitative results
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
  library(tibble)
})

set.seed(123)

# ============================================================
# 1. START FROM VERIFIED H3 PHYLOSEQ OBJECT
# ============================================================

stopifnot(
  exists("ps_C_beta")
)

ps_genus_base <- ps_C_beta


# ============================================================
# 2. REMOVE ASVs WITHOUT GENUS ASSIGNMENT
# ============================================================

tax0 <- as.data.frame(
  tax_table(ps_genus_base),
  stringsAsFactors = FALSE
)

genus0 <- trimws(
  as.character(tax0$Genus)
)

keep_named <- !(
  is.na(genus0) |
    genus0 == "" |
    genus0 %in% c(
      "NA",
      "Unclassified",
      "unclassified"
    )
)

ps_genus_base <- prune_taxa(
  keep_named,
  ps_genus_base
)

cat("\n========================================\n")
cat("NAMED GENUS INPUT\n")
cat("========================================\n")

cat(
  "ASVs retained with genus assignment:",
  ntaxa(ps_genus_base),
  "\n"
)


# ============================================================
# 3. AGGLOMERATE TO GENUS
# ============================================================

ps_genus_base <- tax_glom(
  ps_genus_base,
  taxrank = "Genus",
  NArm = TRUE
)

cat(
  "Genera after aggregation:",
  ntaxa(ps_genus_base),
  "\n"
)


# ============================================================
# 4. FUNCTION TO RUN ONE ECOLOGICAL CONDITION
# ============================================================

# ============================================================
# 4. FUNCTION TO RUN ONE ECOLOGICAL CONDITION
# ============================================================

run_H3_genus_ancom <- function(ps, condition_name) {
  
  cat("\n\n")
  cat("====================================================\n")
  cat("H3 GENUS ANCOM-BC2 —", condition_name, "\n")
  cat("====================================================\n")
  
  
  # ----------------------------------------------------------
  # A. Extract metadata as a TRUE base data.frame
  # ----------------------------------------------------------
  
  sd_all <- methods::as(
    phyloseq::sample_data(ps),
    "data.frame"
  )
  
  sd_all <- data.frame(
    sd_all,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # B. Safely select requested ecological condition
  # ----------------------------------------------------------
  
  keep_ids <- rownames(sd_all)[
    as.character(sd_all$Ecological_condition) == condition_name
  ]
  
  if (length(keep_ids) == 0) {
    stop(
      "No samples found for condition: ",
      condition_name
    )
  }
  
  ps_sub <- phyloseq::prune_samples(
    keep_ids,
    ps
  )
  
  ps_sub <- phyloseq::prune_taxa(
    phyloseq::taxa_sums(ps_sub) > 0,
    ps_sub
  )
  
  
  # ----------------------------------------------------------
  # C. Clean metadata and set factors
  # ----------------------------------------------------------
  
  sd <- methods::as(
    phyloseq::sample_data(ps_sub),
    "data.frame"
  )
  
  sd <- data.frame(
    sd,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  sd$Site <- factor(
    sd$Site,
    levels = c("ACP", "BASS", "TRUG")
  )
  
  sd$Soil_habitat <- factor(
    sd$Soil_habitat,
    levels = c(
      "Bulk soil",
      "Rhizosphere"
    )
  )
  
  sd$PlotID <- factor(sd$PlotID)
  
  
  # ----------------------------------------------------------
  # D. Verify complete BS-RZ pairing using BASE R
  # ----------------------------------------------------------
  
  pair_tab <- table(
    sd$PlotID,
    sd$Soil_habitat
  )
  
  pair_tab <- pair_tab[
    rowSums(pair_tab) > 0,
    ,
    drop = FALSE
  ]
  
  if (
    !all(
      pair_tab[, "Bulk soil"] == 1 &
      pair_tab[, "Rhizosphere"] == 1
    )
  ) {
    print(pair_tab)
    
    stop(
      "Incomplete pairing detected in ",
      condition_name
    )
  }
  
  cat("Samples:", nrow(sd), "\n")
  cat("Paired plots:", nrow(pair_tab), "\n")
  
  cat("\nSamples by Site × Habitat:\n")
  print(
    table(
      sd$Site,
      sd$Soil_habitat
    )
  )
  
  
  # Put clean metadata back into phyloseq
  phyloseq::sample_data(ps_sub) <-
    phyloseq::sample_data(sd)
  
  
  # ----------------------------------------------------------
  # E. Condition-specific 10% prevalence filter
  # ----------------------------------------------------------
  
  otu <- as(
    phyloseq::otu_table(ps_sub),
    "matrix"
  )
  
  if (!phyloseq::taxa_are_rows(ps_sub)) {
    otu <- t(otu)
  }
  
  prevalence <- rowSums(
    otu > 0
  )
  
  prevalence_threshold <- ceiling(
    0.10 * phyloseq::nsamples(ps_sub)
  )
  
  keep_prev <- prevalence >= prevalence_threshold
  
  cat(
    "\n10% prevalence threshold:",
    prevalence_threshold,
    "samples\n"
  )
  
  cat(
    "Genera before prevalence filtering:",
    phyloseq::ntaxa(ps_sub),
    "\n"
  )
  
  ps_sub <- phyloseq::prune_taxa(
    keep_prev,
    ps_sub
  )
  
  cat(
    "Genera entering ANCOM-BC2:",
    phyloseq::ntaxa(ps_sub),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # F. Run paired ANCOM-BC2
  # ----------------------------------------------------------
  
  set.seed(123)
  
  out <- ANCOMBC::ancombc2(
    data = ps_sub,
    
    tax_level = NULL,
    
    fix_formula =
      "Site + Soil_habitat",
    
    rand_formula =
      "(1 | PlotID)",
    
    p_adj_method =
      "BH",
    
    pseudo_sens =
      TRUE,
    
    prv_cut =
      0,
    
    lib_cut =
      0,
    
    group =
      "Soil_habitat",
    
    struc_zero =
      TRUE,
    
    neg_lb =
      FALSE,
    
    alpha =
      0.05,
    
    global =
      FALSE,
    
    pairwise =
      FALSE,
    
    dunnet =
      FALSE,
    
    trend =
      FALSE,
    
    n_cl =
      1,
    
    verbose =
      TRUE
  )
  
  
  # ----------------------------------------------------------
  # G. Extract ANCOM result
  # ----------------------------------------------------------
  
  res <- as.data.frame(
    out$res,
    stringsAsFactors = FALSE
  )
  
  if (!"taxon" %in% colnames(res)) {
    res$taxon <- rownames(res)
  }
  
  lfc_col <- "lfc_Soil_habitatRhizosphere"
  se_col  <- "se_Soil_habitatRhizosphere"
  W_col   <- "W_Soil_habitatRhizosphere"
  p_col   <- "p_Soil_habitatRhizosphere"
  q_col   <- "q_Soil_habitatRhizosphere"
  robust_col <-
    "diff_robust_Soil_habitatRhizosphere"
  
  required_cols <- c(
    lfc_col,
    se_col,
    W_col,
    p_col,
    q_col,
    robust_col
  )
  
  if (!all(required_cols %in% colnames(res))) {
    stop(
      "Expected Rhizosphere-vs-Bulk coefficient columns were not found."
    )
  }
  
  
  # ----------------------------------------------------------
  # H. Diagnose mixed-model output
  # ----------------------------------------------------------
  
  cat("\n========================================\n")
  cat("ANCOM-BC2 MIXED-MODEL DIAGNOSTIC\n")
  cat("========================================\n")
  
  cat(
    "Total genera in result:",
    nrow(res),
    "\n"
  )
  
  cat(
    "NA LFC:",
    sum(is.na(res[[lfc_col]])),
    "\n"
  )
  
  cat(
    "NA p-values:",
    sum(is.na(res[[p_col]])),
    "\n"
  )
  
  cat(
    "NA q-values:",
    sum(is.na(res[[q_col]])),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # I. Structural-zero table
  # ----------------------------------------------------------
  
  zero_tab <- as.data.frame(
    out$zero_ind,
    stringsAsFactors = FALSE
  )
  
  if (!"taxon" %in% colnames(zero_tab)) {
    zero_tab$taxon <- rownames(zero_tab)
  }
  
  zero_cols <- setdiff(
    colnames(zero_tab),
    "taxon"
  )
  
  if (length(zero_cols) > 0) {
    
    zero_tab$Structural_zero_any <- apply(
      zero_tab[
        ,
        zero_cols,
        drop = FALSE
      ],
      1,
      function(x) {
        any(
          as.logical(x),
          na.rm = TRUE
        )
      }
    )
    
  } else {
    
    zero_tab$Structural_zero_any <- FALSE
  }
  
  
  # ----------------------------------------------------------
  # J. Taxonomy
  # ----------------------------------------------------------
  
  tax_sub <- as.data.frame(
    phyloseq::tax_table(ps_sub),
    stringsAsFactors = FALSE
  )
  
  tax_sub$taxon <- rownames(tax_sub)
  
  
  # ----------------------------------------------------------
  # K. Build clean result table WITHOUT dplyr select()
  # ----------------------------------------------------------
  
  result_clean <- data.frame(
    taxon = res$taxon,
    LFC_RZ_vs_BS = res[[lfc_col]],
    SE = res[[se_col]],
    W = res[[W_col]],
    p = res[[p_col]],
    q_BH = res[[q_col]],
    Robust = as.logical(res[[robust_col]]),
    stringsAsFactors = FALSE
  )
  
  tax_match <- match(
    result_clean$taxon,
    tax_sub$taxon
  )
  
  for (
    rank in intersect(
      c(
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species"
      ),
      colnames(tax_sub)
    )
  ) {
    
    result_clean[[rank]] <-
      tax_sub[[rank]][tax_match]
  }
  
  
  zero_match <- match(
    result_clean$taxon,
    zero_tab$taxon
  )
  
  result_clean$Structural_zero_any <-
    zero_tab$Structural_zero_any[
      zero_match
    ]
  
  result_clean$Structural_zero_any[
    is.na(
      result_clean$Structural_zero_any
    )
  ] <- FALSE
  
  
  result_clean$Direction <- ifelse(
    result_clean$LFC_RZ_vs_BS > 0,
    "Higher in Rhizosphere",
    ifelse(
      result_clean$LFC_RZ_vs_BS < 0,
      "Higher in Bulk soil",
      "No direction"
    )
  )
  
  
  # ----------------------------------------------------------
  # L. Quantitative + robust taxa
  # ----------------------------------------------------------
  
  quantitative <- result_clean[
    !result_clean$Structural_zero_any &
      !is.na(result_clean$LFC_RZ_vs_BS),
    ,
    drop = FALSE
  ]
  
  robust_quant <- quantitative[
    !is.na(quantitative$Robust) &
      quantitative$Robust,
    ,
    drop = FALSE
  ]
  
  robust_quant <- robust_quant[
    order(
      abs(robust_quant$LFC_RZ_vs_BS),
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  n_structural <- sum(
    zero_tab$Structural_zero_any,
    na.rm = TRUE
  )
  
  
  # ----------------------------------------------------------
  # M. Print final summary
  # ----------------------------------------------------------
  
  cat("\n========================================\n")
  cat(condition_name, "— SUMMARY\n")
  cat("========================================\n")
  
  cat(
    "Genera entering analysis:",
    phyloseq::ntaxa(ps_sub),
    "\n"
  )
  
  cat(
    "Genera with >=1 structural zero:",
    n_structural,
    "\n"
  )
  
  cat(
    "Genera quantitatively modelled:",
    nrow(quantitative),
    "\n"
  )
  
  cat(
    "Robust quantitative BS-RZ genera:",
    nrow(robust_quant),
    "\n"
  )
  
  if (nrow(robust_quant) > 0) {
    
    cat("\nROBUST GENERA\n")
    
    print(
      robust_quant[
        ,
        c(
          "Genus",
          "Phylum",
          "LFC_RZ_vs_BS",
          "SE",
          "q_BH",
          "Direction"
        ),
        drop = FALSE
      ],
      row.names = FALSE
    )
    
  } else {
    
    cat(
      "\nNo robust quantitative genera detected.\n"
    )
  }
  
  
  # ----------------------------------------------------------
  # N. Return
  # ----------------------------------------------------------
  
  return(
    list(
      condition = condition_name,
      ps = ps_sub,
      ancombc2 = out,
      zero_table = zero_tab,
      full_results = result_clean,
      quantitative = quantitative,
      robust = robust_quant,
      prevalence_threshold =
        prevalence_threshold
    )
  )
}

# ============================================================
# 9. RUN NATURAL / BESE / BARE
# ============================================================

ancom_C_Natural <- run_H3_genus_ancom(
  ps_genus_base,
  "Natural"
)

ancom_C_BESE <- run_H3_genus_ancom(
  ps_genus_base,
  "BESE"
)

ancom_C_BARE <- run_H3_genus_ancom(
  ps_genus_base,
  "BARE"
)


# ============================================================
# 10. COMBINE ROBUST RESULTS
# ============================================================

robust_all_C <- bind_rows(
  ancom_C_Natural$robust %>%
    mutate(
      Ecological_condition = "Natural"
    ),
  
  ancom_C_BESE$robust %>%
    mutate(
      Ecological_condition = "BESE"
    ),
  
  ancom_C_BARE$robust %>%
    mutate(
      Ecological_condition = "BARE"
    )
) %>%
  relocate(
    Ecological_condition
  )


cat("\n\n========================================\n")
cat("H3 — ROBUST GENUS SUMMARY\n")
cat("========================================\n")

if (nrow(robust_all_C) > 0) {
  
  print(
    robust_all_C %>%
      select(
        Ecological_condition,
        Genus,
        Phylum,
        LFC_RZ_vs_BS,
        q_BH,
        Direction
      ),
    n = Inf
  )
  
} else {
  
  cat(
    "No robust quantitative genera detected.\n"
  )
}


# ============================================================
# 11. SAVE EVERYTHING
# ============================================================

dir.create(
  "outputs/H3_ANCOMBC2/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_ANCOMBC2/models",
  recursive = TRUE,
  showWarnings = FALSE
)


write.csv(
  ancom_C_Natural$full_results,
  "outputs/H3_ANCOMBC2/tables/H3_Natural_genus_all.csv",
  row.names = FALSE
)

write.csv(
  ancom_C_BESE$full_results,
  "outputs/H3_ANCOMBC2/tables/H3_BESE_genus_all.csv",
  row.names = FALSE
)

write.csv(
  ancom_C_BARE$full_results,
  "outputs/H3_ANCOMBC2/tables/H3_BARE_genus_all.csv",
  row.names = FALSE
)

write.csv(
  robust_all_C,
  "outputs/H3_ANCOMBC2/tables/H3_genus_robust_BS_RZ.csv",
  row.names = FALSE
)


saveRDS(
  list(
    Natural = ancom_C_Natural,
    BESE = ancom_C_BESE,
    BARE = ancom_C_BARE,
    robust_all = robust_all_C,
    ANCOMBC_version =
      as.character(
        packageVersion("ANCOMBC")
      ),
    sessionInfo = sessionInfo()
  ),
  "outputs/H3_ANCOMBC2/models/H3_genus_paired_ANCOMBC2.rds"
)


cat("\n========================================\n")
cat("H3 GENUS ANCOM-BC2 COMPLETE\n")
cat("========================================\n")

# ============================================================
# H3 — NATURAL GENUS ANCOM-BC2
# SENSITIVITY CHECK:
# PlotID as FIXED blocking factor instead of random effect
#
# Question:
# Are the robust BS-RZ genera consistent with the mixed model?
# ============================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(phyloseq)
})

set.seed(123)

# ------------------------------------------------------------
# 1. Use the exact same Natural input used in the mixed model
# ------------------------------------------------------------

ps_nat_fixed <- ancom_C_Natural$ps

sd_nat <- methods::as(
  phyloseq::sample_data(ps_nat_fixed),
  "data.frame"
)

sd_nat <- data.frame(
  sd_nat,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sd_nat$PlotID <- factor(sd_nat$PlotID)

sd_nat$Soil_habitat <- factor(
  sd_nat$Soil_habitat,
  levels = c(
    "Bulk soil",
    "Rhizosphere"
  )
)

phyloseq::sample_data(ps_nat_fixed) <-
  phyloseq::sample_data(sd_nat)


# ------------------------------------------------------------
# 2. Run ANCOM-BC2 with PlotID as fixed blocking factor
# ------------------------------------------------------------

set.seed(123)

nat_fixed_out <- ANCOMBC::ancombc2(
  data = ps_nat_fixed,
  
  tax_level = NULL,
  
  fix_formula =
    "PlotID + Soil_habitat",
  
  rand_formula =
    NULL,
  
  p_adj_method =
    "BH",
  
  pseudo_sens =
    TRUE,
  
  prv_cut =
    0,
  
  lib_cut =
    0,
  
  group =
    "Soil_habitat",
  
  struc_zero =
    TRUE,
  
  neg_lb =
    FALSE,
  
  alpha =
    0.05,
  
  global =
    FALSE,
  
  pairwise =
    FALSE,
  
  dunnet =
    FALSE,
  
  trend =
    FALSE,
  
  n_cl =
    1,
  
  verbose =
    TRUE
)


# ------------------------------------------------------------
# 3. Extract rhizosphere coefficient
# ------------------------------------------------------------

res_fixed <- as.data.frame(
  nat_fixed_out$res,
  stringsAsFactors = FALSE
)

lfc_col <- "lfc_Soil_habitatRhizosphere"
se_col  <- "se_Soil_habitatRhizosphere"
p_col   <- "p_Soil_habitatRhizosphere"
q_col   <- "q_Soil_habitatRhizosphere"
rob_col <- "diff_robust_Soil_habitatRhizosphere"


# ------------------------------------------------------------
# 4. Taxonomy lookup
# ------------------------------------------------------------

tax_nat <- as.data.frame(
  phyloseq::tax_table(ps_nat_fixed),
  stringsAsFactors = FALSE
)

tax_nat$taxon <- rownames(tax_nat)

tax_match <- match(
  res_fixed$taxon,
  tax_nat$taxon
)


# ------------------------------------------------------------
# 5. Build sensitivity result table
# ------------------------------------------------------------

fixed_tab <- data.frame(
  taxon = res_fixed$taxon,
  Genus = tax_nat$Genus[tax_match],
  Phylum = tax_nat$Phylum[tax_match],
  
  LFC_fixed =
    res_fixed[[lfc_col]],
  
  SE_fixed =
    res_fixed[[se_col]],
  
  p_fixed =
    res_fixed[[p_col]],
  
  q_fixed =
    res_fixed[[q_col]],
  
  Robust_fixed =
    as.logical(
      res_fixed[[rob_col]]
    ),
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 6. Mixed-model robust genera
# ------------------------------------------------------------

mixed_tab <- ancom_C_Natural$robust[
  ,
  c(
    "taxon",
    "Genus",
    "Phylum",
    "LFC_RZ_vs_BS",
    "q_BH"
  ),
  drop = FALSE
]

names(mixed_tab)[
  names(mixed_tab) == "LFC_RZ_vs_BS"
] <- "LFC_mixed"

names(mixed_tab)[
  names(mixed_tab) == "q_BH"
] <- "q_mixed"


# ------------------------------------------------------------
# 7. Compare mixed vs fixed blocking
# ------------------------------------------------------------

comparison <- merge(
  mixed_tab,
  fixed_tab[
    ,
    c(
      "taxon",
      "LFC_fixed",
      "q_fixed",
      "Robust_fixed"
    )
  ],
  by = "taxon",
  all.x = TRUE
)

comparison$Same_direction <-
  sign(comparison$LFC_mixed) ==
  sign(comparison$LFC_fixed)


cat("\n========================================\n")
cat("NATURAL — FIXED-BLOCK SENSITIVITY\n")
cat("========================================\n")

cat(
  "Genera in fixed model:",
  nrow(res_fixed),
  "\n"
)

cat(
  "NA LFC in fixed model:",
  sum(is.na(res_fixed[[lfc_col]])),
  "\n"
)

cat(
  "Robust genera in fixed model:",
  sum(
    res_fixed[[rob_col]] %in% TRUE,
    na.rm = TRUE
  ),
  "\n"
)

cat("\nMIXED-MODEL ROBUST GENERA — SENSITIVITY CHECK\n")

print(
  comparison[
    ,
    c(
      "Genus",
      "LFC_mixed",
      "LFC_fixed",
      "q_mixed",
      "q_fixed",
      "Robust_fixed",
      "Same_direction"
    )
  ],
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Overall overlap
# ------------------------------------------------------------

mixed_robust_taxa <-
  ancom_C_Natural$robust$taxon

fixed_robust_taxa <-
  res_fixed$taxon[
    res_fixed[[rob_col]] %in% TRUE
  ]

cat("\n========================================\n")
cat("ROBUST-TAXON OVERLAP\n")
cat("========================================\n")

cat(
  "Mixed robust:",
  length(mixed_robust_taxa),
  "\n"
)

cat(
  "Fixed-block robust:",
  length(fixed_robust_taxa),
  "\n"
)

cat(
  "Shared robust:",
  length(
    intersect(
      mixed_robust_taxa,
      fixed_robust_taxa
    )
  ),
  "\n"
)

cat(
  "Mixed robust with same LFC direction:",
  sum(
    comparison$Same_direction,
    na.rm = TRUE
  ),
  "of",
  nrow(comparison),
  "\n"
)

# ============================================================
# H3 — GENUS ANCOM-BC2
# RUN BESE + BARE AND COMBINE WITH NATURAL
# ============================================================

# ------------------------------------------------------------
# 1. Run BESE
# ------------------------------------------------------------

ancom_C_BESE <- run_H3_genus_ancom(
  ps_genus_base,
  "BESE"
)


# ------------------------------------------------------------
# 2. Run BARE
# ------------------------------------------------------------

ancom_C_BARE <- run_H3_genus_ancom(
  ps_genus_base,
  "BARE"
)


# ------------------------------------------------------------
# 3. Combine robust genera from all three conditions
# ------------------------------------------------------------

nat_tab <- ancom_C_Natural$robust
nat_tab$Ecological_condition <- "Natural"

bese_tab <- ancom_C_BESE$robust
bese_tab$Ecological_condition <- "BESE"

bare_tab <- ancom_C_BARE$robust
bare_tab$Ecological_condition <- "BARE"

robust_all_C <- rbind(
  nat_tab,
  bese_tab,
  bare_tab
)


# ------------------------------------------------------------
# 4. Arrange useful columns
# ------------------------------------------------------------

keep_cols <- c(
  "Ecological_condition",
  "Genus",
  "Phylum",
  "LFC_RZ_vs_BS",
  "SE",
  "p",
  "q_BH",
  "Direction"
)

robust_summary_C <- robust_all_C[
  ,
  keep_cols,
  drop = FALSE
]

robust_summary_C <- robust_summary_C[
  order(
    robust_summary_C$Ecological_condition,
    -abs(robust_summary_C$LFC_RZ_vs_BS)
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 5. Print summary counts
# ------------------------------------------------------------

cat("\n========================================\n")
cat("H3 — GENUS ANCOM-BC2 SUMMARY\n")
cat("========================================\n")

cat(
  "Natural robust genera:",
  nrow(ancom_C_Natural$robust),
  "\n"
)

cat(
  "BESE robust genera:",
  nrow(ancom_C_BESE$robust),
  "\n"
)

cat(
  "BARE robust genera:",
  nrow(ancom_C_BARE$robust),
  "\n"
)


# ------------------------------------------------------------
# 6. Print all robust genera
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ROBUST BS-RZ GENERA\n")
cat("========================================\n")

print(
  robust_summary_C,
  row.names = FALSE
)


# ------------------------------------------------------------
# 7. Count direction within each condition
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DIRECTION SUMMARY\n")
cat("========================================\n")

print(
  with(
    robust_summary_C,
    table(
      Ecological_condition,
      Direction
    )
  )
)


# ------------------------------------------------------------
# 8. Identify genera shared among conditions
# ------------------------------------------------------------

nat_genera <- unique(
  ancom_C_Natural$robust$Genus
)

bese_genera <- unique(
  ancom_C_BESE$robust$Genus
)

bare_genera <- unique(
  ancom_C_BARE$robust$Genus
)

shared_all <- Reduce(
  intersect,
  list(
    nat_genera,
    bese_genera,
    bare_genera
  )
)

shared_nat_bese <- intersect(
  nat_genera,
  bese_genera
)

shared_nat_bare <- intersect(
  nat_genera,
  bare_genera
)

shared_bese_bare <- intersect(
  bese_genera,
  bare_genera
)

cat("\n========================================\n")
cat("ROBUST GENUS OVERLAP\n")
cat("========================================\n")

cat(
  "Shared by all three:",
  length(shared_all),
  "\n"
)

if (length(shared_all) > 0) {
  print(shared_all)
}

cat(
  "\nNatural + BESE:",
  length(shared_nat_bese),
  "\n"
)

if (length(shared_nat_bese) > 0) {
  print(shared_nat_bese)
}

cat(
  "\nNatural + BARE:",
  length(shared_nat_bare),
  "\n"
)

if (length(shared_nat_bare) > 0) {
  print(shared_nat_bare)
}

cat(
  "\nBESE + BARE:",
  length(shared_bese_bare),
  "\n"
)

if (length(shared_bese_bare) > 0) {
  print(shared_bese_bare)
}


# ------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------

dir.create(
  "outputs/H3_ANCOMBC2/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_ANCOMBC2/models",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  robust_summary_C,
  "outputs/H3_ANCOMBC2/tables/H3_genus_robust_BS_RZ.csv",
  row.names = FALSE
)

write.csv(
  ancom_C_Natural$full_results,
  "outputs/H3_ANCOMBC2/tables/H3_Natural_genus_all.csv",
  row.names = FALSE
)

write.csv(
  ancom_C_BESE$full_results,
  "outputs/H3_ANCOMBC2/tables/H3_BESE_genus_all.csv",
  row.names = FALSE
)

write.csv(
  ancom_C_BARE$full_results,
  "outputs/H3_ANCOMBC2/tables/H3_BARE_genus_all.csv",
  row.names = FALSE
)

saveRDS(
  list(
    Natural = ancom_C_Natural,
    BESE = ancom_C_BESE,
    BARE = ancom_C_BARE,
    robust_summary = robust_summary_C,
    shared_all = shared_all,
    shared_Natural_BESE = shared_nat_bese,
    shared_Natural_BARE = shared_nat_bare,
    shared_BESE_BARE = shared_bese_bare,
    ANCOMBC_version =
      as.character(
        packageVersion("ANCOMBC")
      ),
    sessionInfo = sessionInfo()
  ),
  "outputs/H3_ANCOMBC2/models/H3_genus_paired_ANCOMBC2.rds"
)

cat("\n========================================\n")
cat("H3 GENUS ANCOM-BC2 COMPLETE\n")
cat("========================================\n")

# ============================================================
# H3 — FUNGI BETA DIVERSITY
# STEP 1 — VERIFY COMPLETE BS–RZ PAIRS
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
})

stopifnot(
  exists("ps_fun_C"),
  exists("alpha_fungi_c")
)

# Metadata already restricted to Natural, BESE, BARE
# and complete BS–RZ pairs in the H3 alpha workflow
pair_meta_fun <- alpha_fungi_c %>%
  dplyr::select(
    SampleID,
    Site,
    Ecological_condition,
    Soil_habitat,
    PlotID
  ) %>%
  dplyr::distinct()

# Check that every retained plot contains exactly one BS and one RZ
pair_check_fun <- pair_meta_fun %>%
  dplyr::group_by(PlotID) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_habitats = dplyr::n_distinct(Soil_habitat),
    .groups = "drop"
  )

stopifnot(
  all(pair_check_fun$n_samples == 2),
  all(pair_check_fun$n_habitats == 2)
)

# Check all samples exist in the corrected fungal phyloseq object
missing_fun_ids <- setdiff(
  pair_meta_fun$SampleID,
  sample_names(ps_fun_C)
)

cat("\n========================================\n")
cat("H3 FUNGI — BETA PAIR CHECK\n")
cat("========================================\n")

cat(
  "Samples:",
  nrow(pair_meta_fun),
  "\n"
)

cat(
  "Paired plots:",
  dplyr::n_distinct(pair_meta_fun$PlotID),
  "\n"
)

cat(
  "Missing sample IDs:",
  length(missing_fun_ids),
  "\n\n"
)

cat("Pairs by ecological condition:\n")

print(
  pair_meta_fun %>%
    dplyr::distinct(
      PlotID,
      Site,
      Ecological_condition
    ) %>%
    dplyr::count(
      Ecological_condition,
      name = "Paired_plots"
    )
)

cat("\nPairs by Site × ecological condition:\n")

print(
  pair_meta_fun %>%
    dplyr::distinct(
      PlotID,
      Site,
      Ecological_condition
    ) %>%
    dplyr::count(
      Site,
      Ecological_condition,
      name = "Paired_plots"
    )
)

cat("\nSamples by ecological condition × habitat:\n")

print(
  pair_meta_fun %>%
    dplyr::count(
      Ecological_condition,
      Soil_habitat,
      name = "Samples"
    )
)


# ============================================================
# H3 — FUNGI BETA DIVERSITY
# STEP 2 — BRAY-CURTIS + PCoA FOR VERIFIED PAIRED SAMPLES
#
# Complete pairs:
# Natural = 15
# BESE    = 10
# BARE    = 9
# Total   = 34 paired plots / 68 samples
#
# Bray-Curtis on relative abundance
# Lingoes correction for PCoA
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
})

stopifnot(
  exists("ps_fun_C"),
  exists("pair_meta_fun")
)

# ------------------------------------------------------------
# 1. KEEP ONLY VERIFIED COMPLETE-PAIR SAMPLES
# ------------------------------------------------------------

ps_C_beta_fun <- phyloseq::prune_samples(
  pair_meta_fun$SampleID,
  ps_fun_C
)

# Remove taxa absent from this paired subset
ps_C_beta_fun <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_C_beta_fun) > 0,
  ps_C_beta_fun
)

cat("\n========================================\n")
cat("H3 FUNGI — PAIRED BETA DATA\n")
cat("========================================\n")

cat(
  "Samples retained:",
  phyloseq::nsamples(ps_C_beta_fun),
  "\n"
)

cat(
  "Taxa retained:",
  phyloseq::ntaxa(ps_C_beta_fun),
  "\n"
)

# ------------------------------------------------------------
# 2. RELATIVE ABUNDANCE
# ------------------------------------------------------------

ps_C_rel_fun <- phyloseq::transform_sample_counts(
  ps_C_beta_fun,
  function(x) x / sum(x)
)

# Safety check
sample_sums_fun <- phyloseq::sample_sums(ps_C_rel_fun)

cat(
  "Relative-abundance sample sums range:",
  min(sample_sums_fun),
  "to",
  max(sample_sums_fun),
  "\n"
)

# ------------------------------------------------------------
# 3. BRAY-CURTIS DISTANCE
# ------------------------------------------------------------

bray_C_fun <- phyloseq::distance(
  ps_C_rel_fun,
  method = "bray"
)

cat(
  "Bray-Curtis samples:",
  length(attr(bray_C_fun, "Labels")),
  "\n"
)

# ------------------------------------------------------------
# 4. PCoA WITH LINGOES CORRECTION
# ------------------------------------------------------------

pcoa_C_fun <- vegan::wcmdscale(
  bray_C_fun,
  eig = TRUE,
  add = "lingoes"
)

positive_eig_fun <- pcoa_C_fun$eig[
  pcoa_C_fun$eig > 0
]

axis1_fun <- round(
  100 * positive_eig_fun[1] /
    sum(positive_eig_fun),
  1
)

axis2_fun <- round(
  100 * positive_eig_fun[2] /
    sum(positive_eig_fun),
  1
)

cat("\n========================================\n")
cat("FUNGI PCoA\n")
cat("========================================\n")

cat(
  "PCoA1:",
  axis1_fun,
  "%\n"
)

cat(
  "PCoA2:",
  axis2_fun,
  "%\n"
)

cat(
  "PCoA1 + PCoA2:",
  round(axis1_fun + axis2_fun, 1),
  "%\n"
)

# ============================================================
# FIX — ALIGN FUNGAL METADATA TO BRAY-CURTIS MATRIX
# ============================================================

sample_order_fun <- attr(
  bray_C_fun,
  "Labels"
)

meta_C_fun <- pair_meta_fun %>%
  dplyr::filter(
    SampleID %in% sample_order_fun
  ) %>%
  as.data.frame()

# Reorder metadata to exactly match Bray-Curtis sample order
meta_C_fun <- meta_C_fun[
  match(
    sample_order_fun,
    meta_C_fun$SampleID
  ),
  ,
  drop = FALSE
]

# Set row names only after conversion to data.frame
rownames(meta_C_fun) <-
  meta_C_fun$SampleID

# Factors
meta_C_fun$Site <- factor(
  meta_C_fun$Site,
  levels = c(
    "ACP",
    "BASS",
    "TRUG"
  )
)

meta_C_fun$Ecological_condition <- factor(
  meta_C_fun$Ecological_condition,
  levels = c(
    "Natural",
    "BESE",
    "BARE"
  )
)

meta_C_fun$Soil_habitat <- factor(
  meta_C_fun$Soil_habitat,
  levels = c(
    "Bulk soil",
    "Rhizosphere"
  )
)

meta_C_fun$PlotID <- factor(
  meta_C_fun$PlotID
)

# ============================================================
# SAFETY CHECKS
# ============================================================

cat("\n========================================\n")
cat("FUNGI METADATA ALIGNMENT CHECK\n")
cat("========================================\n")

cat(
  "Metadata rows:",
  nrow(meta_C_fun),
  "\n"
)

cat(
  "Bray-Curtis labels:",
  length(sample_order_fun),
  "\n"
)

cat(
  "Paired plots:",
  dplyr::n_distinct(meta_C_fun$PlotID),
  "\n"
)

cat(
  "Row names identical to Bray order:",
  identical(
    rownames(meta_C_fun),
    sample_order_fun
  ),
  "\n\n"
)

print(
  meta_C_fun %>%
    dplyr::distinct(
      PlotID,
      Site,
      Ecological_condition
    ) %>%
    dplyr::count(
      Site,
      Ecological_condition,
      name = "Paired_plots"
    )
)

stopifnot(
  nrow(meta_C_fun) == 68,
  length(sample_order_fun) == 68,
  dplyr::n_distinct(
    meta_C_fun$PlotID
  ) == 34,
  identical(
    rownames(meta_C_fun),
    sample_order_fun
  )
)

# ============================================================
# H3 — FUNGI BETA DIVERSITY
# STEP 3 — PAIRED PERMANOVA
#
# H3:
# Does fungal BS–RZ compositional differentiation
# vary among Natural, BESE and BARE?
#
# Bray-Curtis on relative abundance
# 9999 permutations restricted within PlotID
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
  library(permute)
})

set.seed(123)

stopifnot(
  exists("bray_C_fun"),
  exists("meta_C_fun")
)

# ============================================================
# 1. RESTRICT PERMUTATIONS WITHIN PLOT
# ============================================================

perm_control_fun <- permute::how(
  blocks = meta_C_fun$PlotID,
  nperm = 9999
)

perm_matrix_fun <- permute::shuffleSet(
  n = nrow(meta_C_fun),
  nset = 9999,
  control = perm_control_fun
)

cat("\n========================================\n")
cat("FUNGAL H3 — PAIRED PERMANOVA\n")
cat("========================================\n")

cat(
  "Samples:",
  nrow(meta_C_fun),
  "\n"
)

cat(
  "Paired plots:",
  dplyr::n_distinct(meta_C_fun$PlotID),
  "\n"
)

# ============================================================
# 2. OVERALL BS vs RZ EFFECT
#
# Site and ecological condition are included first.
# With permutations restricted within PlotID,
# inference focuses on Soil_habitat.
# ============================================================

perm_overall_fun <- vegan::adonis2(
  bray_C_fun ~
    Site +
    Ecological_condition +
    Soil_habitat,
  data = meta_C_fun,
  permutations = perm_matrix_fun,
  by = "terms"
)

cat("\n========================================\n")
cat("1. OVERALL FUNGAL BS vs RZ EFFECT\n")
cat("========================================\n")

print(perm_overall_fun)

# ============================================================
# 3. ECOLOGICAL CONDITION × SOIL HABITAT
#
# Primary H3 test:
# Does the BS–RZ compositional shift differ among
# Natural, BESE and BARE?
# ============================================================

perm_interaction_fun <- vegan::adonis2(
  bray_C_fun ~
    Site +
    Ecological_condition +
    Soil_habitat +
    Ecological_condition:Soil_habitat,
  data = meta_C_fun,
  permutations = perm_matrix_fun,
  by = "terms"
)

cat("\n========================================\n")
cat("2. CONDITION × SOIL-HABITAT INTERACTION\n")
cat("========================================\n")

print(perm_interaction_fun)

# ============================================================
# 4. BS vs RZ WITHIN EACH ECOLOGICAL CONDITION
# ============================================================

bray_matrix_fun <- as.matrix(
  bray_C_fun
)

conditions_fun <- c(
  "Natural",
  "BESE",
  "BARE"
)

within_results_fun <- list()

for (cond in conditions_fun) {
  
  meta_sub_fun <- meta_C_fun %>%
    dplyr::filter(
      Ecological_condition == cond
    ) %>%
    droplevels()
  
  ids_fun <- rownames(
    meta_sub_fun
  )
  
  bray_sub_fun <- as.dist(
    bray_matrix_fun[
      ids_fun,
      ids_fun
    ]
  )
  
  # paired permutations within PlotID
  perm_sub_control_fun <- permute::how(
    blocks = meta_sub_fun$PlotID,
    nperm = 9999
  )
  
  perm_sub_fun <- permute::shuffleSet(
    n = nrow(meta_sub_fun),
    nset = 9999,
    control = perm_sub_control_fun
  )
  
  fit_sub_fun <- vegan::adonis2(
    bray_sub_fun ~
      Site +
      Soil_habitat,
    data = meta_sub_fun,
    permutations = perm_sub_fun,
    by = "terms"
  )
  
  within_results_fun[[cond]] <-
    fit_sub_fun
  
  cat("\n----------------------------------------\n")
  cat(
    cond,
    ": FUNGAL BULK SOIL vs RHIZOSPHERE\n"
  )
  cat("----------------------------------------\n")
  
  print(
    fit_sub_fun
  )
}

# ============================================================
# 5. EXTRACT THE THREE WITHIN-CONDITION TESTS
# ============================================================

within_table_fun <- dplyr::bind_rows(
  lapply(
    conditions_fun,
    function(cond) {
      
      tab_fun <- as.data.frame(
        within_results_fun[[cond]]
      )
      
      habitat_row_fun <- tab_fun[
        rownames(tab_fun) == "Soil_habitat",
        ,
        drop = FALSE
      ]
      
      data.frame(
        Ecological_condition = cond,
        F = habitat_row_fun$F,
        R2 = habitat_row_fun$R2,
        p_raw =
          habitat_row_fun$`Pr(>F)`
      )
    }
  )
)

# Holm correction across the three planned tests
within_table_fun$p_Holm <- p.adjust(
  within_table_fun$p_raw,
  method = "holm"
)

within_table_fun$R2_percent <-
  100 * within_table_fun$R2

cat("\n========================================\n")
cat("3. WITHIN-CONDITION FUNGAL SUMMARY\n")
cat("========================================\n")

print(
  within_table_fun,
  row.names = FALSE
)

# ============================================================
# 6. EXTRACT PRIMARY H3 RESULTS
# ============================================================

overall_df_fun <- as.data.frame(
  perm_overall_fun
)

interaction_df_fun <- as.data.frame(
  perm_interaction_fun
)

primary_results_fun <- data.frame(
  Test = c(
    "Overall BS vs RZ",
    "Ecological condition x Soil habitat"
  ),
  
  F = c(
    overall_df_fun[
      "Soil_habitat",
      "F"
    ],
    interaction_df_fun[
      "Ecological_condition:Soil_habitat",
      "F"
    ]
  ),
  
  R2 = c(
    overall_df_fun[
      "Soil_habitat",
      "R2"
    ],
    interaction_df_fun[
      "Ecological_condition:Soil_habitat",
      "R2"
    ]
  ),
  
  p = c(
    overall_df_fun[
      "Soil_habitat",
      "Pr(>F)"
    ],
    interaction_df_fun[
      "Ecological_condition:Soil_habitat",
      "Pr(>F)"
    ]
  )
)

primary_results_fun$R2_percent <-
  100 * primary_results_fun$R2

cat("\n========================================\n")
cat("PRIMARY FUNGAL H3 RESULTS\n")
cat("========================================\n")

print(
  primary_results_fun,
  row.names = FALSE
)


# ============================================================
# H3 — FUNGI BETA DIVERSITY
# STEP 4 — PERMDISP DIAGNOSTICS
#
# Questions:
# 1. Do the six Condition × Habitat groups differ in dispersion?
# 2. Within Natural, BESE and BARE, do BS and RZ differ
#    in multivariate dispersion?
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(dplyr)
})

set.seed(123)

stopifnot(
  exists("bray_C_fun"),
  exists("meta_C_fun")
)

# ------------------------------------------------------------
# 1. Confirm metadata follows Bray-Curtis order
# ------------------------------------------------------------

sample_order_fun <- attr(
  bray_C_fun,
  "Labels"
)

stopifnot(
  identical(
    rownames(meta_C_fun),
    sample_order_fun
  )
)

# ============================================================
# 2. SIX-GROUP PERMDISP
# Condition × Soil habitat
# ============================================================

meta_C_fun$Condition_habitat <- interaction(
  meta_C_fun$Ecological_condition,
  meta_C_fun$Soil_habitat,
  sep = " | ",
  drop = TRUE
)

bd_all_fun <- vegan::betadisper(
  bray_C_fun,
  group = meta_C_fun$Condition_habitat,
  type = "median",
  bias.adjust = TRUE
)

perm_bd_all_fun <- vegan::permutest(
  bd_all_fun,
  permutations = 9999
)

cat("\n========================================\n")
cat("FUNGAL H3 — GLOBAL PERMDISP\n")
cat("Condition × Soil habitat groups\n")
cat("========================================\n")

print(
  perm_bd_all_fun
)

cat("\nMean distance to group centroid:\n")

dispersion_means_fun <- data.frame(
  Condition_habitat = names(
    tapply(
      bd_all_fun$distances,
      bd_all_fun$group,
      mean
    )
  ),
  Mean_distance = as.numeric(
    tapply(
      bd_all_fun$distances,
      bd_all_fun$group,
      mean
    )
  )
)

print(
  dispersion_means_fun,
  row.names = FALSE
)

# ============================================================
# 3. WITHIN-CONDITION PERMDISP
# BS vs RZ separately in Natural / BESE / BARE
# ============================================================

conditions_fun <- c(
  "Natural",
  "BESE",
  "BARE"
)

bray_matrix_fun <- as.matrix(
  bray_C_fun
)

permdisp_results_fun <- list()
permdisp_summary_fun <- list()

for (cond in conditions_fun) {
  
  meta_sub_fun <- meta_C_fun %>%
    dplyr::filter(
      Ecological_condition == cond
    ) %>%
    droplevels()
  
  ids_fun <- rownames(
    meta_sub_fun
  )
  
  bray_sub_fun <- as.dist(
    bray_matrix_fun[
      ids_fun,
      ids_fun
    ]
  )
  
  bd_sub_fun <- vegan::betadisper(
    bray_sub_fun,
    group = meta_sub_fun$Soil_habitat,
    type = "median",
    bias.adjust = TRUE
  )
  
  perm_sub_fun <- vegan::permutest(
    bd_sub_fun,
    permutations = 9999
  )
  
  permdisp_results_fun[[cond]] <- list(
    betadisper = bd_sub_fun,
    permutest = perm_sub_fun
  )
  
  perm_tab_fun <- as.data.frame(
    perm_sub_fun$tab
  )
  
  p_value_fun <- perm_tab_fun[
    1,
    "Pr(>F)"
  ]
  
  F_value_fun <- perm_tab_fun[
    1,
    "F"
  ]
  
  mean_disp_fun <- tapply(
    bd_sub_fun$distances,
    bd_sub_fun$group,
    mean
  )
  
  permdisp_summary_fun[[cond]] <- data.frame(
    Ecological_condition = cond,
    F = F_value_fun,
    p_raw = p_value_fun,
    BS_mean_dispersion =
      unname(
        mean_disp_fun["Bulk soil"]
      ),
    RZ_mean_dispersion =
      unname(
        mean_disp_fun["Rhizosphere"]
      )
  )
  
  cat("\n----------------------------------------\n")
  cat(
    cond,
    ": FUNGAL BS vs RZ PERMDISP\n"
  )
  cat("----------------------------------------\n")
  
  print(
    perm_sub_fun
  )
  
  cat("\nMean distance to centroid:\n")
  
  print(
    mean_disp_fun
  )
}

# ============================================================
# 4. COMBINE WITHIN-CONDITION RESULTS
# ============================================================

permdisp_table_fun <- dplyr::bind_rows(
  permdisp_summary_fun
)

permdisp_table_fun$p_Holm <- p.adjust(
  permdisp_table_fun$p_raw,
  method = "holm"
)

permdisp_table_fun <-
  permdisp_table_fun %>%
  dplyr::mutate(
    Difference_RZ_minus_BS =
      RZ_mean_dispersion -
      BS_mean_dispersion
  )

cat("\n========================================\n")
cat("FUNGAL H3 — WITHIN-CONDITION PERMDISP SUMMARY\n")
cat("========================================\n")

print(
  permdisp_table_fun,
  row.names = FALSE
)



# ============================================================
# H3 — FUNGI BETA DIVERSITY
# STEP 5 — MAGNITUDE OF WITHIN-PLOT BS–RZ DIFFERENTIATION
#
# Response:
# Bray-Curtis distance between paired BS and RZ samples
#
# H3 question:
# Does the magnitude of BS–RZ differentiation differ among
# Natural, BESE and BARE after accounting for Site?
#
# n = 34 paired plots
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(emmeans)
})

set.seed(123)

stopifnot(
  exists("bray_C_fun"),
  exists("meta_C_fun")
)

bray_matrix_fun <- as.matrix(
  bray_C_fun
)

# ============================================================
# 1. CALCULATE ONE BS–RZ BRAY DISTANCE PER PLOT
# ============================================================

plot_ids_fun <- unique(
  meta_C_fun$PlotID
)

paired_bray_fun_list <- lapply(
  plot_ids_fun,
  function(pid) {
    
    tmp <- meta_C_fun[
      meta_C_fun$PlotID == pid,
      ,
      drop = FALSE
    ]
    
    # Exactly one BS and one RZ required
    if (
      nrow(tmp) != 2 ||
      !all(
        c(
          "Bulk soil",
          "Rhizosphere"
        ) %in% tmp$Soil_habitat
      )
    ) {
      stop(
        "Plot ",
        pid,
        " does not contain exactly one BS and one RZ sample."
      )
    }
    
    bs_id <- rownames(tmp)[
      tmp$Soil_habitat == "Bulk soil"
    ]
    
    rz_id <- rownames(tmp)[
      tmp$Soil_habitat == "Rhizosphere"
    ]
    
    data.frame(
      PlotID = as.character(pid),
      Site = as.character(
        tmp$Site[1]
      ),
      Ecological_condition =
        as.character(
          tmp$Ecological_condition[1]
        ),
      Bray_BS_RZ =
        bray_matrix_fun[
          bs_id,
          rz_id
        ]
    )
  }
)

paired_bray_fun <- dplyr::bind_rows(
  paired_bray_fun_list
)

# ============================================================
# 2. FACTORS
# ============================================================

paired_bray_fun <- paired_bray_fun %>%
  dplyr::mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Ecological_condition = factor(
      Ecological_condition,
      levels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    )
  )

# ============================================================
# 3. VERIFY DESIGN + RAW DESCRIPTIVES
# ============================================================

cat("\n========================================\n")
cat("FUNGAL PAIRED BS–RZ BRAY DISTANCE\n")
cat("========================================\n")

cat(
  "Total paired plots:",
  nrow(paired_bray_fun),
  "\n\n"
)

print(
  paired_bray_fun %>%
    dplyr::count(
      Site,
      Ecological_condition,
      name = "Paired_plots"
    )
)

cat("\nRaw Bray–Curtis distance by condition:\n")

print(
  paired_bray_fun %>%
    dplyr::group_by(
      Ecological_condition
    ) %>%
    dplyr::summarise(
      n = dplyr::n(),
      Mean = mean(Bray_BS_RZ),
      SD = sd(Bray_BS_RZ),
      Median = median(Bray_BS_RZ),
      .groups = "drop"
    )
)

stopifnot(
  nrow(paired_bray_fun) == 34
)

# ============================================================
# 4. SITE-ADJUSTED MODEL
# ============================================================

mod_bray_magnitude_fun <- lm(
  Bray_BS_RZ ~
    Site +
    Ecological_condition,
  data = paired_bray_fun
)

cat("\n========================================\n")
cat("FUNGAL BRAY MAGNITUDE — ADJUSTED TERM TESTS\n")
cat("========================================\n")

anova_bray_fun <- drop1(
  mod_bray_magnitude_fun,
  test = "F"
)

print(
  anova_bray_fun
)

# ============================================================
# 5. RESIDUAL DIAGNOSTIC
# ============================================================

shapiro_bray_fun <- shapiro.test(
  residuals(
    mod_bray_magnitude_fun
  )
)

cat("\n========================================\n")
cat("RESIDUAL NORMALITY\n")
cat("========================================\n")

print(
  shapiro_bray_fun
)

# ============================================================
# 6. SITE-RESTRICTED RESIDUAL PERMUTATION TEST
#
# Reduced model = Site
# Test ecological condition after accounting for Site
# ============================================================

mod_reduced_fun <- lm(
  Bray_BS_RZ ~ Site,
  data = paired_bray_fun
)

obs_F_fun <- drop1(
  mod_bray_magnitude_fun,
  test = "F"
)[
  "Ecological_condition",
  "F value"
]

reduced_fitted_fun <- fitted(
  mod_reduced_fun
)

reduced_resid_fun <- residuals(
  mod_reduced_fun
)

nperm_fun <- 9999

perm_F_fun <- numeric(
  nperm_fun
)

for (i in seq_len(nperm_fun)) {
  
  perm_resid_fun <-
    reduced_resid_fun
  
  for (
    s in levels(
      paired_bray_fun$Site
    )
  ) {
    
    idx <- which(
      paired_bray_fun$Site == s
    )
    
    perm_resid_fun[idx] <- sample(
      reduced_resid_fun[idx],
      length(idx),
      replace = FALSE
    )
  }
  
  y_perm_fun <-
    reduced_fitted_fun +
    perm_resid_fun
  
  perm_dat_fun <-
    paired_bray_fun
  
  perm_dat_fun$Bray_perm <-
    y_perm_fun
  
  fit_perm_fun <- lm(
    Bray_perm ~
      Site +
      Ecological_condition,
    data = perm_dat_fun
  )
  
  perm_F_fun[i] <- drop1(
    fit_perm_fun,
    test = "F"
  )[
    "Ecological_condition",
    "F value"
  ]
}

p_perm_fun <- (
  1 +
    sum(
      perm_F_fun >= obs_F_fun
    )
) / (
  nperm_fun + 1
)

cat("\n========================================\n")
cat("SITE-RESTRICTED PERMUTATION TEST\n")
cat("========================================\n")

cat(
  "Observed F =",
  round(
    obs_F_fun,
    4
  ),
  "\n"
)

cat(
  "Permutation p =",
  format(
    p_perm_fun,
    digits = 4
  ),
  "\n"
)

# ============================================================
# 7. SITE-ADJUSTED ESTIMATED MEANS
# ============================================================

emm_bray_fun <- emmeans::emmeans(
  mod_bray_magnitude_fun,
  ~ Ecological_condition
)

cat("\n========================================\n")
cat("SITE-ADJUSTED ESTIMATED MEANS\n")
cat("========================================\n")

print(
  summary(
    emm_bray_fun,
    infer = c(
      TRUE,
      FALSE
    )
  )
)

# ============================================================
# 8. PAIRWISE CONDITION COMPARISONS
# Holm adjustment
# ============================================================

cat("\n========================================\n")
cat("PAIRWISE CONDITION COMPARISONS\n")
cat("Holm adjusted\n")
cat("========================================\n")

print(
  summary(
    pairs(
      emm_bray_fun,
      adjust = "holm"
    ),
    infer = c(
      TRUE,
      TRUE
    )
  )
)

# ============================================================
# 9. PLANNED CONTRAST
# Natural vs mean(BESE + BARE)
# ============================================================

planned_bray_fun <- emmeans::contrast(
  emm_bray_fun,
  method = list(
    "Natural vs mean(BESE + BARE)" =
      c(
        Natural = 1,
        BESE = -0.5,
        BARE = -0.5
      )
  )
)

cat("\n========================================\n")
cat("PLANNED CONTRAST\n")
cat("Natural vs mean(BESE + BARE)\n")
cat("========================================\n")

print(
  summary(
    planned_bray_fun,
    infer = c(
      TRUE,
      TRUE
    )
  )
)


