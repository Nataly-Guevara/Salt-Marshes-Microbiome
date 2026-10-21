# ============================================================
# H2 — RHIZOSPHERE REASSEMBLY
# ALPHA DIVERSITY
#
# Bacteria + Fungi
#
# Strategy:
# 1. Chao1 + Shannon from verified historical phyloseq objects
# 2. Test Site × Condition interaction first
# 3. If interaction significant:
#       retain full model and examine condition within Site
# 4. If interaction not significant:
#       use additive Site + Condition model
#       and examine overall condition differences
# 5. Holm-adjust pairwise comparisons
#
# No rarefaction; same analysis-ready dataset used for
# the established Chapter 1 alpha-diversity workflow.
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(car)
  library(emmeans)
})

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

# ============================================================
# 1. LOAD VERIFIED OBJECTS
# ============================================================

ps_bacteria <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)


# ============================================================
# 2. H2 SUBSETS
# ============================================================

ps_B_bac <- subset_samples(
  ps_bacteria,
  SoilPortion == "RZ" &
    Treatment %in% c("NATURAL", "BESE", "BARE")
)

ps_B_bac <- prune_taxa(
  taxa_sums(ps_B_bac) > 0,
  ps_B_bac
)


ps_B_fun <- subset_samples(
  ps_fungi,
  SoilPortion == "RZ" &
    Treatment %in% c("NATURAL", "BESE", "BARE")
)

ps_B_fun <- prune_taxa(
  taxa_sums(ps_B_fun) > 0,
  ps_B_fun
)


# ============================================================
# 3. FUNCTION TO BUILD ALPHA-DIVERSITY DATA
# ============================================================

make_alpha_data <- function(ps) {
  
  rich <- estimate_richness(
    ps,
    measures = c("Chao1", "Shannon")
  )
  
  # phyloseq sometimes adds X before numeric sample IDs
  rownames(rich) <- sub(
    "^X(?=\\d)",
    "",
    rownames(rich),
    perl = TRUE
  )
  
  meta <- as(
    sample_data(ps),
    "data.frame"
  )
  
  # Match explicitly by sample ID
  rich <- rich[
    rownames(meta),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(rich),
      rownames(meta)
    )
  )
  
  out <- data.frame(
    SampleID = rownames(meta),
    
    Site = factor(
      meta$Site,
      levels = c("ACP", "BASS", "TRUG")
    ),
    
    Condition = factor(
      meta$Treatment,
      levels = c("NATURAL", "BESE", "BARE")
    ),
    
    Chao1 = rich$Chao1,
    Shannon = rich$Shannon,
    
    LibrarySize = sample_sums(ps),
    
    stringsAsFactors = FALSE
  )
  
  rownames(out) <- out$SampleID
  
  out
}


alpha_bac <- make_alpha_data(ps_B_bac)
alpha_fun <- make_alpha_data(ps_B_fun)


# ============================================================
# 4. BASIC ALPHA SUMMARY
# ============================================================

cat("\n========================================\n")
cat("H2 — ALPHA DATA\n")
cat("========================================\n")

cat("\nBacteria:\n")
print(
  alpha_bac %>%
    group_by(Site, Condition) %>%
    summarise(
      n = n(),
      Chao1_mean = mean(Chao1),
      Chao1_SD = sd(Chao1),
      Shannon_mean = mean(Shannon),
      Shannon_SD = sd(Shannon),
      .groups = "drop"
    )
)

cat("\nFungi:\n")
print(
  alpha_fun %>%
    group_by(Site, Condition) %>%
    summarise(
      n = n(),
      Chao1_mean = mean(Chao1),
      Chao1_SD = sd(Chao1),
      Shannon_mean = mean(Shannon),
      Shannon_SD = sd(Shannon),
      .groups = "drop"
    )
)


# ============================================================
# 5. ANALYSIS FUNCTION
#
# First tests Site × Condition interaction.
#
# If interaction p < 0.05:
#   retain full model
#   compare conditions within Site
#
# If interaction p >= 0.05:
#   fit additive Site + Condition model
#   compare overall condition means
# ============================================================

analyse_alpha <- function(data, metric, kingdom) {
  
  cat("\n\n========================================\n")
  cat(kingdom, "—", metric, "\n")
  cat("========================================\n")
  
  dat <- data[
    complete.cases(
      data[, c(metric, "Site", "Condition")]
    ),
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # FULL INTERACTION MODEL
  # ----------------------------------------------------------
  
  formula_full <- as.formula(
    paste0(
      metric,
      " ~ Site * Condition"
    )
  )
  
  mod_full <- lm(
    formula_full,
    data = dat
  )
  
  # Interaction test
  full_anova <- car::Anova(
    mod_full,
    type = 2
  )
  
  cat("\nFULL MODEL: Site × Condition\n")
  print(full_anova)
  
  interaction_p <- full_anova[
    "Site:Condition",
    "Pr(>F)"
  ]
  
  cat(
    "\nInteraction p =",
    format.pval(
      interaction_p,
      digits = 4
    ),
    "\n"
  )
  
  # ----------------------------------------------------------
  # SIGNIFICANT INTERACTION
  # ----------------------------------------------------------
  
  if (!is.na(interaction_p) &&
      interaction_p < 0.05) {
    
    cat("\n>>> Significant Site × Condition interaction\n")
    cat(">>> Retaining full interaction model\n")
    
    emm <- emmeans(
      mod_full,
      ~ Condition | Site
    )
    
    pairwise_raw <- pairs(
      emm,
      adjust = "none"
    )
    
    # Adjust across ALL site-specific pairwise contrasts
    pairwise_results <- as.data.frame(
      summary(
        pairwise_raw,
        by = NULL,
        adjust = "holm"
      )
    )
    
    cat("\nEstimated marginal means by Site:\n")
    print(emm)
    
    cat("\nCondition comparisons within Site — Holm adjusted:\n")
    print(
      pairwise_results,
      row.names = FALSE
    )
    
    final_model <- mod_full
    model_type <- "Site_x_Condition"
    
  } else {
    
    # --------------------------------------------------------
    # NON-SIGNIFICANT INTERACTION
    # --------------------------------------------------------
    
    cat("\n>>> No detectable Site × Condition interaction\n")
    cat(">>> Using additive Site + Condition model\n")
    
    formula_add <- as.formula(
      paste0(
        metric,
        " ~ Site + Condition"
      )
    )
    
    mod_add <- lm(
      formula_add,
      data = dat
    )
    
    add_anova <- car::Anova(
      mod_add,
      type = 2
    )
    
    cat("\nADDITIVE MODEL:\n")
    print(add_anova)
    
    emm <- emmeans(
      mod_add,
      ~ Condition
    )
    
    pairwise_results <- as.data.frame(
      pairs(
        emm,
        adjust = "holm"
      )
    )
    
    cat("\nEstimated marginal means:\n")
    print(emm)
    
    cat("\nOverall condition comparisons — Holm adjusted:\n")
    print(
      pairwise_results,
      row.names = FALSE
    )
    
    final_model <- mod_add
    model_type <- "Site_plus_Condition"
  }
  
  # ----------------------------------------------------------
  # LIBRARY-SIZE DIAGNOSTIC
  # ----------------------------------------------------------
  
  depth_cor <- cor.test(
    dat[[metric]],
    log10(dat$LibrarySize),
    method = "spearman",
    exact = FALSE
  )
  
  cat("\nLibrary-size diagnostic:\n")
  
  cat(
    "Spearman rho =",
    round(
      unname(depth_cor$estimate),
      3
    ),
    "; p =",
    format.pval(
      depth_cor$p.value,
      digits = 4
    ),
    "\n"
  )
  
  # ----------------------------------------------------------
  # BASIC MODEL DIAGNOSTICS
  # ----------------------------------------------------------
  
  shapiro <- shapiro.test(
    residuals(final_model)
  )
  
  levene <- car::leveneTest(
    dat[[metric]] ~ interaction(
      dat$Site,
      dat$Condition
    )
  )
  
  cat("\nResidual normality — Shapiro-Wilk:\n")
  print(shapiro)
  
  cat("\nVariance homogeneity — Levene test:\n")
  print(levene)
  
  # ----------------------------------------------------------
  # RETURN EVERYTHING
  # ----------------------------------------------------------
  
  list(
    final_model = final_model,
    model_type = model_type,
    full_interaction_model = mod_full,
    interaction_p = interaction_p,
    emmeans = emm,
    pairwise = pairwise_results,
    library_depth_test = depth_cor,
    shapiro = shapiro,
    levene = levene
  )
}


# ============================================================
# 6. RUN FOUR H2 ALPHA ANALYSES
# ============================================================

B_bac_Chao1 <- analyse_alpha(
  alpha_bac,
  "Chao1",
  "BACTERIA"
)

B_bac_Shannon <- analyse_alpha(
  alpha_bac,
  "Shannon",
  "BACTERIA"
)

B_fun_Chao1 <- analyse_alpha(
  alpha_fun,
  "Chao1",
  "FUNGI"
)

B_fun_Shannon <- analyse_alpha(
  alpha_fun,
  "Shannon",
  "FUNGI"
)


# ============================================================
# 7. FINAL COMPACT SUMMARY
# ============================================================

alpha_summary <- tibble(
  Kingdom = c(
    "Bacteria",
    "Bacteria",
    "Fungi",
    "Fungi"
  ),
  
  Metric = c(
    "Chao1",
    "Shannon",
    "Chao1",
    "Shannon"
  ),
  
  Interaction_p = c(
    B_bac_Chao1$interaction_p,
    B_bac_Shannon$interaction_p,
    B_fun_Chao1$interaction_p,
    B_fun_Shannon$interaction_p
  ),
  
  Final_model = c(
    B_bac_Chao1$model_type,
    B_bac_Shannon$model_type,
    B_fun_Chao1$model_type,
    B_fun_Shannon$model_type
  )
)

cat("\n\n========================================\n")
cat("H2 ALPHA — FINAL MODEL SELECTION\n")
cat("========================================\n")

print(
  alpha_summary,
  n = Inf
)


# ============================================================
# 8. SAVE RESULTS
# ============================================================

outdir <- "outputs/H2_alpha"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  list(
    bacteria_Chao1 = B_bac_Chao1,
    bacteria_Shannon = B_bac_Shannon,
    fungi_Chao1 = B_fun_Chao1,
    fungi_Shannon = B_fun_Shannon
  ),
  file.path(
    outdir,
    "H2_alpha_models.rds"
  )
)

write.csv(
  alpha_bac,
  file.path(
    outdir,
    "H2_bacteria_alpha_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  alpha_fun,
  file.path(
    outdir,
    "H2_fungi_alpha_data.csv"
  ),
  row.names = FALSE
)

write.csv(
  alpha_summary,
  file.path(
    outdir,
    "H2_alpha_model_selection.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "sessionInfo_H2_alpha.txt"
  )
)

cat("\n========================================\n")
cat("H2 ALPHA COMPLETE\n")
cat("========================================\n")

# ============================================================
# H2 — SAVE SUPPLEMENTARY TABLE S2
# Alpha-diversity model statistics
# No models are rerun
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(car)
})

outdir <- "outputs/H2_alpha"

# Function to extract Type-II ANOVA from the full
# Site × Condition model
extract_anova <- function(model_object, kingdom, metric) {
  
  a <- as.data.frame(
    car::Anova(
      model_object$full_interaction_model,
      type = 2
    )
  )
  
  a$Effect <- rownames(a)
  rownames(a) <- NULL
  
  a %>%
    filter(Effect != "Residuals") %>%
    transmute(
      Kingdom = kingdom,
      Metric = metric,
      Effect = recode(
        Effect,
        "Site" = "Site",
        "Condition" = "Ecological condition",
        "Site:Condition" = "Site × ecological condition"
      ),
      Df = Df,
      F = `F value`,
      p = `Pr(>F)`
    )
}


# Combine all four analyses
Supplementary_Table_S2 <- bind_rows(
  
  extract_anova(
    B_bac_Chao1,
    "Bacteria",
    "Chao1"
  ),
  
  extract_anova(
    B_bac_Shannon,
    "Bacteria",
    "Shannon"
  ),
  
  extract_anova(
    B_fun_Chao1,
    "Fungi",
    "Chao1"
  ),
  
  extract_anova(
    B_fun_Shannon,
    "Fungi",
    "Shannon"
  )
)


# Display
cat("\n========================================\n")
cat("SUPPLEMENTARY TABLE S2\n")
cat("H2 rhizosphere alpha diversity\n")
cat("========================================\n")

print(
  Supplementary_Table_S2,
  n = Inf
)


# Save
write.csv(
  Supplementary_Table_S2,
  file.path(
    outdir,
    "Supplementary_Table_S2_H2_alpha_diversity.csv"
  ),
  row.names = FALSE
)

cat("\nSaved:\n")
cat(
  normalizePath(
    file.path(
      outdir,
      "Supplementary_Table_S2_H2_alpha_diversity.csv"
    )
  ),
  "\n"
)

# ============================================================
# H2 — RHIZOSPHERE COMMUNITY REASSEMBLY
# PCoA + PERMANOVA + PERMDISP
#
# Kingdoms:
#   Bacteria
#   Fungi
#
# Conditions:
#   Natural
#   BESE
#   BARE
#
# Community data:
#   ASV-level relative abundance
#   Bray-Curtis dissimilarity
#
# PCoA:
#   descriptive only
#   Lingoes correction
#
# PERMANOVA:
#   Bray ~ Site + Condition
#   Condition tested AFTER accounting for Site
#   9,999 permutations restricted within Site
#
# PERMDISP:
#   median centroid
#   bias correction
#   same Site-restricted permutations
#
# Pairwise:
#   Natural vs BESE
#   Natural vs BARE
#   BESE vs BARE
#   Holm correction across the 3 comparisons
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(permute)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))


# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

outdir <- "outputs/H2_beta"
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 3. LOAD VERIFIED HISTORICAL PHYLOSEQ OBJECTS
# ============================================================

ps_bacteria <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)


# ============================================================
# 4. SUBSET H2
# Rhizosphere only: NATURAL, BESE, BARE
# ============================================================

ps_B_bac <- subset_samples(
  ps_bacteria,
  SoilPortion == "RZ" &
    Treatment %in% c("NATURAL", "BESE", "BARE")
)

ps_B_bac <- prune_taxa(
  taxa_sums(ps_B_bac) > 0,
  ps_B_bac
)


ps_B_fun <- subset_samples(
  ps_fungi,
  SoilPortion == "RZ" &
    Treatment %in% c("NATURAL", "BESE", "BARE")
)

ps_B_fun <- prune_taxa(
  taxa_sums(ps_B_fun) > 0,
  ps_B_fun
)


stopifnot(
  nsamples(ps_B_bac) == 45,
  nsamples(ps_B_fun) == 41
)


# ============================================================
# 5. PALETTE + SITE SHAPES
# ============================================================

condition_cols <- c(
  NATURAL = "#9BC1BC",
  BESE    = "#FF5A5F",
  BARE    = "#FFA630"
)

site_shapes <- c(
  ACP  = 21,
  BASS = 22,
  TRUG = 24
)


# ============================================================
# 6. PREPARE COMMUNITY DATA
#
# Returns:
#   raw counts: samples x ASVs
#   relative abundance
#   metadata
#   Bray-Curtis matrix
# ============================================================

prepare_beta <- function(ps) {
  
  counts <- as(
    otu_table(ps),
    "matrix"
  )
  
  if (taxa_are_rows(ps)) {
    counts <- t(counts)
  }
  
  counts <- counts[
    sample_names(ps),
    ,
    drop = FALSE
  ]
  
  # Relative abundance
  rel <- counts / rowSums(counts)
  
  # Metadata
  meta <- as(
    sample_data(ps),
    "data.frame"
  )
  
  meta <- meta[
    rownames(rel),
    ,
    drop = FALSE
  ]
  
  meta$SampleID <- rownames(meta)
  
  meta$Site <- factor(
    meta$Site,
    levels = c("ACP", "BASS", "TRUG")
  )
  
  meta$Condition <- factor(
    meta$Treatment,
    levels = c("NATURAL", "BESE", "BARE")
  )
  
  stopifnot(
    identical(
      rownames(rel),
      rownames(meta)
    )
  )
  
  # Bray-Curtis
  bray <- vegan::vegdist(
    rel,
    method = "bray"
  )
  
  list(
    counts = counts,
    rel = rel,
    meta = meta,
    bray = bray
  )
}


bac <- prepare_beta(ps_B_bac)
fun <- prepare_beta(ps_B_fun)


# ============================================================
# 7. PERMUTATION MATRIX
#
# Samples can be permuted ONLY within Site.
# ============================================================

make_site_permutations <- function(meta, nperm = 9999) {
  
  control <- permute::how(
    nperm = nperm,
    blocks = meta$Site
  )
  
  permute::shuffleSet(
    n = nrow(meta),
    nset = nperm,
    control = control
  )
}


# ============================================================
# 8. PCoA FUNCTION
#
# Bray-Curtis + Lingoes correction
# ============================================================

run_pcoa <- function(beta_object, kingdom) {
  
  ord <- vegan::wcmdscale(
    beta_object$bray,
    k = 2,
    eig = TRUE,
    add = "lingoes"
  )
  
  eig <- ord$eig
  
  positive_eig <- eig[
    eig > 0
  ]
  
  var1 <- 100 * eig[1] / sum(positive_eig)
  var2 <- 100 * eig[2] / sum(positive_eig)
  
  coords <- as.data.frame(
    ord$points[, 1:2, drop = FALSE]
  )
  
  colnames(coords) <- c(
    "PCoA1",
    "PCoA2"
  )
  
  coords$SampleID <- rownames(coords)
  
  coords <- coords %>%
    left_join(
      beta_object$meta %>%
        tibble::rownames_to_column("RowID"),
      by = c("SampleID" = "RowID")
    )
  
  coords$Condition <- factor(
    coords$Condition,
    levels = c("NATURAL", "BESE", "BARE")
  )
  
  coords$Site <- factor(
    coords$Site,
    levels = c("ACP", "BASS", "TRUG")
  )
  
  p <- ggplot(
    coords,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = Condition,
      shape = Site
    )
  ) +
    
    geom_point(
      size = 3.5,
      stroke = 0.75,
      colour = "black",
      alpha = 0.9
    ) +
    
    scale_fill_manual(
      values = condition_cols,
      labels = c(
        NATURAL = "Natural",
        BESE = "BESE",
        BARE = "BARE"
      )
    ) +
    
    scale_shape_manual(
      values = site_shapes
    ) +
    
    labs(
      x = paste0(
        "PCoA1 (",
        round(var1, 1),
        "%)"
      ),
      
      y = paste0(
        "PCoA2 (",
        round(var2, 1),
        "%)"
      ),
      
      fill = "Ecological condition",
      shape = "Site"
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      axis.title = element_text(
        size = 11,
        colour = "black"
      ),
      
      axis.text = element_text(
        size = 10,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.55
      ),
      
      legend.title = element_text(
        size = 10
      ),
      
      legend.text = element_text(
        size = 9
      ),
      
      plot.margin = margin(
        6, 6, 6, 6
      )
    )
  
  list(
    ordination = ord,
    coordinates = coords,
    variance = c(
      PCoA1 = var1,
      PCoA2 = var2
    ),
    plot = p
  )
}


pcoa_bac <- run_pcoa(
  bac,
  "Bacteria"
)

pcoa_fun <- run_pcoa(
  fun,
  "Fungi"
)


# ============================================================
# 9. GLOBAL PERMANOVA
#
# Condition is interpreted AFTER accounting for Site.
#
# Site R2 can be reported descriptively.
# Site p-value should NOT be interpreted because
# permutations are restricted within Site.
# ============================================================

run_permanova <- function(beta_object, kingdom) {
  
  meta <- beta_object$meta
  
  perms <- make_site_permutations(
    meta,
    nperm = 9999
  )
  
  fit <- vegan::adonis2(
    beta_object$bray ~ Site + Condition,
    data = meta,
    permutations = perms,
    by = "terms"
  )
  
  tab <- as.data.frame(fit)
  
  tab$Effect <- rownames(tab)
  rownames(tab) <- NULL
  
  # Clean names
  tab$Effect_clean <- ifelse(
    tab$Effect == "Site",
    "Site",
    ifelse(
      tab$Effect == "Condition",
      "Ecological condition",
      tab$Effect
    )
  )
  
  # Whole-model R2
  term_rows <- tab$Effect %in% c(
    "Site",
    "Condition"
  )
  
  model_R2 <- sum(
    tab$R2[term_rows],
    na.rm = TRUE
  )
  
  condition_row <- tab[
    tab$Effect == "Condition",
    ,
    drop = FALSE
  ]
  
  site_row <- tab[
    tab$Effect == "Site",
    ,
    drop = FALSE
  ]
  
  summary <- tibble(
    Kingdom = kingdom,
    
    Samples = nrow(meta),
    
    Site_R2 = site_row$R2,
    
    Condition_F = condition_row$F,
    
    Condition_R2 = condition_row$R2,
    
    Condition_p = condition_row$`Pr(>F)`,
    
    Whole_model_R2 = model_R2
  )
  
  list(
    fit = fit,
    table = tab,
    summary = summary,
    permutations = perms
  )
}


perm_bac <- run_permanova(
  bac,
  "Bacteria"
)

perm_fun <- run_permanova(
  fun,
  "Fungi"
)


# ============================================================
# 10. GLOBAL PERMDISP
# ============================================================

run_permdisp <- function(beta_object, kingdom) {
  
  meta <- beta_object$meta
  
  perms <- make_site_permutations(
    meta,
    nperm = 9999
  )
  
  bd <- vegan::betadisper(
    beta_object$bray,
    group = meta$Condition,
    type = "median",
    bias.adjust = TRUE
  )
  
  test <- vegan::permutest(
    bd,
    permutations = perms
  )
  
  tab <- as.data.frame(
    test$tab
  )
  
  # Mean distance to centroid by condition
  distances <- tibble(
    SampleID = names(bd$distances),
    Distance_to_centroid = as.numeric(
      bd$distances
    ),
    Condition = meta$Condition
  )
  
  means <- distances %>%
    group_by(Condition) %>%
    summarise(
      n = n(),
      Mean_distance = mean(Distance_to_centroid),
      SD_distance = sd(Distance_to_centroid),
      .groups = "drop"
    )
  
  summary <- tibble(
    Kingdom = kingdom,
    F = tab$F[1],
    p = tab$`Pr(>F)`[1]
  )
  
  list(
    model = bd,
    test = test,
    table = tab,
    distances = distances,
    means = means,
    summary = summary
  )
}


disp_bac <- run_permdisp(
  bac,
  "Bacteria"
)

disp_fun <- run_permdisp(
  fun,
  "Fungi"
)


# ============================================================
# 11. PAIRWISE PERMANOVA + PERMDISP
# ============================================================

run_pairwise_beta <- function(
    beta_object,
    kingdom
) {
  
  pairs <- list(
    c("NATURAL", "BESE"),
    c("NATURAL", "BARE"),
    c("BESE", "BARE")
  )
  
  results <- lapply(
    pairs,
    function(pair) {
      
      keep <- beta_object$meta$Condition %in% pair
      
      meta_sub <- droplevels(
        beta_object$meta[
          keep,
          ,
          drop = FALSE
        ]
      )
      
      rel_sub <- beta_object$rel[
        keep,
        ,
        drop = FALSE
      ]
      
      # Remove ASVs absent from this pair
      rel_sub <- rel_sub[
        ,
        colSums(rel_sub) > 0,
        drop = FALSE
      ]
      
      bray_sub <- vegan::vegdist(
        rel_sub,
        method = "bray"
      )
      
      perms <- make_site_permutations(
        meta_sub,
        nperm = 9999
      )
      
      
      # ------------------------------------------------------
      # PAIRWISE PERMANOVA
      # ------------------------------------------------------
      
      ad <- vegan::adonis2(
        bray_sub ~ Site + Condition,
        data = meta_sub,
        permutations = perms,
        by = "terms"
      )
      
      ad_df <- as.data.frame(ad)
      
      condition_row <- ad_df[
        rownames(ad_df) == "Condition",
        ,
        drop = FALSE
      ]
      
      
      # ------------------------------------------------------
      # PAIRWISE PERMDISP
      # ------------------------------------------------------
      
      bd <- vegan::betadisper(
        bray_sub,
        group = meta_sub$Condition,
        type = "median",
        bias.adjust = TRUE
      )
      
      pd <- vegan::permutest(
        bd,
        permutations = perms
      )
      
      pd_tab <- as.data.frame(
        pd$tab
      )
      
      
      # ------------------------------------------------------
      # RESULTS
      # ------------------------------------------------------
      
      tibble(
        Kingdom = kingdom,
        
        Contrast = paste(
          pair[1],
          "vs",
          pair[2]
        ),
        
        n = nrow(meta_sub),
        
        PERMANOVA_F = condition_row$F,
        
        PERMANOVA_R2 = condition_row$R2,
        
        PERMANOVA_p = condition_row$`Pr(>F)`,
        
        PERMDISP_F = pd_tab$F[1],
        
        PERMDISP_p = pd_tab$`Pr(>F)`[1]
      )
    }
  )
  
  bind_rows(results) %>%
    mutate(
      PERMANOVA_p_Holm = p.adjust(
        PERMANOVA_p,
        method = "holm"
      ),
      
      PERMDISP_p_Holm = p.adjust(
        PERMDISP_p,
        method = "holm"
      )
    )
}


pairwise_bac <- run_pairwise_beta(
  bac,
  "Bacteria"
)

pairwise_fun <- run_pairwise_beta(
  fun,
  "Fungi"
)


# ============================================================
# 12. GLOBAL SUMMARY
# ============================================================

global_PERMANOVA <- bind_rows(
  perm_bac$summary,
  perm_fun$summary
)

global_PERMDISP <- bind_rows(
  disp_bac$summary,
  disp_fun$summary
)

pairwise_beta <- bind_rows(
  pairwise_bac,
  pairwise_fun
)


# ============================================================
# 13. PRINT RESULTS
# ============================================================

cat("\n\n========================================\n")
cat("H2 — GLOBAL PERMANOVA\n")
cat("========================================\n")

print(
  global_PERMANOVA,
  n = Inf
)


cat("\n\n========================================\n")
cat("H2 — GLOBAL PERMDISP\n")
cat("========================================\n")

print(
  global_PERMDISP,
  n = Inf
)


cat("\n\n========================================\n")
cat("H2 — PAIRWISE PERMANOVA + PERMDISP\n")
cat("========================================\n")

print(
  pairwise_beta,
  n = Inf
)


cat("\n\n========================================\n")
cat("MEAN DISTANCE TO CENTROID — BACTERIA\n")
cat("========================================\n")

print(
  disp_bac$means,
  n = Inf
)


cat("\n\n========================================\n")
cat("MEAN DISTANCE TO CENTROID — FUNGI\n")
cat("========================================\n")

print(
  disp_fun$means,
  n = Inf
)


# ============================================================
# 14. COMBINE PCoA FIGURE
#
# Treatment colour + Site shape
# No statistical inference shown on ordination itself.
# ============================================================

p_bac <- pcoa_bac$plot +
  labs(
    title = "Bacteria"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    )
  )

p_fun <- pcoa_fun$plot +
  labs(
    title = "Fungi"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    )
  )


fig_H2_PCoA <- (
  p_bac +
    p_fun
) +
  
  plot_layout(
    guides = "collect"
  ) +
  
  plot_annotation(
    tag_levels = "a"
  ) &
  
  theme(
    legend.position = "bottom",
    
    plot.tag = element_text(
      face = "bold",
      size = 12
    )
  )


print(fig_H2_PCoA)


# ============================================================
# 15. SAVE FIGURE
# ============================================================

ggsave(
  file.path(
    figdir,
    "H2_Rhizosphere_PCoA_Bacteria_Fungi.pdf"
  ),
  fig_H2_PCoA,
  width = 8.2,
  height = 4.3
)

ggsave(
  file.path(
    figdir,
    "H2_Rhizosphere_PCoA_Bacteria_Fungi.png"
  ),
  fig_H2_PCoA,
  width = 8.2,
  height = 4.3,
  dpi = 600
)


# ============================================================
# 16. SAVE PCoA COORDINATES
# ============================================================

write.csv(
  pcoa_bac$coordinates,
  file.path(
    tabdir,
    "H2_Bacteria_PCoA_coordinates.csv"
  ),
  row.names = FALSE
)

write.csv(
  pcoa_fun$coordinates,
  file.path(
    tabdir,
    "H2_Fungi_PCoA_coordinates.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. SAVE PERMANOVA TABLES
# ============================================================

write.csv(
  perm_bac$table,
  file.path(
    tabdir,
    "H2_Bacteria_PERMANOVA_full.csv"
  ),
  row.names = FALSE
)

write.csv(
  perm_fun$table,
  file.path(
    tabdir,
    "H2_Fungi_PERMANOVA_full.csv"
  ),
  row.names = FALSE
)

write.csv(
  global_PERMANOVA,
  file.path(
    tabdir,
    "H2_PERMANOVA_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. SAVE PERMDISP TABLES
# ============================================================

write.csv(
  global_PERMDISP,
  file.path(
    tabdir,
    "H2_PERMDISP_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  disp_bac$means,
  file.path(
    tabdir,
    "H2_Bacteria_PERMDISP_distances_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  disp_fun$means,
  file.path(
    tabdir,
    "H2_Fungi_PERMDISP_distances_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. SAVE PAIRWISE RESULTS
# ============================================================

write.csv(
  pairwise_beta,
  file.path(
    tabdir,
    "H2_pairwise_PERMANOVA_PERMDISP.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 20. SAVE FULL MODEL OBJECTS
# ============================================================

saveRDS(
  list(
    bacteria = list(
      PCoA = pcoa_bac,
      PERMANOVA = perm_bac,
      PERMDISP = disp_bac
    ),
    
    fungi = list(
      PCoA = pcoa_fun,
      PERMANOVA = perm_fun,
      PERMDISP = disp_fun
    ),
    
    pairwise = pairwise_beta
  ),
  
  file.path(
    outdir,
    "H2_beta_full_objects.rds"
  )
)


# ============================================================
# 21. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "sessionInfo_H2_beta.txt"
  )
)


# ============================================================
# 22. FINAL COMPACT SUMMARY
# ============================================================

cat("\n\n========================================\n")
cat("H2 BETA — FINAL SUMMARY\n")
cat("========================================\n")

cat("\nGLOBAL PERMANOVA:\n")
print(
  global_PERMANOVA,
  n = Inf
)

cat("\nGLOBAL PERMDISP:\n")
print(
  global_PERMDISP,
  n = Inf
)

cat("\nPAIRWISE TESTS:\n")
print(
  pairwise_beta,
  n = Inf
)

cat("\nPCoA variance explained:\n")

cat(
  "Bacteria:",
  "PCoA1 =",
  round(pcoa_bac$variance["PCoA1"], 1),
  "%;",
  "PCoA2 =",
  round(pcoa_bac$variance["PCoA2"], 1),
  "%\n"
)

cat(
  "Fungi:",
  "PCoA1 =",
  round(pcoa_fun$variance["PCoA1"], 1),
  "%;",
  "PCoA2 =",
  round(pcoa_fun$variance["PCoA2"], 1),
  "%\n"
)

cat("\nResults saved in:\n")
cat(normalizePath(outdir), "\n")

cat("\n========================================\n")
cat("H2 BETA COMPLETE\n")
cat("========================================\n")





#DISTANCE TO NATURAL
# ============================================================
# H2 — RHIZOSPHERE REASSEMBLY
# DISTANCE TO NATURAL
#
# Research question:
# Are BESE rhizosphere communities more Natural-like than
# BARE rhizosphere communities?
#
# Community data:
#   ASV-level relative abundance
#   Bray-Curtis dissimilarity
#
# Natural reference:
#   same-site Natural samples
#
# For Natural samples:
#   leave-one-out mean distance to the other Natural
#   samples from the same Site
#
# PRIMARY INFERENCE:
#   BESE vs BARE within matched restoration blocks
#
# Natural is shown descriptively as a reference baseline,
# but is NOT treated as a randomized restoration treatment.
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))


# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

outdir <- "outputs/H2_beta/distance_to_Natural"
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 3. LOAD VERIFIED PHYLOSEQ OBJECTS
# ============================================================

ps_bacteria <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)


# ============================================================
# 4. SUBSET H2
# ============================================================

ps_B_bac <- subset_samples(
  ps_bacteria,
  SoilPortion == "RZ" &
    Treatment %in% c("NATURAL", "BESE", "BARE")
)

ps_B_bac <- prune_taxa(
  taxa_sums(ps_B_bac) > 0,
  ps_B_bac
)


ps_B_fun <- subset_samples(
  ps_fungi,
  SoilPortion == "RZ" &
    Treatment %in% c("NATURAL", "BESE", "BARE")
)

ps_B_fun <- prune_taxa(
  taxa_sums(ps_B_fun) > 0,
  ps_B_fun
)

stopifnot(
  nsamples(ps_B_bac) == 45,
  nsamples(ps_B_fun) == 41
)


# ============================================================
# 5. FUNCTION — DISTANCE TO SAME-SITE NATURAL
# ============================================================

calculate_distance_to_natural <- function(ps, kingdom) {
  
  # ----------------------------------------------------------
  # Raw ASV matrix
  # samples x ASVs
  # ----------------------------------------------------------
  
  counts <- as(
    otu_table(ps),
    "matrix"
  )
  
  if (taxa_are_rows(ps)) {
    counts <- t(counts)
  }
  
  counts <- counts[
    sample_names(ps),
    ,
    drop = FALSE
  ]
  
  # ----------------------------------------------------------
  # Relative abundance
  # ----------------------------------------------------------
  
  rel <- counts / rowSums(counts)
  
  # ----------------------------------------------------------
  # Bray-Curtis
  # ----------------------------------------------------------
  
  bray <- as.matrix(
    vegan::vegdist(
      rel,
      method = "bray"
    )
  )
  
  # ----------------------------------------------------------
  # Metadata
  # ----------------------------------------------------------
  
  meta <- as(
    sample_data(ps),
    "data.frame"
  )
  
  meta <- meta[
    rownames(rel),
    ,
    drop = FALSE
  ]
  
  meta$SampleID <- rownames(meta)
  
  meta$Site <- factor(
    meta$Site,
    levels = c("ACP", "BASS", "TRUG")
  )
  
  meta$Condition <- factor(
    meta$Treatment,
    levels = c("NATURAL", "BESE", "BARE")
  )
  
  # ----------------------------------------------------------
  # Check Natural references per Site
  # ----------------------------------------------------------
  
  natural_counts <- meta %>%
    filter(Condition == "NATURAL") %>%
    count(Site, name = "Natural_n")
  
  cat("\n========================================\n")
  cat(kingdom, "— NATURAL REFERENCE COUNTS\n")
  cat("========================================\n")
  
  print(natural_counts)
  
  if (any(natural_counts$Natural_n < 2)) {
    stop(
      "At least two Natural samples per Site are required."
    )
  }
  
  # ----------------------------------------------------------
  # Calculate mean Bray-Curtis distance to same-site Natural
  # ----------------------------------------------------------
  
  distance <- numeric(
    nrow(meta)
  )
  
  names(distance) <- meta$SampleID
  
  for (i in seq_len(nrow(meta))) {
    
    sample_i <- meta$SampleID[i]
    site_i   <- meta$Site[i]
    cond_i   <- meta$Condition[i]
    
    natural_ids <- meta$SampleID[
      meta$Site == site_i &
        meta$Condition == "NATURAL"
    ]
    
    # Natural samples:
    # exclude self from the reference set
    if (cond_i == "NATURAL") {
      
      natural_ids <- setdiff(
        natural_ids,
        sample_i
      )
    }
    
    if (length(natural_ids) == 0) {
      stop(
        paste(
          "No valid Natural references for sample",
          sample_i
        )
      )
    }
    
    distance[i] <- mean(
      bray[
        sample_i,
        natural_ids,
        drop = TRUE
      ]
    )
  }
  
  # ----------------------------------------------------------
  # Final sample-level table
  # ----------------------------------------------------------
  
  out <- meta %>%
    mutate(
      Kingdom = kingdom,
      
      Distance_to_Natural =
        distance[SampleID]
    )
  
  out
}


# ============================================================
# 6. CALCULATE DISTANCES
# ============================================================

dist_bac <- calculate_distance_to_natural(
  ps_B_bac,
  "Bacteria"
)

dist_fun <- calculate_distance_to_natural(
  ps_B_fun,
  "Fungi"
)


# ============================================================
# 7. DESCRIPTIVE SUMMARY
# ============================================================

make_descriptive <- function(dat) {
  
  dat %>%
    group_by(
      Kingdom,
      Condition
    ) %>%
    summarise(
      n = n(),
      
      Mean = mean(
        Distance_to_Natural
      ),
      
      SD = sd(
        Distance_to_Natural
      ),
      
      Median = median(
        Distance_to_Natural
      ),
      
      Q1 = quantile(
        Distance_to_Natural,
        0.25
      ),
      
      Q3 = quantile(
        Distance_to_Natural,
        0.75
      ),
      
      .groups = "drop"
    )
}


descriptive <- bind_rows(
  make_descriptive(dist_bac),
  make_descriptive(dist_fun)
)

cat("\n========================================\n")
cat("DISTANCE TO NATURAL — DESCRIPTIVE SUMMARY\n")
cat("========================================\n")

print(
  descriptive,
  n = Inf
)


# ============================================================
# 8. IDENTIFY COMPLETE BESE–BARE BLOCK PAIRS
#
# Natural has no restoration Block and is NOT part of
# this inferential analysis.
# ============================================================

prepare_paired_data <- function(dat, kingdom) {
  
  restoration <- dat %>%
    filter(
      Condition %in% c(
        "BESE",
        "BARE"
      ),
      !is.na(Block)
    ) %>%
    mutate(
      PairID = paste(
        Site,
        Block,
        sep = "_"
      )
    )
  
  pair_check <- restoration %>%
    group_by(
      Site,
      Block,
      PairID
    ) %>%
    summarise(
      n_samples = n(),
      
      n_conditions =
        n_distinct(Condition),
      
      Has_BESE =
        any(Condition == "BESE"),
      
      Has_BARE =
        any(Condition == "BARE"),
      
      .groups = "drop"
    )
  
  complete_pairs <- pair_check %>%
    filter(
      Has_BESE,
      Has_BARE,
      n_conditions == 2
    )
  
  paired <- restoration %>%
    semi_join(
      complete_pairs,
      by = c(
        "Site",
        "Block",
        "PairID"
      )
    )
  
  cat("\n========================================\n")
  cat(kingdom, "— BLOCK PAIRING\n")
  cat("========================================\n")
  
  cat(
    "Potential restoration blocks:",
    nrow(pair_check),
    "\n"
  )
  
  cat(
    "Complete BESE–BARE pairs:",
    nrow(complete_pairs),
    "\n"
  )
  
  cat(
    "Samples entering paired analysis:",
    nrow(paired),
    "\n"
  )
  
  cat("\nComplete pairs by Site:\n")
  
  print(
    table(
      complete_pairs$Site
    )
  )
  
  list(
    paired = paired,
    pair_check = pair_check,
    complete_pairs = complete_pairs
  )
}


paired_bac <- prepare_paired_data(
  dist_bac,
  "BACTERIA"
)

paired_fun <- prepare_paired_data(
  dist_fun,
  "FUNGI"
)


# ============================================================
# 9. PAIRED INFERENTIAL ANALYSIS
#
# Equivalent to testing the within-block BESE - BARE
# differences.
#
# Positive difference:
#   BESE is FARTHER from Natural than BARE
#
# Negative difference:
#   BESE is CLOSER to Natural than BARE
# ============================================================

run_paired_test <- function(
    paired_object,
    kingdom
) {
  
  dat <- paired_object$paired
  
  wide <- dat %>%
    select(
      Site,
      Block,
      PairID,
      Condition,
      Distance_to_Natural
    ) %>%
    pivot_wider(
      names_from = Condition,
      values_from = Distance_to_Natural
    ) %>%
    filter(
      !is.na(BESE),
      !is.na(BARE)
    ) %>%
    mutate(
      Difference_BESE_minus_BARE =
        BESE - BARE
    )
  
  # ----------------------------------------------------------
  # Primary paired t-test
  # ----------------------------------------------------------
  
  test <- t.test(
    wide$BESE,
    wide$BARE,
    paired = TRUE
  )
  
  # ----------------------------------------------------------
  # Difference normality diagnostic
  # ----------------------------------------------------------
  
  shap <- shapiro.test(
    wide$Difference_BESE_minus_BARE
  )
  
  # ----------------------------------------------------------
  # Wilcoxon paired sensitivity analysis
  # ----------------------------------------------------------
  
  wilcox <- wilcox.test(
    wide$BESE,
    wide$BARE,
    paired = TRUE,
    exact = FALSE
  )
  
  summary <- tibble(
    Kingdom = kingdom,
    
    Complete_pairs = nrow(wide),
    
    Mean_BESE =
      mean(wide$BESE),
    
    Mean_BARE =
      mean(wide$BARE),
    
    Mean_difference_BESE_minus_BARE =
      mean(
        wide$Difference_BESE_minus_BARE
      ),
    
    SE_difference =
      sd(
        wide$Difference_BESE_minus_BARE
      ) /
      sqrt(nrow(wide)),
    
    t =
      unname(test$statistic),
    
    df =
      unname(test$parameter),
    
    p =
      test$p.value,
    
    CI_low =
      test$conf.int[1],
    
    CI_high =
      test$conf.int[2],
    
    Shapiro_W =
      unname(shap$statistic),
    
    Shapiro_p =
      shap$p.value,
    
    Wilcoxon_V =
      unname(wilcox$statistic),
    
    Wilcoxon_p =
      wilcox$p.value
  )
  
  cat("\n========================================\n")
  cat(kingdom, "— DISTANCE TO NATURAL\n")
  cat("BESE vs BARE — PAIRED TEST\n")
  cat("========================================\n")
  
  print(
    summary,
    n = Inf,
    width = Inf
  )
  
  cat("\nWithin-block differences by Site:\n")
  
  print(
    wide %>%
      group_by(Site) %>%
      summarise(
        n = n(),
        
        Mean_difference =
          mean(
            Difference_BESE_minus_BARE
          ),
        
        SD_difference =
          sd(
            Difference_BESE_minus_BARE
          ),
        
        .groups = "drop"
      )
  )
  
  list(
    wide = wide,
    test = test,
    shapiro = shap,
    wilcox = wilcox,
    summary = summary
  )
}


test_bac <- run_paired_test(
  paired_bac,
  "Bacteria"
)

test_fun <- run_paired_test(
  paired_fun,
  "Fungi"
)


# ============================================================
# 10. COMBINE TEST SUMMARY
# ============================================================

test_summary <- bind_rows(
  test_bac$summary,
  test_fun$summary
)

cat("\n\n========================================\n")
cat("H2 — DISTANCE TO NATURAL FINAL SUMMARY\n")
cat("========================================\n")

print(
  test_summary,
  n = Inf,
  width = Inf
)


# ============================================================
# 11. COLOURS
# ============================================================

condition_cols <- c(
  NATURAL = "#9BC1BC",
  BESE    = "#FF5A5F",
  BARE    = "#FFA630"
)


# ============================================================
# 12. PLOT FUNCTION
#
# Natural is descriptive reference.
# Significance bracket is ONLY BESE vs BARE.
# ============================================================

make_distance_plot <- function(
    dat,
    test_object,
    kingdom
) {
  
  dat$Condition <- factor(
    dat$Condition,
    levels = c(
      "NATURAL",
      "BESE",
      "BARE"
    )
  )
  
  pval <- test_object$summary$p
  
  sig_label <- ifelse(
    pval < 0.001,
    "***",
    ifelse(
      pval < 0.01,
      "**",
      ifelse(
        pval < 0.05,
        "*",
        "ns"
      )
    )
  )
  
  ymax <- max(
    dat$Distance_to_Natural,
    na.rm = TRUE
  )
  
  ymin <- min(
    dat$Distance_to_Natural,
    na.rm = TRUE
  )
  
  yrange <- ymax - ymin
  
  bracket_y <-
    ymax + 0.08 * yrange
  
  p <- ggplot(
    dat,
    aes(
      x = Condition,
      y = Distance_to_Natural,
      fill = Condition
    )
  ) +
    
    geom_boxplot(
      width = 0.58,
      linewidth = 0.65,
      alpha = 0.75,
      outlier.shape = NA,
      colour = "black"
    ) +
    
    geom_jitter(
      width = 0.09,
      shape = 21,
      size = 2.5,
      stroke = 0.45,
      colour = "black",
      alpha = 0.9
    ) +
    
    # BESE vs BARE bracket only
    geom_segment(
      aes(
        x = 2,
        xend = 3,
        y = bracket_y,
        yend = bracket_y
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    
    geom_segment(
      aes(
        x = 2,
        xend = 2,
        y = bracket_y,
        yend =
          bracket_y -
          0.025 * yrange
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    
    geom_segment(
      aes(
        x = 3,
        xend = 3,
        y = bracket_y,
        yend =
          bracket_y -
          0.025 * yrange
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    
    annotate(
      "text",
      x = 2.5,
      y =
        bracket_y +
        0.035 * yrange,
      label = sig_label,
      size = 4.2
    ) +
    
    scale_fill_manual(
      values = condition_cols,
      labels = c(
        NATURAL = "Natural",
        BESE = "BESE",
        BARE = "BARE"
      )
    ) +
    
    scale_x_discrete(
      labels = c(
        NATURAL = "Natural",
        BESE = "BESE",
        BARE = "BARE"
      )
    ) +
    
    scale_y_continuous(
      expand = expansion(
        mult = c(
          0.04,
          0.18
        )
      )
    ) +
    
    labs(
      title = kingdom,
      x = NULL,
      y = "Bray–Curtis distance to Natural"
    ) +
    
    theme_classic(
      base_size = 10.5
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 11
      ),
      
      axis.title = element_text(
        size = 11
      ),
      
      axis.text = element_text(
        size = 10,
        colour = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.55
      ),
      
      legend.position = "none"
    )
  
  p
}


# ============================================================
# 13. CREATE PLOTS
# ============================================================

p_bac_dist <- make_distance_plot(
  dist_bac,
  test_bac,
  "Bacteria"
)

p_fun_dist <- make_distance_plot(
  dist_fun,
  test_fun,
  "Fungi"
)


fig_distance <- (
  p_bac_dist +
    p_fun_dist
) +
  
  plot_annotation(
    tag_levels = "a"
  ) &
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 12
    )
  )


print(fig_distance)


# ============================================================
# 14. SAVE FIGURE
# ============================================================

ggsave(
  file.path(
    figdir,
    "H2_Distance_to_Natural.pdf"
  ),
  fig_distance,
  width = 6.7,
  height = 3.8
)

ggsave(
  file.path(
    figdir,
    "H2_Distance_to_Natural.png"
  ),
  fig_distance,
  width = 6.7,
  height = 3.8,
  dpi = 600
)


# ============================================================
# 15. SAVE SAMPLE-LEVEL DISTANCES
# ============================================================

write.csv(
  dist_bac,
  file.path(
    tabdir,
    "H2_Bacteria_distance_to_Natural_samples.csv"
  ),
  row.names = FALSE
)

write.csv(
  dist_fun,
  file.path(
    tabdir,
    "H2_Fungi_distance_to_Natural_samples.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. SAVE DESCRIPTIVE SUMMARY
# ============================================================

write.csv(
  descriptive,
  file.path(
    tabdir,
    "H2_distance_to_Natural_descriptive_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 17. SAVE PAIRED DATA
# ============================================================

write.csv(
  test_bac$wide,
  file.path(
    tabdir,
    "H2_Bacteria_BESE_BARE_complete_pairs.csv"
  ),
  row.names = FALSE
)

write.csv(
  test_fun$wide,
  file.path(
    tabdir,
    "H2_Fungi_BESE_BARE_complete_pairs.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 18. SAVE INFERENTIAL SUMMARY
# ============================================================

write.csv(
  test_summary,
  file.path(
    tabdir,
    "H2_distance_to_Natural_BESE_vs_BARE_tests.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 19. SAVE OBJECTS
# ============================================================

saveRDS(
  list(
    bacteria = list(
      sample_distances = dist_bac,
      paired_analysis = test_bac
    ),
    
    fungi = list(
      sample_distances = dist_fun,
      paired_analysis = test_fun
    )
  ),
  
  file.path(
    outdir,
    "H2_distance_to_Natural_full_objects.rds"
  )
)


# ============================================================
# 20. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file = file.path(
    outdir,
    "sessionInfo_H2_distance_to_Natural.txt"
  )
)


# ============================================================
# 21. FINAL INTERPRETATION GUIDE
# ============================================================

cat("\n\n========================================\n")
cat("INTERPRETATION\n")
cat("========================================\n")

cat(
  "\nMean difference = BESE - BARE\n"
)

cat(
  "\nNegative value:\n",
  "BESE is closer to the Natural reference than BARE.\n",
  sep = ""
)

cat(
  "\nPositive value:\n",
  "BARE is closer to the Natural reference than BESE.\n",
  sep = ""
)

cat(
  "\nNon-significant result:\n",
  "No evidence that BESE and BARE differ in their ",
  "compositional resemblance to Natural.\n",
  sep = ""
)

cat(
  "\nIMPORTANT:\n",
  "Because there is no Degraded rhizosphere reference, ",
  "this analysis compares BESE with BARE but does not, ",
  "by itself, demonstrate recovery toward Natural.\n",
  sep = ""
)

cat("\n========================================\n")
cat("H2 DISTANCE TO NATURAL COMPLETE\n")
cat("========================================\n")


# ============================================================
# H2 ANCOM-BC2
# RECOVER BACTERIA GENUS + FIX EMPTY-CONTRAST BUG
# THEN CONTINUE ASV / FUNGI ANALYSES
#
# IMPORTANT:
# Bacteria genus ANCOM-BC2 is NOT rerun.
# The already-saved fitted object is recovered.
# ============================================================

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# ============================================================
# 1. FIX PAIRWISE EXTRACTION FUNCTION
#
# Main fix:
# force Higher to CHARACTER even when there are zero rows.
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
        "Missing ANCOM-BC2 pairwise columns for",
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
      taxon = as.character(taxon),
      
      Contrast = as.character(
        contrast_label
      ),
      
      LFC = as.numeric(
        .data[[lfc_col]]
      ),
      
      q = as.numeric(
        .data[[q_col]]
      ),
      
      Higher = as.character(
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
    arrange(q)
  
  tibble::as_tibble(out)
}


# ============================================================
# 2. RECOVER SAVED BACTERIA-GENUS MODEL
# ============================================================

bac_genus_dir <-
  "outputs/H2_ANCOMBC2/bacteria_genus"

fit <- readRDS(
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_ANCOMBC2_full_object.rds"
  )
)

map <- read.csv(
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_taxon_map.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

map <- tibble::as_tibble(map)

cat("\n========================================\n")
cat("RECOVERED SAVED BACTERIA GENUS MODEL\n")
cat("========================================\n")

cat(
  "Quantitatively modelled genera:",
  nrow(fit$res),
  "\n"
)


# ============================================================
# 3. PRIMARY RESULTS
# ============================================================

primary <- as.data.frame(
  fit$res
)

primary_named <- tibble::as_tibble(
  primary
) %>%
  left_join(
    map,
    by = "taxon"
  )

write.csv(
  primary_named,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_primary_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 4. GLOBAL RESULTS
# ============================================================

global <- as.data.frame(
  fit$res_global
)

global_named <- tibble::as_tibble(
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
    bac_genus_dir,
    "H2_Bacteria_Genus_global_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  global_robust,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_global_ROBUST.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 5. STRUCTURAL ZEROS
# ============================================================

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
    function(x) any(x %in% TRUE)
  )
  
} else {
  
  has_structural_zero <- rep(
    FALSE,
    nrow(zero_ind)
  )
}

zero_named <- tibble::as_tibble(
  zero_ind
) %>%
  left_join(
    map,
    by = "taxon"
  )

structural_zero_taxa <- zero_named[
  has_structural_zero,
  ,
  drop = FALSE
]

write.csv(
  zero_named,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_structural_zero_all.csv"
  ),
  row.names = FALSE
)

write.csv(
  structural_zero_taxa,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_structural_zero_taxa.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. PAIRWISE RESULTS
# ============================================================

pair <- as.data.frame(
  fit$res_pair
)

pair_named <- tibble::as_tibble(
  pair
) %>%
  left_join(
    map,
    by = "taxon"
  )

write.csv(
  pair_named,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_pairwise_all.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 7. EXTRACT ROBUST CONTRASTS
# ============================================================

BESE_vs_NATURAL <- extract_robust_contrast(
  pair_table = pair,
  map = map,
  term = "ConditionBESE",
  contrast_label = "BESE - NATURAL",
  positive_label = "BESE",
  negative_label = "NATURAL"
)

BARE_vs_NATURAL <- extract_robust_contrast(
  pair_table = pair,
  map = map,
  term = "ConditionBARE",
  contrast_label = "BARE - NATURAL",
  positive_label = "BARE",
  negative_label = "NATURAL"
)

BARE_vs_BESE <- extract_robust_contrast(
  pair_table = pair,
  map = map,
  term = "ConditionBARE_ConditionBESE",
  contrast_label = "BARE - BESE",
  positive_label = "BARE",
  negative_label = "BESE"
)


# ============================================================
# 8. SAVE ROBUST TABLES
# ============================================================

write.csv(
  BESE_vs_NATURAL,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_BESE_vs_NATURAL_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  BARE_vs_NATURAL,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_BARE_vs_NATURAL_ROBUST.csv"
  ),
  row.names = FALSE
)

write.csv(
  BARE_vs_BESE,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_BARE_vs_BESE_ROBUST.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. PAIRWISE SUMMARY
# ============================================================

pairwise_summary <- tibble(
  Kingdom = "Bacteria",
  
  Level = "Genus",
  
  Contrast = c(
    "BESE - NATURAL",
    "BARE - NATURAL",
    "BARE - BESE"
  ),
  
  Robust_taxa = c(
    nrow(BESE_vs_NATURAL),
    nrow(BARE_vs_NATURAL),
    nrow(BARE_vs_BESE)
  ),
  
  Higher_NATURAL = c(
    sum(
      BESE_vs_NATURAL$Higher == "NATURAL",
      na.rm = TRUE
    ),
    
    sum(
      BARE_vs_NATURAL$Higher == "NATURAL",
      na.rm = TRUE
    ),
    
    0L
  ),
  
  Higher_BESE = c(
    sum(
      BESE_vs_NATURAL$Higher == "BESE",
      na.rm = TRUE
    ),
    
    0L,
    
    sum(
      BARE_vs_BESE$Higher == "BESE",
      na.rm = TRUE
    )
  ),
  
  Higher_BARE = c(
    0L,
    
    sum(
      BARE_vs_NATURAL$Higher == "BARE",
      na.rm = TRUE
    ),
    
    sum(
      BARE_vs_BESE$Higher == "BARE",
      na.rm = TRUE
    )
  )
)


# ============================================================
# 10. COMBINE ROBUST RESULTS
#
# This is where the old script failed.
# It now works because Higher is always character.
# ============================================================

all_pairwise_robust <- bind_rows(
  BESE_vs_NATURAL,
  BARE_vs_NATURAL,
  BARE_vs_BESE
)

write.csv(
  all_pairwise_robust,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_all_ROBUST_pairwise_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  pairwise_summary,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_pairwise_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 11. STRUCTURAL-ZERO SUMMARY
# ============================================================

structural_zero_summary <- tibble(
  Kingdom = "Bacteria",
  Level = "Genus",
  
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

write.csv(
  structural_zero_summary,
  file.path(
    bac_genus_dir,
    "H2_Bacteria_Genus_structural_zero_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# H1 — RECOVER FINAL ROBUST GENUS TABLES
# AND CHECK BACTERIA vs ARCHAEA
#
# NO ANALYSIS IS RERUN
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})


# ============================================================
# 1. USE ONLY FINAL H1 PAIRWISE TABLES
# ============================================================

H1_pairwise_files <- robust_H1[
  grepl(
    "bacteria_genus_FINAL",
    robust_H1,
    fixed = TRUE
  ) &
    !grepl(
      "global",
      robust_H1,
      ignore.case = TRUE
    )
]


cat("\n========================================\n")
cat("FINAL H1 PAIRWISE FILES\n")
cat("========================================\n")

print(H1_pairwise_files)


# ============================================================
# 2. READ NON-EMPTY TABLES
# ============================================================

H1_list <- list()

cat("\n========================================\n")
cat("H1 PAIRWISE ROW COUNTS\n")
cat("========================================\n")

for (f in H1_pairwise_files) {
  
  x <- read.csv(
    f,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  cat(
    basename(f),
    ":",
    nrow(x),
    "rows\n"
  )
  
  # Empty file = valid zero robust taxa result
  if (nrow(x) == 0) {
    next
  }
  
  # Standardise relevant types
  if ("taxon" %in% names(x)) {
    x$taxon <- as.character(x$taxon)
  }
  
  if ("Genus" %in% names(x)) {
    x$Genus <- as.character(x$Genus)
  }
  
  if ("Contrast" %in% names(x)) {
    x$Contrast <- as.character(x$Contrast)
  }
  
  if ("Higher" %in% names(x)) {
    x$Higher <- as.character(x$Higher)
  }
  
  x$Source_file <- basename(f)
  
  H1_list[[length(H1_list) + 1]] <- x
}


# ============================================================
# 3. COMBINE
# ============================================================

H1_robust <- bind_rows(
  H1_list
) %>%
  as_tibble()


cat("\n========================================\n")
cat("H1 ROBUST CONTRAST ROWS\n")
cat("========================================\n")

cat(
  "Total robust contrast rows:",
  nrow(H1_robust),
  "\n"
)


# ============================================================
# 4. UNIQUE ROBUST GENERA
# ============================================================

H1_genera <- H1_robust %>%
  filter(
    !is.na(Genus),
    trimws(Genus) != ""
  ) %>%
  distinct(
    Genus
  )


cat(
  "Unique robust genera:",
  nrow(H1_genera),
  "\n"
)


# ============================================================
# 5. MAP TO FULL TAXONOMY
#
# tax16 already exists from the earlier domain audit
# ============================================================

H1_domain <- tax16 %>%
  filter(
    Genus %in% H1_genera$Genus
  ) %>%
  distinct(
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus
  ) %>%
  arrange(
    Kingdom,
    Genus
  ) %>%
  as_tibble()


cat("\n========================================\n")
cat("H1 ROBUST GENERA — FULL TAXONOMY\n")
cat("========================================\n")

print(
  H1_domain,
  n = Inf,
  width = Inf
)


# ============================================================
# 6. KINGDOM COUNTS
# ============================================================

H1_kingdom_counts <- H1_domain %>%
  count(
    Kingdom,
    name = "Robust_genera"
  )


cat("\n========================================\n")
cat("H1 ROBUST GENERA — KINGDOM COUNTS\n")
cat("========================================\n")

print(
  H1_kingdom_counts,
  n = Inf,
  width = Inf
)


# ============================================================
# 7. SHOW ARCHAEAL H1 GENERA, IF ANY
# ============================================================

cat("\n========================================\n")
cat("H1 ROBUST ARCHAEAL GENERA\n")
cat("========================================\n")

print(
  H1_domain %>%
    filter(
      Kingdom == "Archaea"
    ),
  n = Inf,
  width = Inf
)


# ============================================================
# 8. SAVE AUDIT
# ============================================================

write.csv(
  H1_domain,
  file.path(
    audit_dir,
    "H1_robust_genus_domain.csv"
  ),
  row.names = FALSE
)

write.csv(
  H1_kingdom_counts,
  file.path(
    audit_dir,
    "H1_robust_genus_kingdom_counts.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("H1 DOMAIN CHECK COMPLETE\n")
cat("========================================\n")



# ============================================================
# 16S DOMAIN ABUNDANCE AUDIT
#
# PURPOSE:
# 1. Quantify Bacteria vs Archaea by ASV richness and reads
# 2. Do this for:
#       Full historical 16S dataset
#       H1 = bulk soil
#       H2 = rhizosphere
# 3. Quantify archaeal read abundance by ecological condition
#
# NO ECOLOGICAL ANALYSIS IS RERUN
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
})


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

audit_dir <- "outputs/16S_domain_audit"

dir.create(
  audit_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. LOAD HISTORICAL 16S OBJECT
# ============================================================

ps16 <- readRDS(
  "phyloseq_16S_corrected_metadata copy.rds"
)


# ============================================================
# 3. DEFINE H1 AND H2
# ============================================================

# H1 = bulk soil
ps_A <- subset_samples(
  ps16,
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


# H2 = rhizosphere
ps_B <- subset_samples(
  ps16,
  SoilPortion == "RZ" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE"
    )
)

ps_B <- prune_taxa(
  taxa_sums(ps_B) > 0,
  ps_B
)


# ============================================================
# 4. FUNCTION:
# BACTERIA / ARCHAEA — ASVs AND READS
# ============================================================

domain_summary <- function(ps, dataset_name) {
  
  counts <- as(
    otu_table(ps),
    "matrix"
  )
  
  if (taxa_are_rows(ps)) {
    counts <- t(counts)
  }
  
  counts <- counts[
    sample_names(ps),
    taxa_names(ps),
    drop = FALSE
  ]
  
  
  tax <- as.data.frame(
    tax_table(ps),
    stringsAsFactors = FALSE
  )
  
  tax <- tax[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  
  stopifnot(
    "Kingdom" %in% colnames(tax)
  )
  
  
  # Total reads per ASV
  taxon_reads <- colSums(counts)
  
  
  out <- tibble(
    Kingdom = tax$Kingdom,
    Reads = as.numeric(taxon_reads)
  ) %>%
    
    filter(
      Kingdom %in% c(
        "Bacteria",
        "Archaea"
      )
    ) %>%
    
    group_by(Kingdom) %>%
    
    summarise(
      ASVs = n(),
      Reads = sum(Reads),
      .groups = "drop"
    ) %>%
    
    mutate(
      Dataset = dataset_name,
      
      ASV_percent =
        100 * ASVs / sum(ASVs),
      
      Read_percent =
        100 * Reads / sum(Reads)
    ) %>%
    
    select(
      Dataset,
      Kingdom,
      ASVs,
      ASV_percent,
      Reads,
      Read_percent
    )
  
  out
}


# ============================================================
# 5. RUN DOMAIN SUMMARY
# ============================================================

domain_full <- domain_summary(
  ps16,
  "Full 16S dataset"
)

domain_A <- domain_summary(
  ps_A,
  "H1 - Bulk soil"
)

domain_B <- domain_summary(
  ps_B,
  "H2 - Rhizosphere"
)


domain_all <- bind_rows(
  domain_full,
  domain_A,
  domain_B
)


cat("\n========================================\n")
cat("BACTERIA / ARCHAEA — ASVs AND READS\n")
cat("========================================\n")

print(
  domain_all,
  n = Inf,
  width = Inf
)


write.csv(
  domain_all,
  file.path(
    audit_dir,
    "16S_Bacteria_Archaea_ASV_and_read_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. FUNCTION:
# SAMPLE-LEVEL ARCHAEAL READ PERCENTAGE
# ============================================================

sample_archaea <- function(ps, model_name) {
  
  counts <- as(
    otu_table(ps),
    "matrix"
  )
  
  if (taxa_are_rows(ps)) {
    counts <- t(counts)
  }
  
  counts <- counts[
    sample_names(ps),
    taxa_names(ps),
    drop = FALSE
  ]
  
  
  tax <- as.data.frame(
    tax_table(ps),
    stringsAsFactors = FALSE
  )
  
  tax <- tax[
    colnames(counts),
    ,
    drop = FALSE
  ]
  
  
  archaeal_taxa <- rownames(tax)[
    tax$Kingdom == "Archaea"
  ]
  
  bacterial_taxa <- rownames(tax)[
    tax$Kingdom == "Bacteria"
  ]
  
  
  archaea_reads <- rowSums(
    counts[
      ,
      archaeal_taxa,
      drop = FALSE
    ]
  )
  
  bacteria_reads <- rowSums(
    counts[
      ,
      bacterial_taxa,
      drop = FALSE
    ]
  )
  
  
  meta <- as.data.frame(
    sample_data(ps),
    stringsAsFactors = FALSE
  )
  
  meta <- meta[
    rownames(counts),
    ,
    drop = FALSE
  ]
  
  
  tibble(
    Model = model_name,
    SampleID = rownames(counts),
    Site = meta$Site,
    Condition = meta$Treatment,
    Archaea_reads = as.numeric(archaea_reads),
    Bacteria_reads = as.numeric(bacteria_reads)
  ) %>%
    
    mutate(
      Prokaryotic_reads =
        Archaea_reads + Bacteria_reads,
      
      Archaea_percent =
        100 *
        Archaea_reads /
        Prokaryotic_reads
    )
}


# ============================================================
# 7. CALCULATE SAMPLE-LEVEL VALUES
# ============================================================

archaea_A <- sample_archaea(
  ps_A,
  "H1 - Bulk soil"
)

archaea_B <- sample_archaea(
  ps_B,
  "H2 - Rhizosphere"
)


archaea_samples <- bind_rows(
  archaea_A,
  archaea_B
)


# ============================================================
# 8. ARCHAEAL ABUNDANCE BY CONDITION
#
# Reports BOTH:
#
# Mean_archaea_percent:
#   average sample-level archaeal percentage
#
# Pooled_archaea_percent:
#   archaeal reads / total prokaryotic reads across all
#   samples in that condition
# ============================================================

archaea_condition <- archaea_samples %>%
  
  group_by(
    Model,
    Condition
  ) %>%
  
  summarise(
    n = n(),
    
    Mean_archaea_percent =
      mean(
        Archaea_percent,
        na.rm = TRUE
      ),
    
    SD_archaea_percent =
      sd(
        Archaea_percent,
        na.rm = TRUE
      ),
    
    Median_archaea_percent =
      median(
        Archaea_percent,
        na.rm = TRUE
      ),
    
    Total_archaea_reads =
      sum(
        Archaea_reads,
        na.rm = TRUE
      ),
    
    Total_bacteria_reads =
      sum(
        Bacteria_reads,
        na.rm = TRUE
      ),
    
    Pooled_archaea_percent =
      100 *
      Total_archaea_reads /
      (
        Total_archaea_reads +
          Total_bacteria_reads
      ),
    
    .groups = "drop"
  )


cat("\n========================================\n")
cat("ARCHAEAL READ ABUNDANCE BY CONDITION\n")
cat("========================================\n")

print(
  archaea_condition,
  n = Inf,
  width = Inf
)


# ============================================================
# 9. SAVE
# ============================================================

write.csv(
  archaea_samples,
  file.path(
    audit_dir,
    "16S_sample_level_archaeal_read_percent.csv"
  ),
  row.names = FALSE
)

write.csv(
  archaea_condition,
  file.path(
    audit_dir,
    "16S_archaeal_read_abundance_by_condition.csv"
  ),
  row.names = FALSE
)


cat("\n========================================\n")
cat("DOMAIN ABUNDANCE AUDIT COMPLETE\n")
cat("========================================\n")

cat(
  "Outputs saved in:\n",
  normalizePath(audit_dir),
  "\n"
)