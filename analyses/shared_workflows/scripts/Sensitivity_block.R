# ============================================================
# RESTORATION-ONLY BLOCK SENSITIVITY
# H1 — BULK SOIL
#
# Treatments:
#   BESE, BARE, DEGRADED
#
# Site = fixed effect
# SiteBlock = Site x Block
#
# Alpha / distance:
#   response ~ Site + Treatment + (1 | SiteBlock)
#
# PERMANOVA:
#   Bray ~ Site + Treatment
#   permutations restricted within SiteBlock
#
# PERMDISP:
#   treatment dispersions
#   permutations restricted within SiteBlock
#
# Natural is NOT included in the treatment model.
# Natural samples are used only as references for
# distance-to-Natural.
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(vegan)
  library(permute)
})

set.seed(123)

outdir <- "outputs/restoration_block_sensitivity/H1"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

restoration_levels <- c(
  "BESE",
  "BARE",
  "DEGRADED"
)

site_levels <- c(
  "ACP",
  "BASS",
  "TRUG"
)


# ============================================================
# 1. PREPARE H1 RESTORATION-ONLY DATA
# ============================================================

prepare_H1_blocked <- function(ps, kingdom) {
  
  ps_sub <- subset_samples(
    ps,
    SoilPortion == "BS" &
      Treatment %in% restoration_levels
  )
  
  ps_sub <- prune_taxa(
    taxa_sums(ps_sub) > 0,
    ps_sub
  )
  
  md <- as.data.frame(
    sample_data(ps_sub),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"Block" %in% names(md)) {
    stop("Block column not found.")
  }
  
  if (anyNA(md$Block)) {
    stop(
      kingdom,
      ": missing Block values in restoration samples."
    )
  }
  
  md$Site <- factor(
    md$Site,
    levels = site_levels
  )
  
  md$Treatment <- factor(
    md$Treatment,
    levels = restoration_levels
  )
  
  # IMPORTANT:
  # Block 1 at ACP is NOT Block 1 at BASS etc.
  md$SiteBlock <- interaction(
    md$Site,
    md$Block,
    drop = TRUE,
    sep = "_Block"
  )
  
  sample_data(ps_sub) <- sample_data(md)
  
  cat("\n\n========================================\n")
  cat(
    kingdom,
    "— H1 BLOCK DESIGN\n"
  )
  cat("========================================\n")
  
  cat(
    "Samples:",
    nsamples(ps_sub),
    "\n"
  )
  
  cat(
    "Site-specific blocks:",
    nlevels(md$SiteBlock),
    "\n\n"
  )
  
  print(
    table(
      md$Site,
      md$Treatment
    )
  )
  
  cat("\nSAMPLES PER SITE-BLOCK x TREATMENT\n")
  
  block_table <- table(
    md$SiteBlock,
    md$Treatment
  )
  
  print(
    block_table
  )
  
  complete_blocks <- rowSums(
    block_table > 0
  ) == 3
  
  cat(
    "\nBlocks containing all 3 treatments:",
    sum(complete_blocks),
    "of",
    nrow(block_table),
    "\n"
  )
  
  list(
    ps = ps_sub,
    meta = md
  )
}


# ============================================================
# 2. ALPHA DIVERSITY
# ============================================================

run_blocked_alpha <- function(dat, kingdom) {
  
  ps <- dat$ps
  md <- dat$meta
  
  rich <- estimate_richness(
    ps,
    measures = c(
      "Chao1",
      "Shannon"
    )
  )
  
  rich$SampleID_raw <- rownames(rich)
  
  rich$SampleID <- sub(
    "^X(?=\\d)",
    "",
    rich$SampleID_raw,
    perl = TRUE
  )
  
  md2 <- md
  
  md2$SampleID_raw <- rownames(md2)
  
  md2$SampleID <- sub(
    "^X(?=\\d)",
    "",
    md2$SampleID_raw,
    perl = TRUE
  )
  
  alpha <- merge(
    rich[
      ,
      c(
        "SampleID",
        "Chao1",
        "Shannon"
      )
    ],
    md2[
      ,
      c(
        "SampleID",
        "Site",
        "Treatment",
        "SiteBlock"
      )
    ],
    by = "SampleID",
    all.x = TRUE,
    sort = FALSE
  )
  
  if (anyNA(alpha$Site)) {
    stop(
      kingdom,
      ": alpha metadata join failed."
    )
  }
  
  
  # ----------------------------------------------------------
  # Chao1
  # ----------------------------------------------------------
  
  fit_chao <- lmerTest::lmer(
    Chao1 ~
      Site +
      Treatment +
      (1 | SiteBlock),
    data = alpha,
    REML = TRUE
  )
  
  test_chao <- as.data.frame(
    lmerTest::anova(
      fit_chao,
      ddf = "Satterthwaite"
    )
  ) %>%
    rownames_to_column(
      "Term"
    )
  
  emm_chao <- emmeans(
    fit_chao,
    ~ Treatment
  )
  
  pair_chao <- as.data.frame(
    pairs(
      emm_chao,
      adjust = "holm"
    )
  )
  
  
  # ----------------------------------------------------------
  # Shannon
  # ----------------------------------------------------------
  
  fit_shannon <- lmerTest::lmer(
    Shannon ~
      Site +
      Treatment +
      (1 | SiteBlock),
    data = alpha,
    REML = TRUE
  )
  
  test_shannon <- as.data.frame(
    lmerTest::anova(
      fit_shannon,
      ddf = "Satterthwaite"
    )
  ) %>%
    rownames_to_column(
      "Term"
    )
  
  emm_shannon <- emmeans(
    fit_shannon,
    ~ Treatment
  )
  
  pair_shannon <- as.data.frame(
    pairs(
      emm_shannon,
      adjust = "holm"
    )
  )
  
  
  cat("\n========================================\n")
  cat(
    kingdom,
    "— BLOCKED ALPHA\n"
  )
  cat("========================================\n")
  
  cat("\nCHAO1 — TREATMENT EFFECT\n")
  
  print(
    test_chao[
      test_chao$Term == "Treatment",
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  cat("\nCHAO1 — HOLM PAIRWISE\n")
  
  print(
    pair_chao,
    row.names = FALSE
  )
  
  cat("\nSHANNON — TREATMENT EFFECT\n")
  
  print(
    test_shannon[
      test_shannon$Term == "Treatment",
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  cat("\nSHANNON — HOLM PAIRWISE\n")
  
  print(
    pair_shannon,
    row.names = FALSE
  )
  
  cat(
    "\nChao1 singular fit:",
    lme4::isSingular(
      fit_chao
    ),
    "\n"
  )
  
  cat(
    "Shannon singular fit:",
    lme4::isSingular(
      fit_shannon
    ),
    "\n"
  )
  
  
  write.csv(
    test_chao,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_Chao1_ANOVA.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    pair_chao,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_Chao1_pairwise.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    test_shannon,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_Shannon_ANOVA.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    pair_shannon,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_Shannon_pairwise.csv"
      )
    ),
    row.names = FALSE
  )
  
  list(
    alpha = alpha,
    Chao1_model = fit_chao,
    Shannon_model = fit_shannon
  )
}


# ============================================================
# 3. BETA DIVERSITY
# ============================================================

run_blocked_beta <- function(dat, kingdom) {
  
  ps <- dat$ps
  
  ps_rel <- transform_sample_counts(
    ps,
    function(x) {
      x / sum(x)
    }
  )
  
  comm <- as(
    otu_table(ps_rel),
    "matrix"
  )
  
  if (taxa_are_rows(ps_rel)) {
    comm <- t(comm)
  }
  
  md <- as.data.frame(
    sample_data(ps_rel),
    stringsAsFactors = FALSE
  )
  
  md <- md[
    rownames(comm),
    ,
    drop = FALSE
  ]
  
  bray <- vegan::vegdist(
    comm,
    method = "bray"
  )
  
  stopifnot(
    identical(
      attr(
        bray,
        "Labels"
      ),
      rownames(md)
    )
  )
  
  
  # ----------------------------------------------------------
  # Restricted permutations within Site x Block
  # ----------------------------------------------------------
  
  perm_block <- permute::how(
    nperm = 9999,
    blocks = md$SiteBlock
  )
  
  
  # ----------------------------------------------------------
  # PERMANOVA
  # ----------------------------------------------------------
  
  set.seed(123)
  
  permanova <- vegan::adonis2(
    bray ~
      Site +
      Treatment,
    data = md,
    permutations = perm_block,
    by = "margin"
  )
  
  
  # ----------------------------------------------------------
  # PERMDISP
  # ----------------------------------------------------------
  
  dispersion <- vegan::betadisper(
    bray,
    group = md$Treatment,
    type = "median",
    bias.adjust = TRUE
  )
  
  set.seed(123)
  
  permdisp <- vegan::permutest(
    dispersion,
    permutations = perm_block
  )
  
  
  # ----------------------------------------------------------
  # Pairwise blocked PERMANOVA
  # ----------------------------------------------------------
  
  comparisons <- list(
    c(
      "BESE",
      "BARE"
    ),
    c(
      "BESE",
      "DEGRADED"
    ),
    c(
      "BARE",
      "DEGRADED"
    )
  )
  
  pair_perm <- lapply(
    comparisons,
    function(groups) {
      
      ids <- rownames(md)[
        md$Treatment %in% groups
      ]
      
      md2 <- droplevels(
        md[
          ids,
          ,
          drop = FALSE
        ]
      )
      
      d2 <- as.dist(
        as.matrix(bray)[
          ids,
          ids
        ]
      )
      
      perm2 <- permute::how(
        nperm = 9999,
        blocks = md2$SiteBlock
      )
      
      set.seed(123)
      
      fit <- vegan::adonis2(
        d2 ~
          Site +
          Treatment,
        data = md2,
        permutations = perm2,
        by = "margin"
      )
      
      data.frame(
        Contrast =
          paste(
            groups,
            collapse = " vs "
          ),
        
        n =
          nrow(md2),
        
        R2 =
          fit[
            "Treatment",
            "R2"
          ],
        
        pseudo_F =
          fit[
            "Treatment",
            "F"
          ],
        
        p_value =
          fit[
            "Treatment",
            "Pr(>F)"
          ]
      )
    }
  ) %>%
    bind_rows() %>%
    mutate(
      p_Holm =
        p.adjust(
          p_value,
          method = "holm"
        )
    )
  
  
  # ----------------------------------------------------------
  # Pairwise blocked PERMDISP
  # ----------------------------------------------------------
  
  pair_disp <- lapply(
    comparisons,
    function(groups) {
      
      ids <- rownames(md)[
        md$Treatment %in% groups
      ]
      
      md2 <- droplevels(
        md[
          ids,
          ,
          drop = FALSE
        ]
      )
      
      d2 <- as.dist(
        as.matrix(bray)[
          ids,
          ids
        ]
      )
      
      disp2 <- vegan::betadisper(
        d2,
        group =
          md2$Treatment,
        type = "median",
        bias.adjust = TRUE
      )
      
      perm2 <- permute::how(
        nperm = 9999,
        blocks = md2$SiteBlock
      )
      
      set.seed(123)
      
      test2 <- vegan::permutest(
        disp2,
        permutations = perm2
      )
      
      data.frame(
        Contrast =
          paste(
            groups,
            collapse = " vs "
          ),
        
        F =
          test2$tab[
            1,
            "F"
          ],
        
        p_value =
          test2$tab[
            1,
            "Pr(>F)"
          ]
      )
    }
  ) %>%
    bind_rows() %>%
    mutate(
      p_Holm =
        p.adjust(
          p_value,
          method = "holm"
        )
    )
  
  
  cat("\n========================================\n")
  cat(
    kingdom,
    "— BLOCKED BETA\n"
  )
  cat("========================================\n")
  
  cat("\nPERMANOVA — TREATMENT\n")
  
  print(
    permanova[
      "Treatment",
      ,
      drop = FALSE
    ]
  )
  
  cat("\nPAIRWISE PERMANOVA\n")
  
  print(
    pair_perm,
    row.names = FALSE
  )
  
  cat("\nPERMDISP\n")
  
  print(
    permdisp$tab,
    row.names = TRUE
  )
  
  cat("\nPAIRWISE PERMDISP\n")
  
  print(
    pair_disp,
    row.names = FALSE
  )
  
  
  write.csv(
    as.data.frame(
      permanova
    ),
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_PERMANOVA.csv"
      )
    )
  )
  
  write.csv(
    pair_perm,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_pairwise_PERMANOVA.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    as.data.frame(
      permdisp$tab
    ),
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_PERMDISP.csv"
      )
    )
  )
  
  write.csv(
    pair_disp,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_pairwise_PERMDISP.csv"
      )
    ),
    row.names = FALSE
  )
  
  list(
    bray = bray,
    meta = md,
    permanova = permanova,
    permdisp = permdisp,
    pairwise_permanova = pair_perm,
    pairwise_permdisp = pair_disp
  )
}


# ============================================================
# 4. DISTANCE TO SITE-SPECIFIC NATURAL REFERENCE
#
# Natural is used ONLY to construct the reference distance.
# Statistical model includes only BESE / BARE / DEGRADED.
# ============================================================

run_blocked_distance_to_natural <- function(
    ps_original,
    kingdom
) {
  
  ps_full <- subset_samples(
    ps_original,
    SoilPortion == "BS" &
      Treatment %in% c(
        "NATURAL",
        "BESE",
        "BARE",
        "DEGRADED"
      )
  )
  
  ps_full <- prune_taxa(
    taxa_sums(ps_full) > 0,
    ps_full
  )
  
  ps_rel <- transform_sample_counts(
    ps_full,
    function(x) {
      x / sum(x)
    }
  )
  
  comm <- as(
    otu_table(ps_rel),
    "matrix"
  )
  
  if (taxa_are_rows(ps_rel)) {
    comm <- t(comm)
  }
  
  md <- as.data.frame(
    sample_data(ps_rel),
    stringsAsFactors = FALSE
  )
  
  md <- md[
    rownames(comm),
    ,
    drop = FALSE
  ]
  
  bray <- as.matrix(
    vegan::vegdist(
      comm,
      method = "bray"
    )
  )
  
  
  # ----------------------------------------------------------
  # Mean distance of each restoration sample to
  # same-Site Natural reference samples
  # ----------------------------------------------------------
  
  restoration_ids <- rownames(md)[
    md$Treatment %in%
      restoration_levels
  ]
  
  distance_df <- lapply(
    restoration_ids,
    function(id) {
      
      refs <- rownames(md)[
        md$Site == md[
          id,
          "Site"
        ] &
          md$Treatment ==
          "NATURAL"
      ]
      
      if (length(refs) == 0) {
        stop(
          "No Natural reference available for ",
          id
        )
      }
      
      data.frame(
        SampleID = id,
        
        Site =
          as.character(
            md[
              id,
              "Site"
            ]
          ),
        
        Treatment =
          as.character(
            md[
              id,
              "Treatment"
            ]
          ),
        
        Block =
          as.character(
            md[
              id,
              "Block"
            ]
          ),
        
        Distance_to_Natural =
          mean(
            bray[
              id,
              refs
            ]
          ),
        
        stringsAsFactors = FALSE
      )
    }
  ) %>%
    bind_rows()
  
  
  distance_df$Site <- factor(
    distance_df$Site,
    levels = site_levels
  )
  
  distance_df$Treatment <- factor(
    distance_df$Treatment,
    levels = restoration_levels
  )
  
  distance_df$SiteBlock <- interaction(
    distance_df$Site,
    distance_df$Block,
    drop = TRUE,
    sep = "_Block"
  )
  
  
  # ----------------------------------------------------------
  # Blocked mixed model
  # ----------------------------------------------------------
  
  fit <- lmerTest::lmer(
    Distance_to_Natural ~
      Site +
      Treatment +
      (1 | SiteBlock),
    data = distance_df,
    REML = TRUE
  )
  
  test <- as.data.frame(
    lmerTest::anova(
      fit,
      ddf = "Satterthwaite"
    )
  ) %>%
    rownames_to_column(
      "Term"
    )
  
  emm <- emmeans(
    fit,
    ~ Treatment
  )
  
  pairs_out <- as.data.frame(
    pairs(
      emm,
      adjust = "holm"
    )
  )
  
  
  cat("\n========================================\n")
  cat(
    kingdom,
    "— BLOCKED DISTANCE TO NATURAL\n"
  )
  cat("========================================\n")
  
  print(
    test[
      test$Term == "Treatment",
      ,
      drop = FALSE
    ],
    row.names = FALSE
  )
  
  cat("\nHOLM PAIRWISE\n")
  
  print(
    pairs_out,
    row.names = FALSE
  )
  
  cat(
    "\nSingular fit:",
    lme4::isSingular(
      fit
    ),
    "\n"
  )
  
  
  write.csv(
    distance_df,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_distance_to_Natural_raw.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    test,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_distance_to_Natural_ANOVA.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    pairs_out,
    file.path(
      outdir,
      paste0(
        kingdom,
        "_H1_blocked_distance_to_Natural_pairwise.csv"
      )
    ),
    row.names = FALSE
  )
  
  list(
    data = distance_df,
    model = fit,
    test = test,
    pairs = pairs_out
  )
}


# ============================================================
# 5. RUN H1 — PROKARYOTES
# ============================================================

H1_block_prok <- prepare_H1_blocked(
  ps_bacteria,
  "Prokaryotes"
)

H1_block_prok_alpha <- run_blocked_alpha(
  H1_block_prok,
  "Prokaryotes"
)

H1_block_prok_beta <- run_blocked_beta(
  H1_block_prok,
  "Prokaryotes"
)

H1_block_prok_distance <-
  run_blocked_distance_to_natural(
    ps_bacteria,
    "Prokaryotes"
  )


# ============================================================
# 6. RUN H1 — FUNGI
# ============================================================

H1_block_fungi <- prepare_H1_blocked(
  ps_fungi,
  "Fungi"
)

H1_block_fungi_alpha <- run_blocked_alpha(
  H1_block_fungi,
  "Fungi"
)

H1_block_fungi_beta <- run_blocked_beta(
  H1_block_fungi,
  "Fungi"
)

H1_block_fungi_distance <-
  run_blocked_distance_to_natural(
    ps_fungi,
    "Fungi"
  )


# ============================================================
# 7. SAVE MODELS
# ============================================================

saveRDS(
  list(
    prokaryotes = list(
      alpha = H1_block_prok_alpha,
      beta = H1_block_prok_beta,
      distance = H1_block_prok_distance
    ),
    
    fungi = list(
      alpha = H1_block_fungi_alpha,
      beta = H1_block_fungi_beta,
      distance = H1_block_fungi_distance
    ),
    
    sessionInfo = sessionInfo()
  ),
  
  file.path(
    outdir,
    "H1_restoration_only_block_sensitivity.rds"
  )
)

cat("\n\n========================================\n")
cat("H1 BLOCK SENSITIVITY COMPLETE\n")
cat("========================================\n")

# ============================================================
# H1 — RESTORATION-ONLY BLOCKED ANCOM-BC2
# GENUS LEVEL
#
# Samples:
#   BESE + BARE + DEGRADED
#   Bulk soil only
#
# Site      = fixed effect
# Treatment = fixed effect
# SiteBlock = random intercept
#
# Model:
#   Site + Treatment + (1 | SiteBlock)
#
# Pairwise contrasts:
#   BESE vs BARE
#   BESE vs DEGRADED
#   BARE vs DEGRADED
#
# Robust calls:
#   diff_robust == TRUE
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

set.seed(123)

outdir_ancom <-
  "outputs/restoration_block_sensitivity/H1/ANCOMBC2"

dir.create(
  outdir_ancom,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. PREPARE GENUS-LEVEL INPUT
# ============================================================

prepare_blocked_genus <- function(
    block_object,
    kingdom_label
) {
  
  ps <- block_object$ps
  
  # ----------------------------------------------------------
  # Check metadata
  # ----------------------------------------------------------
  
  md <- as.data.frame(
    sample_data(ps),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  required <- c(
    "Site",
    "Treatment",
    "SiteBlock"
  )
  
  missing_vars <- setdiff(
    required,
    names(md)
  )
  
  if (length(missing_vars) > 0) {
    stop(
      kingdom_label,
      ": missing metadata variables: ",
      paste(
        missing_vars,
        collapse = ", "
      )
    )
  }
  
  md$Site <- factor(
    as.character(md$Site),
    levels = c(
      "ACP",
      "BASS",
      "TRUG"
    )
  )
  
  # BESE is reference level
  md$Treatment <- factor(
    as.character(md$Treatment),
    levels = c(
      "BESE",
      "BARE",
      "DEGRADED"
    )
  )
  
  md$SiteBlock <- factor(
    as.character(md$SiteBlock)
  )
  
  sample_data(ps) <-
    sample_data(md)
  
  
  # ----------------------------------------------------------
  # Remove unresolved genus assignments
  # ----------------------------------------------------------
  
  tax <- as.data.frame(
    tax_table(ps),
    stringsAsFactors = FALSE
  )
  
  if (!"Genus" %in% colnames(tax)) {
    stop(
      kingdom_label,
      ": taxonomy has no Genus column."
    )
  }
  
  genus <- trimws(
    as.character(
      tax$Genus
    )
  )
  
  keep_named <- !(
    is.na(genus) |
      genus == "" |
      genus %in% c(
        "NA",
        "Unclassified",
        "unclassified"
      )
  )
  
  ps <- prune_taxa(
    keep_named,
    ps
  )
  
  
  # ----------------------------------------------------------
  # Agglomerate to Genus
  # ----------------------------------------------------------
  
  ps <- tax_glom(
    ps,
    taxrank = "Genus",
    NArm = TRUE
  )
  
  genera_before <- ntaxa(ps)
  
  
  # ----------------------------------------------------------
  # 10% prevalence filter
  # restoration-only sample set
  # ----------------------------------------------------------
  
  otu <- as(
    otu_table(ps),
    "matrix"
  )
  
  if (!taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  
  prevalence <- rowSums(
    otu > 0
  )
  
  prevalence_threshold <- ceiling(
    0.10 *
      nsamples(ps)
  )
  
  keep_prev <- prevalence >=
    prevalence_threshold
  
  ps <- prune_taxa(
    keep_prev,
    ps
  )
  
  ps <- prune_samples(
    sample_sums(ps) > 0,
    ps
  )
  
  
  # ----------------------------------------------------------
  # Print input summary
  # ----------------------------------------------------------
  
  cat("\n========================================\n")
  cat(
    kingdom_label,
    "— BLOCKED GENUS INPUT\n"
  )
  cat("========================================\n")
  
  cat(
    "Samples:",
    nsamples(ps),
    "\n"
  )
  
  cat(
    "SiteBlocks:",
    nlevels(
      factor(
        sample_data(ps)$SiteBlock
      )
    ),
    "\n"
  )
  
  cat(
    "Resolved genera before prevalence filter:",
    genera_before,
    "\n"
  )
  
  cat(
    "10% prevalence threshold:",
    prevalence_threshold,
    "samples\n"
  )
  
  cat(
    "Genera entering ANCOM-BC2:",
    ntaxa(ps),
    "\n"
  )
  
  cat("\nTreatment counts\n")
  
  print(
    table(
      sample_data(ps)$Treatment
    )
  )
  
  list(
    ps = ps,
    prevalence_threshold =
      prevalence_threshold,
    genera_before =
      genera_before
  )
}


# ============================================================
# 2. RUN BLOCKED ANCOM-BC2
# ============================================================

run_blocked_ancom <- function(
    prepared,
    kingdom_label
) {
  
  ps <- prepared$ps
  
  cat("\n\n========================================\n")
  cat(
    "RUNNING BLOCKED ANCOM-BC2 —",
    kingdom_label,
    "\n"
  )
  cat("========================================\n")
  
  
  # ----------------------------------------------------------
  # ANCOM-BC2
  # ----------------------------------------------------------
  
  set.seed(123)
  
  out <- ANCOMBC::ancombc2(
    
    data = ps,
    
    tax_level = NULL,
    
    fix_formula =
      "Site + Treatment",
    
    rand_formula =
      "(1 | SiteBlock)",
    
    p_adj_method =
      "BH",
    
    pseudo_sens =
      TRUE,
    
    prv_cut =
      0,
    
    lib_cut =
      0,
    
    group =
      "Treatment",
    
    struc_zero =
      TRUE,
    
    neg_lb =
      FALSE,
    
    alpha =
      0.05,
    
    global =
      TRUE,
    
    pairwise =
      TRUE,
    
    dunnet =
      FALSE,
    
    trend =
      FALSE,
    
    n_cl =
      1,
    
    verbose =
      TRUE
  )
  
  
  # ==========================================================
  # 3. BASIC OUTPUT ARCHITECTURE
  # ==========================================================
  
  cat("\n========================================\n")
  cat(
    kingdom_label,
    "— ANCOM-BC2 OUTPUT\n"
  )
  cat("========================================\n")
  
  cat(
    "Output components:\n"
  )
  
  print(
    names(out)
  )
  
  
  # ==========================================================
  # 4. STRUCTURAL ZEROS
  # ==========================================================
  
  zero_tab <- as.data.frame(
    out$zero_ind,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"taxon" %in% names(zero_tab)) {
    zero_tab$taxon <-
      rownames(zero_tab)
  }
  
  zero_cols <- setdiff(
    names(zero_tab),
    "taxon"
  )
  
  if (length(zero_cols) > 0) {
    
    zero_tab$Structural_zero_any <-
      apply(
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
    
    zero_tab$Structural_zero_any <-
      FALSE
  }
  
  n_structural <- sum(
    zero_tab$Structural_zero_any,
    na.rm = TRUE
  )
  
  
  # ==========================================================
  # 5. GLOBAL RESULT
  # ==========================================================
  
  global <- as.data.frame(
    out$res_global,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"taxon" %in% names(global)) {
    global$taxon <-
      rownames(global)
  }
  
  if ("diff_robust_abn" %in% names(global)) {
    
    robust_global <- global[
      !is.na(
        global$diff_robust_abn
      ) &
        global$diff_robust_abn,
      ,
      drop = FALSE
    ]
    
  } else {
    
    robust_global <-
      global[0, , drop = FALSE]
  }
  
  
  # ==========================================================
  # 6. PAIRWISE RESULT
  # ==========================================================
  
  pair <- as.data.frame(
    out$res_pair,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (!"taxon" %in% names(pair)) {
    pair$taxon <-
      rownames(pair)
  }
  
  cat("\nPAIRWISE COLUMN NAMES\n")
  
  print(
    names(pair)
  )
  
  
  # ==========================================================
  # 7. TAXONOMY LOOKUP
  # ==========================================================
  
  tax <- as.data.frame(
    tax_table(ps),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  tax$taxon <-
    rownames(tax)
  
  tax_keep <- intersect(
    c(
      "taxon",
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus"
    ),
    names(tax)
  )
  
  tax <- tax[
    ,
    tax_keep,
    drop = FALSE
  ]
  
  
  # ==========================================================
  # 8. EXTRACT ROBUST PAIRWISE CONTRASTS DYNAMICALLY
  #
  # Avoid assuming exact ANCOM-BC2 column suffixes.
  # ==========================================================
  
  robust_cols <- grep(
    "^diff_robust_",
    names(pair),
    value = TRUE
  )
  
  if (length(robust_cols) == 0) {
    
    stop(
      kingdom_label,
      ": no diff_robust pairwise columns found."
    )
  }
  
  
  robust_list <- list()
  
  counter <- 1
  
  
  for (robust_col in robust_cols) {
    
    suffix <- sub(
      "^diff_robust_",
      "",
      robust_col
    )
    
    lfc_col <- paste0(
      "lfc_",
      suffix
    )
    
    se_col <- paste0(
      "se_",
      suffix
    )
    
    p_col <- paste0(
      "p_",
      suffix
    )
    
    q_col <- paste0(
      "q_",
      suffix
    )
    
    
    if (!lfc_col %in% names(pair)) {
      next
    }
    
    
    keep <- !is.na(
      pair[[robust_col]]
    ) &
      as.logical(
        pair[[robust_col]]
      )
    
    
    if (!any(keep)) {
      next
    }
    
    
    tmp <- data.frame(
      
      taxon =
        pair$taxon[keep],
      
      Contrast_code =
        suffix,
      
      LFC =
        pair[[lfc_col]][keep],
      
      SE =
        if (
          se_col %in%
          names(pair)
        ) {
          pair[[se_col]][keep]
        } else {
          NA_real_
        },
      
      p =
        if (
          p_col %in%
          names(pair)
        ) {
          pair[[p_col]][keep]
        } else {
          NA_real_
        },
      
      q_BH =
        if (
          q_col %in%
          names(pair)
        ) {
          pair[[q_col]][keep]
        } else {
          NA_real_
        },
      
      stringsAsFactors = FALSE
    )
    
    
    robust_list[[counter]] <-
      tmp
    
    counter <- counter + 1
  }
  
  
  if (length(robust_list) > 0) {
    
    robust_pair <- bind_rows(
      robust_list
    )
    
    robust_pair <- left_join(
      robust_pair,
      tax,
      by = "taxon"
    )
    
    robust_pair <- robust_pair %>%
      arrange(
        Contrast_code,
        desc(
          abs(LFC)
        )
      )
    
  } else {
    
    robust_pair <- data.frame(
      taxon = character(0),
      Contrast_code = character(0),
      LFC = numeric(0),
      SE = numeric(0),
      p = numeric(0),
      q_BH = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  
  # ==========================================================
  # 9. ROBUST PAIRWISE COUNTS
  # ==========================================================
  
  pair_counts <- data.frame(
    Contrast_code =
      robust_cols,
    Robust_taxa =
      sapply(
        robust_cols,
        function(x) {
          sum(
            as.logical(
              pair[[x]]
            ),
            na.rm = TRUE
          )
        }
      ),
    stringsAsFactors = FALSE
  )
  
  pair_counts$Contrast_code <-
    sub(
      "^diff_robust_",
      "",
      pair_counts$Contrast_code
    )
  
  
  # ==========================================================
  # 10. PRINT FINAL SUMMARY
  # ==========================================================
  
  cat("\n========================================\n")
  cat(
    kingdom_label,
    "— BLOCKED ANCOM-BC2 SUMMARY\n"
  )
  cat("========================================\n")
  
  cat(
    "Genera entering ANCOM-BC2:",
    ntaxa(ps),
    "\n"
  )
  
  cat(
    "Structural-zero genera:",
    n_structural,
    "\n"
  )
  
  cat(
    "Robust global genera:",
    nrow(
      robust_global
    ),
    "\n"
  )
  
  cat("\nROBUST PAIRWISE COUNTS\n")
  
  print(
    pair_counts,
    row.names = FALSE
  )
  
  
  cat("\nROBUST PAIRWISE GENERA\n")
  
  if (nrow(robust_pair) == 0) {
    
    cat(
      "No robust pairwise genera detected.\n"
    )
    
  } else {
    
    print(
      robust_pair,
      row.names = FALSE
    )
  }
  
  
  # ==========================================================
  # 11. SAVE
  # ==========================================================
  
  kingdom_file <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    kingdom_label
  )
  
  write.csv(
    global,
    file.path(
      outdir_ancom,
      paste0(
        kingdom_file,
        "_blocked_ANCOMBC2_global.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    pair,
    file.path(
      outdir_ancom,
      paste0(
        kingdom_file,
        "_blocked_ANCOMBC2_pairwise_all.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    robust_pair,
    file.path(
      outdir_ancom,
      paste0(
        kingdom_file,
        "_blocked_ANCOMBC2_pairwise_robust.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    zero_tab,
    file.path(
      outdir_ancom,
      paste0(
        kingdom_file,
        "_blocked_ANCOMBC2_structural_zeros.csv"
      )
    ),
    row.names = FALSE
  )
  
  
  list(
    ps =
      ps,
    
    output =
      out,
    
    global =
      global,
    
    robust_global =
      robust_global,
    
    pairwise =
      pair,
    
    pairwise_counts =
      pair_counts,
    
    robust_pairwise =
      robust_pair,
    
    structural_zeros =
      zero_tab
  )
}


# ============================================================
# 12. PREPARE BOTH DATASETS
# ============================================================

blocked_genus_prok <-
  prepare_blocked_genus(
    H1_block_prok,
    "Prokaryotes"
  )

blocked_genus_fungi <-
  prepare_blocked_genus(
    H1_block_fungi,
    "Fungi"
  )


# ============================================================
# 13. RUN PROKARYOTES
# ============================================================

H1_blocked_ANCOM_prok <-
  run_blocked_ancom(
    blocked_genus_prok,
    "Prokaryotes"
  )


# ============================================================
# 14. RUN FUNGI
# ============================================================

H1_blocked_ANCOM_fungi <-
  run_blocked_ancom(
    blocked_genus_fungi,
    "Fungi"
  )


# ============================================================
# 15. SAVE COMPLETE OBJECTS
# ============================================================

saveRDS(
  list(
    Prokaryotes =
      H1_blocked_ANCOM_prok,
    
    Fungi =
      H1_blocked_ANCOM_fungi,
    
    ANCOMBC_version =
      as.character(
        packageVersion(
          "ANCOMBC"
        )
      ),
    
    sessionInfo =
      sessionInfo()
  ),
  
  file.path(
    outdir_ancom,
    "H1_restoration_only_blocked_genus_ANCOMBC2.rds"
  )
)


cat("\n========================================\n")
cat("H1 BLOCKED ANCOM-BC2 COMPLETE\n")
cat("========================================\n")
