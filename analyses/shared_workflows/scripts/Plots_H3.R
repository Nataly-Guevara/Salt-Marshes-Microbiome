# ============================================================
# H3 — FUNGI — CHAO1 FIGURE
# Uses VERIFIED paired dataset: alpha_fungi_c
#
# Bulk soil   = circle
# Rhizosphere = triangle
# Lines connect paired samples
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Prepare plotting data
# ------------------------------------------------------------

plot_fun_chao <- alpha_fungi_c %>%
  dplyr::mutate(
    Ecological_condition = factor(
      Ecological_condition,
      levels = c("Natural", "BESE", "BARE")
    ),
    Soil_habitat = factor(
      Soil_habitat,
      levels = c("Bulk soil", "Rhizosphere")
    )
  )

# Verify
cat("\nSamples plotted:", nrow(plot_fun_chao), "\n")
cat("Paired plots:", dplyr::n_distinct(plot_fun_chao$PlotID), "\n")

print(
  plot_fun_chao %>%
    dplyr::distinct(PlotID, Ecological_condition) %>%
    dplyr::count(Ecological_condition)
)


# ------------------------------------------------------------
# 2. Colours
# ------------------------------------------------------------

treatment_cols <- c(
  "Natural" = "#9BC1BC",
  "BESE"    = "#FF5A5F",
  "BARE"    = "#FFA630"
)

# Filled shapes with black borders:
# 21 = circle
# 24 = triangle
habitat_shapes <- c(
  "Bulk soil"   = 21,
  "Rhizosphere" = 24
)


# ------------------------------------------------------------
# 3. Figure
# ------------------------------------------------------------

p_fun_chao_C <- ggplot(
  plot_fun_chao,
  aes(
    x = Soil_habitat,
    y = Chao1,
    group = PlotID
  )
) +
  
  # connect paired BS-RZ samples
  geom_line(
    colour = "grey65",
    linewidth = 0.55,
    alpha = 0.65
  ) +
  
  # individual samples
  geom_point(
    aes(
      shape = Soil_habitat,
      fill = Ecological_condition
    ),
    colour = "black",
    size = 3.4,
    stroke = 0.8,
    alpha = 0.95
  ) +
  
  facet_wrap(
    ~ Ecological_condition,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = treatment_cols,
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = habitat_shapes,
    name = "Soil habitat"
  ) +
  
  scale_y_log10() +
  
  labs(
    x = NULL,
    y = expression(Chao[1]~richness),
    shape = "Soil habitat"
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
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 11
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    )
  )


# ------------------------------------------------------------
# 4. Show figure
# ------------------------------------------------------------

print(p_fun_chao_C)


# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------

dir.create(
  "outputs/H3_alpha/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "outputs/H3_alpha/figures/H3_fungi_Chao1_paired.png",
  p_fun_chao_C,
  width = 8,
  height = 4.3,
  dpi = 600
)

ggsave(
  "outputs/H3_alpha/figures/H3_fungi_Chao1_paired.pdf",
  p_fun_chao_C,
  width = 8,
  height = 4.3
)

# ============================================================
# H3 — FUNGI — CHAO1
# Paired rhizosphere response plot
#
# Response = ln(RZ Chao1 / BS Chao1)
#
# > 0  = higher richness in rhizosphere
# = 0  = no difference
# < 0  = higher richness in bulk soil
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(emmeans)
})

# ============================================================
# 1. USE VERIFIED H3 PAIRED DATA
# ============================================================

stopifnot(
  exists("fungi_C_chao_logratio"),
  exists("model_chao_fun_C_logratio")
)

plot_dat <- fungi_C_chao_logratio %>%
  dplyr::mutate(
    Ecological_condition = factor(
      Ecological_condition,
      levels = c("Natural", "BESE", "BARE")
    )
  )


# ============================================================
# 2. MODEL-BASED ESTIMATES
# ============================================================

emm_plot <- emmeans::emmeans(
  model_chao_fun_C_logratio,
  ~ Ecological_condition
)

emm_df <- as.data.frame(
  summary(
    emm_plot,
    infer = c(TRUE, FALSE)
  )
)


# ============================================================
# 3. WITHIN-CONDITION TESTS AGAINST ZERO
# Holm correction
# ============================================================

test_df <- as.data.frame(
  summary(
    test(
      emm_plot,
      null = 0,
      adjust = "holm"
    )
  )
) %>%
  dplyr::rename(
    p_Holm = p.value
  ) %>%
  dplyr::mutate(
    sig = dplyr::case_when(
      p_Holm < 0.001 ~ "***",
      p_Holm < 0.01  ~ "**",
      p_Holm < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  )


# ============================================================
# 4. JOIN ESTIMATES + ADJUSTED TESTS
# ============================================================

emm_df <- emm_df %>%
  dplyr::left_join(
    test_df %>%
      dplyr::select(
        Ecological_condition,
        p_Holm,
        sig
      ),
    by = "Ecological_condition"
  )


# ============================================================
# 5. CHECK EXACT VALUES USED IN FIGURE
# ============================================================

cat("\n========================================\n")
cat("MODEL ESTIMATES USED IN FIGURE\n")
cat("========================================\n")

print(
  emm_df[
    ,
    c(
      "Ecological_condition",
      "emmean",
      "lower.CL",
      "upper.CL",
      "p_Holm",
      "sig"
    )
  ],
  row.names = FALSE
)

# ============================================================
# FIGURE — FUNGI CHAO1 PAIRED RESPONSE
# ============================================================

treatment_cols <- c(
  "Natural" = "#9BC1BC",
  "BESE"    = "#FF5A5F",
  "BARE"    = "#FFA630"
)

# Only annotate significant condition
sig_plot <- emm_df %>%
  dplyr::filter(
    p_Holm < 0.05
  )

top_y <- max(
  c(
    plot_dat$Log_ratio_Chao1,
    emm_df$upper.CL
  ),
  na.rm = TRUE
)

sig_plot$sig_y <- top_y + 0.25


p_fun_chao_response <- ggplot(
  plot_dat,
  aes(
    x = Ecological_condition,
    y = Log_ratio_Chao1
  )
) +
  
  # zero = no RZ-BS difference
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.65,
    colour = "grey45"
  ) +
  
  # raw paired responses
  geom_point(
    aes(
      fill = Ecological_condition
    ),
    shape = 21,
    colour = "black",
    stroke = 0.45,
    size = 2.7,
    alpha = 0.55,
    position = position_jitter(
      width = 0.09,
      height = 0,
      seed = 123
    )
  ) +
  
  # model 95% CI
  geom_errorbar(
    data = emm_df,
    aes(
      x = Ecological_condition,
      ymin = lower.CL,
      ymax = upper.CL
    ),
    inherit.aes = FALSE,
    width = 0.10,
    linewidth = 0.85,
    colour = "black"
  ) +
  
  # model mean
  geom_point(
    data = emm_df,
    aes(
      x = Ecological_condition,
      y = emmean,
      fill = Ecological_condition
    ),
    inherit.aes = FALSE,
    shape = 21,
    colour = "black",
    stroke = 1,
    size = 4.5
  ) +
  
  # significance — BESE only
  geom_text(
    data = sig_plot,
    aes(
      x = Ecological_condition,
      y = sig_y,
      label = sig
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  
  scale_fill_manual(
    values = treatment_cols,
    guide = "none"
  ) +
  
  labs(
    x = NULL,
    y = expression(
      "Rhizosphere richness response  ln(RZ / BS)"
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text.x = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 11
    ),
    axis.line = element_line(
      linewidth = 0.7
    ),
    axis.ticks = element_line(
      linewidth = 0.6
    ),
    plot.margin = margin(
      10, 15, 8, 10
    )
  )

print(p_fun_chao_response)


#BOXPLOT
# ============================================================
# H3 — FUNGI CHAO1
# CLEAN BOXplot of paired ln(RZ / BS) response
# ============================================================

p_fun_chao_box <- ggplot(
  plot_dat,
  aes(
    x = Ecological_condition,
    y = Log_ratio_Chao1,
    fill = Ecological_condition
  )
) +
  
  # zero = no rhizosphere effect
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.65,
    colour = "grey45"
  ) +
  
  # paired-response distributions
  geom_boxplot(
    width = 0.48,
    colour = "black",
    linewidth = 0.8,
    alpha = 0.85,
    outlier.shape = 21,
    outlier.size = 2.3,
    outlier.stroke = 0.5
  ) +
  
  # significance of within-condition response
  geom_text(
    data = sig_plot,
    aes(
      x = Ecological_condition,
      y = sig_y,
      label = sig
    ),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  
  scale_fill_manual(
    values = treatment_cols,
    guide = "none"
  ) +
  
  labs(
    x = NULL,
    y = expression(
      "Rhizosphere richness response  ln(RZ / BS)"
    )
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      size = 11,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 11
    ),
    axis.line = element_line(
      linewidth = 0.7
    ),
    axis.ticks = element_line(
      linewidth = 0.6
    )
  )

print(p_fun_chao_box)

# ============================================================
# H3 — PROKARYOTIC GENUS ANCOM-BC2 PLOT
#
# SAME VISUAL FORMAT AS H1 / H2
#
# Panel order:
#   BESE | BARE | Natural
#
# Negative LFC = higher estimated abundance in Bulk soil
# Positive LFC = higher estimated abundance in Rhizosphere
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# 1. PREPARE DATA
# ============================================================

plot_C <- robust_summary_C %>%
  mutate(
    Ecological_condition = factor(
      Ecological_condition,
      levels = c(
        "BESE",
        "BARE",
        "Natural"
      )
    ),
    
    # Archaeal genera
    # In our robust H3 results, Halobacterota = Archaea
    Genus_label = ifelse(
      Phylum == "Halobacterota",
      paste0(Genus, " (A)"),
      Genus
    ),
    
    # Unique taxon identifier for independent ordering
    # inside each facet
    Panel_taxon = paste(
      Ecological_condition,
      Genus_label,
      sep = "___"
    )
  )


# ============================================================
# 2. ORDER GENERA WITHIN EACH PANEL
#
# Highest positive LFC at top;
# most negative at bottom
# ============================================================

taxon_levels <- plot_C %>%
  arrange(
    Ecological_condition,
    LFC_RZ_vs_BS
  ) %>%
  pull(
    Panel_taxon
  ) %>%
  unique()

plot_C$Panel_taxon <- factor(
  plot_C$Panel_taxon,
  levels = taxon_levels
)


# ============================================================
# 3. COLOURS
# ============================================================

condition_cols <- c(
  "Natural" = "#9BC1BC",
  "BESE"    = "#FF5A5F",
  "BARE"    = "#FFA630"
)


# ============================================================
# 4. H3 ANCOM-BC2 FIGURE
# ============================================================

p_H3_ANCOM <- ggplot(
  plot_C,
  aes(
    x = LFC_RZ_vs_BS,
    y = Panel_taxon
  )
) +
  
  # ----------------------------------------------------------
# Zero reference
# ----------------------------------------------------------

geom_vline(
  xintercept = 0,
  linetype = "dashed",
  linewidth = 0.65,
  colour = "grey50"
) +
  
  # ----------------------------------------------------------
# Stem from zero to LFC
# ----------------------------------------------------------

geom_segment(
  aes(
    x = 0,
    xend = LFC_RZ_vs_BS,
    y = Panel_taxon,
    yend = Panel_taxon
  ),
  colour = "grey65",
  linewidth = 0.75
) +
  
  # ----------------------------------------------------------
# Robust ANCOM-BC2 estimate
# ----------------------------------------------------------

geom_point(
  aes(
    fill = Ecological_condition
  ),
  shape = 21,
  size = 3.5,
  stroke = 0.8,
  colour = "black"
) +
  
  # ----------------------------------------------------------
# Panels
# restoration treatments first, Natural last
# ----------------------------------------------------------

facet_wrap(
  ~ Ecological_condition,
  nrow = 1,
  scales = "free_y"
) +
  
  # ----------------------------------------------------------
# Common LFC scale
# ----------------------------------------------------------

scale_x_continuous(
  limits = c(-6, 6),
  breaks = c(
    -6,
    -4,
    -2,
    0,
    2,
    4,
    6
  )
) +
  
  # Remove condition prefix from y-axis labels
  scale_y_discrete(
    labels = function(x) {
      sub(
        ".*___",
        "",
        x
      )
    }
  ) +
  
  scale_fill_manual(
    values = condition_cols,
    guide = "none"
  ) +
  
  # ----------------------------------------------------------
# Labels
# ----------------------------------------------------------

labs(
  title =
    "Root-associated prokaryotic genera",
  
  subtitle =
    "Robust ANCOM-BC2 bulk-soil–rhizosphere contrasts after accounting for Site",
  
  x =
    "Bias-corrected log-fold change",
  
  y = NULL,
  
  caption =
    "(A) = archaeal genus; negative values indicate higher estimated abundance in bulk soil; positive values indicate higher estimated abundance in rhizosphere."
) +
  
  # ----------------------------------------------------------
# SAME STYLE AS H1 / H2
# ----------------------------------------------------------

theme_bw(
  base_size = 12
) +
  
  theme(
    
    # Title
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0
    ),
    
    plot.subtitle = element_text(
      size = 11,
      hjust = 0,
      margin = margin(
        b = 12
      )
    ),
    
    # Grey panel strips
    strip.background = element_rect(
      fill = "grey85",
      colour = "grey40",
      linewidth = 0.6
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 12,
      colour = "black"
    ),
    
    # Italic genera
    axis.text.y = element_text(
      face = "italic",
      size = 9.5,
      colour = "grey30"
    ),
    
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),
    
    # H1/H2 grid style
    panel.grid.major.y = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_line(
      colour = "grey90",
      linewidth = 0.45
    ),
    
    # Space between facets
    panel.spacing = unit(
      1.1,
      "lines"
    ),
    
    # Caption
    plot.caption = element_text(
      size = 9,
      hjust = 0,
      margin = margin(
        t = 10
      )
    )
  )


# ============================================================
# 5. DISPLAY
# ============================================================

print(
  p_H3_ANCOM
)


# ============================================================
# 6. SAVE
# ============================================================

dir.create(
  "outputs/H3_ANCOMBC2/figures",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename =
    "outputs/H3_ANCOMBC2/figures/H3_prokaryotic_genus_ANCOMBC2.png",
  plot =
    p_H3_ANCOM,
  width =
    12,
  height =
    7.2,
  dpi =
    600
)

ggsave(
  filename =
    "outputs/H3_ANCOMBC2/figures/H3_prokaryotic_genus_ANCOMBC2.pdf",
  plot =
    p_H3_ANCOM,
  width =
    12,
  height =
    7.2
)


# ============================================================
# H3 — FUNGI — GENUS-LEVEL ANCOM-BC2
#
# Paired Bulk soil vs Rhizosphere
#
# Separate models:
#   Natural
#   BESE
#   BARE
#
# Fixed:
#   Site + Soil_habitat
#
# Random:
#   (1 | PlotID)
#
# Positive LFC = higher estimated abundance in Rhizosphere
# Negative LFC = higher estimated abundance in Bulk soil
#
# Robust = diff_robust == TRUE
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
})

set.seed(123)


# ============================================================
# 1. CHECK VERIFIED H3 OBJECTS
# ============================================================

stopifnot(
  exists("ps_fun_C"),
  exists("alpha_fungi_c")
)


# ============================================================
# 2. KEEP ONLY VERIFIED COMPLETE BS-RZ PAIRS
# ============================================================

pair_meta_fun <- as.data.frame(
  alpha_fungi_c,
  stringsAsFactors = FALSE
)

pair_ids_fun <- unique(
  as.character(
    pair_meta_fun$SampleID
  )
)

missing_ids <- setdiff(
  pair_ids_fun,
  sample_names(ps_fun_C)
)

if (length(missing_ids) > 0) {
  stop(
    "Missing paired fungal SampleIDs: ",
    paste(
      missing_ids,
      collapse = ", "
    )
  )
}

ps_fun_pair <- prune_samples(
  pair_ids_fun,
  ps_fun_C
)

ps_fun_pair <- prune_taxa(
  taxa_sums(ps_fun_pair) > 0,
  ps_fun_pair
)


# ============================================================
# 3. ATTACH VERIFIED H3 METADATA
# ============================================================

meta_match <- match(
  sample_names(ps_fun_pair),
  pair_meta_fun$SampleID
)

stopifnot(
  !anyNA(meta_match)
)

sd_fun <- as.data.frame(
  sample_data(ps_fun_pair),
  stringsAsFactors = FALSE
)

sd_fun$Site <-
  pair_meta_fun$Site[
    meta_match
  ]

sd_fun$Ecological_condition <-
  pair_meta_fun$Ecological_condition[
    meta_match
  ]

sd_fun$Soil_habitat <-
  pair_meta_fun$Soil_habitat[
    meta_match
  ]

sd_fun$PlotID <-
  pair_meta_fun$PlotID[
    meta_match
  ]

sd_fun$Site <- factor(
  sd_fun$Site,
  levels = c(
    "ACP",
    "BASS",
    "TRUG"
  )
)

sd_fun$Ecological_condition <- factor(
  sd_fun$Ecological_condition,
  levels = c(
    "Natural",
    "BESE",
    "BARE"
  )
)

sd_fun$Soil_habitat <- factor(
  sd_fun$Soil_habitat,
  levels = c(
    "Bulk soil",
    "Rhizosphere"
  )
)

sd_fun$PlotID <- factor(
  sd_fun$PlotID
)

sample_data(ps_fun_pair) <-
  sample_data(sd_fun)


# ============================================================
# 4. VERIFY PAIRING
# ============================================================

cat("\n========================================\n")
cat("H3 FUNGI — PAIRING\n")
cat("========================================\n")

pair_check_fun <- unique(
  sd_fun[
    ,
    c(
      "Site",
      "Ecological_condition",
      "PlotID"
    )
  ]
)

print(
  table(
    pair_check_fun$Ecological_condition
  )
)

cat(
  "\nTotal paired plots:",
  nrow(pair_check_fun),
  "\n"
)

cat(
  "Total samples:",
  nsamples(ps_fun_pair),
  "\n"
)

stopifnot(
  nrow(pair_check_fun) == 34,
  nsamples(ps_fun_pair) == 68
)


# ============================================================
# 5. REMOVE ASVs WITHOUT GENUS ASSIGNMENT
# ============================================================

tax_fun <- as.data.frame(
  tax_table(ps_fun_pair),
  stringsAsFactors = FALSE
)

if (!"Genus" %in% colnames(tax_fun)) {
  stop(
    "No Genus column in fungal taxonomy."
  )
}

genus_fun <- trimws(
  as.character(
    tax_fun$Genus
  )
)

keep_named_fun <- !(
  is.na(genus_fun) |
    genus_fun == "" |
    genus_fun %in% c(
      "NA",
      "Unclassified",
      "unclassified"
    )
)

ps_fun_genus <- prune_taxa(
  keep_named_fun,
  ps_fun_pair
)


# ============================================================
# 6. AGGLOMERATE TO GENUS
# ============================================================

ps_fun_genus <- tax_glom(
  ps_fun_genus,
  taxrank = "Genus",
  NArm = TRUE
)

cat("\n========================================\n")
cat("H3 FUNGI — GENUS INPUT\n")
cat("========================================\n")

cat(
  "Genera before condition-specific prevalence filtering:",
  ntaxa(ps_fun_genus),
  "\n"
)


# ============================================================
# 7. FUNCTION — ONE ECOLOGICAL CONDITION
# ============================================================

run_H3_fungal_genus <- function(
    ps,
    condition_name
) {
  
  cat("\n\n")
  cat("====================================================\n")
  cat(
    "H3 FUNGI ANCOM-BC2 —",
    condition_name,
    "\n"
  )
  cat("====================================================\n")
  
  
  # ----------------------------------------------------------
  # A. Safe condition subset
  # ----------------------------------------------------------
  
  sd_all <- methods::as(
    sample_data(ps),
    "data.frame"
  )
  
  sd_all <- data.frame(
    sd_all,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  keep_ids <- rownames(sd_all)[
    as.character(
      sd_all$Ecological_condition
    ) == condition_name
  ]
  
  ps_sub <- prune_samples(
    keep_ids,
    ps
  )
  
  ps_sub <- prune_taxa(
    taxa_sums(ps_sub) > 0,
    ps_sub
  )
  
  
  # ----------------------------------------------------------
  # B. Metadata
  # ----------------------------------------------------------
  
  sd <- methods::as(
    sample_data(ps_sub),
    "data.frame"
  )
  
  sd <- data.frame(
    sd,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  sd$Site <- factor(
    sd$Site,
    levels = c(
      "ACP",
      "BASS",
      "TRUG"
    )
  )
  
  sd$Soil_habitat <- factor(
    sd$Soil_habitat,
    levels = c(
      "Bulk soil",
      "Rhizosphere"
    )
  )
  
  sd$PlotID <- factor(
    sd$PlotID
  )
  
  
  # ----------------------------------------------------------
  # C. Verify complete pairs
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
  
  stopifnot(
    all(
      pair_tab[, "Bulk soil"] == 1
    ),
    all(
      pair_tab[, "Rhizosphere"] == 1
    )
  )
  
  cat(
    "Samples:",
    nrow(sd),
    "\n"
  )
  
  cat(
    "Paired plots:",
    nrow(pair_tab),
    "\n"
  )
  
  print(
    table(
      sd$Site,
      sd$Soil_habitat
    )
  )
  
  sample_data(ps_sub) <-
    sample_data(sd)
  
  
  # ----------------------------------------------------------
  # D. Condition-specific 10% prevalence filter
  # ----------------------------------------------------------
  
  otu <- as(
    otu_table(ps_sub),
    "matrix"
  )
  
  if (!taxa_are_rows(ps_sub)) {
    otu <- t(otu)
  }
  
  prevalence <- rowSums(
    otu > 0
  )
  
  threshold <- ceiling(
    0.10 *
      nsamples(ps_sub)
  )
  
  keep_prev <-
    prevalence >= threshold
  
  cat(
    "\n10% prevalence threshold:",
    threshold,
    "samples\n"
  )
  
  cat(
    "Genera before prevalence filter:",
    ntaxa(ps_sub),
    "\n"
  )
  
  ps_sub <- prune_taxa(
    keep_prev,
    ps_sub
  )
  
  cat(
    "Genera entering ANCOM-BC2:",
    ntaxa(ps_sub),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # E. ANCOM-BC2
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
  # F. Quantitative output
  # ----------------------------------------------------------
  
  res <- as.data.frame(
    out$res,
    stringsAsFactors = FALSE
  )
  
  if (!"taxon" %in% colnames(res)) {
    res$taxon <- rownames(res)
  }
  
  lfc_col <-
    "lfc_Soil_habitatRhizosphere"
  
  se_col <-
    "se_Soil_habitatRhizosphere"
  
  p_col <-
    "p_Soil_habitatRhizosphere"
  
  q_col <-
    "q_Soil_habitatRhizosphere"
  
  robust_col <-
    "diff_robust_Soil_habitatRhizosphere"
  
  
  # ----------------------------------------------------------
  # G. Structural zeros
  # ----------------------------------------------------------
  
  zero_tab <- as.data.frame(
    out$zero_ind,
    stringsAsFactors = FALSE
  )
  
  if (!"taxon" %in% colnames(zero_tab)) {
    zero_tab$taxon <-
      rownames(zero_tab)
  }
  
  zero_cols <- setdiff(
    colnames(zero_tab),
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
  
  
  # ----------------------------------------------------------
  # H. Taxonomy
  # ----------------------------------------------------------
  
  tax_sub <- as.data.frame(
    tax_table(ps_sub),
    stringsAsFactors = FALSE
  )
  
  tax_sub$taxon <-
    rownames(tax_sub)
  
  
  # ----------------------------------------------------------
  # I. Clean result
  # ----------------------------------------------------------
  
  result_clean <- data.frame(
    taxon =
      res$taxon,
    
    LFC_RZ_vs_BS =
      res[[lfc_col]],
    
    SE =
      res[[se_col]],
    
    p =
      res[[p_col]],
    
    q_BH =
      res[[q_col]],
    
    Robust =
      as.logical(
        res[[robust_col]]
      ),
    
    stringsAsFactors = FALSE
  )
  
  
  # Taxonomy matching
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
      tax_sub[[rank]][
        tax_match
      ]
  }
  
  
  # Structural-zero matching
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
  # J. Quantitative + robust
  # ----------------------------------------------------------
  
  quantitative <- result_clean[
    !result_clean$Structural_zero_any &
      !is.na(
        result_clean$LFC_RZ_vs_BS
      ),
    ,
    drop = FALSE
  ]
  
  robust <- quantitative[
    !is.na(
      quantitative$Robust
    ) &
      quantitative$Robust,
    ,
    drop = FALSE
  ]
  
  robust <- robust[
    order(
      abs(
        robust$LFC_RZ_vs_BS
      ),
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # K. Summary
  # ----------------------------------------------------------
  
  cat("\n========================================\n")
  cat(
    condition_name,
    "— SUMMARY\n"
  )
  cat("========================================\n")
  
  cat(
    "Genera entering analysis:",
    ntaxa(ps_sub),
    "\n"
  )
  
  cat(
    "Genera with >=1 structural zero:",
    sum(
      zero_tab$Structural_zero_any,
      na.rm = TRUE
    ),
    "\n"
  )
  
  cat(
    "Genera quantitatively modelled:",
    nrow(quantitative),
    "\n"
  )
  
  cat(
    "Robust quantitative BS-RZ genera:",
    nrow(robust),
    "\n"
  )
  
  if (nrow(robust) > 0) {
    
    cat("\nROBUST GENERA\n")
    
    print(
      robust[
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
  # L. Return
  # ----------------------------------------------------------
  
  return(
    list(
      condition =
        condition_name,
      
      ps =
        ps_sub,
      
      ancombc2 =
        out,
      
      zero_table =
        zero_tab,
      
      full_results =
        result_clean,
      
      quantitative =
        quantitative,
      
      robust =
        robust,
      
      prevalence_threshold =
        threshold
    )
  )
}


# ============================================================
# 8. RUN NATURAL / BESE / BARE
# ============================================================

fun_C_Natural <- run_H3_fungal_genus(
  ps_fun_genus,
  "Natural"
)

fun_C_BESE <- run_H3_fungal_genus(
  ps_fun_genus,
  "BESE"
)

fun_C_BARE <- run_H3_fungal_genus(
  ps_fun_genus,
  "BARE"
)


# ============================================================
# 9. COMBINE ROBUST RESULTS
# ============================================================

# ============================================================
# 9. COMBINE ROBUST RESULTS
# Handles conditions with ZERO robust genera
# ============================================================

add_condition <- function(x, condition_name) {
  
  x$Ecological_condition <- rep(
    condition_name,
    nrow(x)
  )
  
  x
}

nat <- add_condition(
  fun_C_Natural$robust,
  "Natural"
)

bese <- add_condition(
  fun_C_BESE$robust,
  "BESE"
)

bare <- add_condition(
  fun_C_BARE$robust,
  "BARE"
)

robust_fun_C <- dplyr::bind_rows(
  nat,
  bese,
  bare
)

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

robust_fun_C <- robust_fun_C[
  ,
  intersect(
    keep_cols,
    colnames(robust_fun_C)
  ),
  drop = FALSE
]


# ============================================================
# 10. SUMMARY
# ============================================================

cat("\n========================================\n")
cat("H3 FUNGI — GENUS ANCOM-BC2 SUMMARY\n")
cat("========================================\n")

cat(
  "Natural robust genera:",
  nrow(fun_C_Natural$robust),
  "\n"
)

cat(
  "BESE robust genera:",
  nrow(fun_C_BESE$robust),
  "\n"
)

cat(
  "BARE robust genera:",
  nrow(fun_C_BARE$robust),
  "\n"
)

cat("\nROBUST GENERA\n")

if (nrow(robust_fun_C) > 0) {
  
  print(
    robust_fun_C,
    row.names = FALSE
  )
  
} else {
  
  cat(
    "No robust fungal genera detected in any condition.\n"
  )
}


# ============================================================
# 11. OVERLAP
# ============================================================

nat_g <- unique(
  as.character(
    fun_C_Natural$robust$Genus
  )
)

bese_g <- unique(
  as.character(
    fun_C_BESE$robust$Genus
  )
)

bare_g <- unique(
  as.character(
    fun_C_BARE$robust$Genus
  )
)

cat("\n========================================\n")
cat("ROBUST GENUS OVERLAP\n")
cat("========================================\n")

cat(
  "Shared all three:",
  length(
    Reduce(
      intersect,
      list(
        nat_g,
        bese_g,
        bare_g
      )
    )
  ),
  "\n"
)

cat(
  "Natural + BESE:",
  length(
    intersect(
      nat_g,
      bese_g
    )
  ),
  "\n"
)

cat(
  "Natural + BARE:",
  length(
    intersect(
      nat_g,
      bare_g
    )
  ),
  "\n"
)

cat(
  "BESE + BARE:",
  length(
    intersect(
      bese_g,
      bare_g
    )
  ),
  "\n"
)


# ============================================================
# 10. SUMMARY
# ============================================================

cat("\n\n========================================\n")
cat("H3 FUNGI — GENUS ANCOM-BC2 SUMMARY\n")
cat("========================================\n")

cat(
  "Natural robust genera:",
  nrow(fun_C_Natural$robust),
  "\n"
)

cat(
  "BESE robust genera:",
  nrow(fun_C_BESE$robust),
  "\n"
)

cat(
  "BARE robust genera:",
  nrow(fun_C_BARE$robust),
  "\n"
)

cat("\nROBUST GENERA\n")

print(
  robust_fun_C,
  row.names = FALSE
)


# ============================================================
# 11. OVERLAP
# ============================================================

nat_g <- unique(
  fun_C_Natural$robust$Genus
)

bese_g <- unique(
  fun_C_BESE$robust$Genus
)

bare_g <- unique(
  fun_C_BARE$robust$Genus
)

cat("\n========================================\n")
cat("ROBUST GENUS OVERLAP\n")
cat("========================================\n")

cat(
  "Shared all three:",
  length(
    Reduce(
      intersect,
      list(
        nat_g,
        bese_g,
        bare_g
      )
    )
  ),
  "\n"
)

cat(
  "Natural + BESE:",
  length(
    intersect(
      nat_g,
      bese_g
    )
  ),
  "\n"
)

cat(
  "Natural + BARE:",
  length(
    intersect(
      nat_g,
      bare_g
    )
  ),
  "\n"
)

cat(
  "BESE + BARE:",
  length(
    intersect(
      bese_g,
      bare_g
    )
  ),
  "\n"
)


# ============================================================
# 12. SAVE
# ============================================================

dir.create(
  "outputs/H3_ANCOMBC2/fungi/tables",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H3_ANCOMBC2/fungi/models",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  robust_fun_C,
  "outputs/H3_ANCOMBC2/fungi/tables/H3_fungi_genus_robust_BS_RZ.csv",
  row.names = FALSE
)

write.csv(
  fun_C_Natural$full_results,
  "outputs/H3_ANCOMBC2/fungi/tables/H3_fungi_Natural_genus_all.csv",
  row.names = FALSE
)

write.csv(
  fun_C_BESE$full_results,
  "outputs/H3_ANCOMBC2/fungi/tables/H3_fungi_BESE_genus_all.csv",
  row.names = FALSE
)

write.csv(
  fun_C_BARE$full_results,
  "outputs/H3_ANCOMBC2/fungi/tables/H3_fungi_BARE_genus_all.csv",
  row.names = FALSE
)

saveRDS(
  list(
    Natural =
      fun_C_Natural,
    
    BESE =
      fun_C_BESE,
    
    BARE =
      fun_C_BARE,
    
    robust_summary =
      robust_fun_C,
    
    ANCOMBC_version =
      as.character(
        packageVersion(
          "ANCOMBC"
        )
      ),
    
    sessionInfo =
      sessionInfo()
  ),
  "outputs/H3_ANCOMBC2/fungi/models/H3_fungi_genus_paired_ANCOMBC2.rds"
)

cat("\n========================================\n")
cat("H3 FUNGI ANCOM-BC2 COMPLETE\n")
cat("========================================\n")

