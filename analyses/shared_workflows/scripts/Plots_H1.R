# ============================================================
# H1 — BULK-SOIL ALPHA DIVERSITY FIGURE
# FINAL PLOT: CHAO1 RICHNESS
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ============================================================
# 1. LOAD VERIFIED PHYLOSEQ OBJECTS
# ============================================================

ps_bacteria <- readRDS("phyloseq_16S_corrected_metadata copy.rds")
ps_fungi    <- readRDS("phyloseq_ITS_corrected_metadata copy.rds")

stopifnot(
  nsamples(ps_bacteria) == 103,
  ntaxa(ps_bacteria) == 32530,
  nsamples(ps_fungi) == 87,
  ntaxa(ps_fungi) == 11176
)

# ============================================================
# 2. SUBSET H1 — BULK SOIL
# ============================================================

ps_bac_A <- subset_samples(
  ps_bacteria,
  SoilPortion == "BS" &
    Treatment %in% c("NATURAL", "BESE", "BARE", "DEGRADED")
)
ps_bac_A <- prune_taxa(taxa_sums(ps_bac_A) > 0, ps_bac_A)

ps_fun_A <- subset_samples(
  ps_fungi,
  SoilPortion == "BS" &
    Treatment %in% c("NATURAL", "BESE", "BARE", "DEGRADED")
)
ps_fun_A <- prune_taxa(taxa_sums(ps_fun_A) > 0, ps_fun_A)

stopifnot(
  nsamples(ps_bac_A) == 58,
  ntaxa(ps_bac_A) == 28358,
  nsamples(ps_fun_A) == 46,
  ntaxa(ps_fun_A) == 6495
)

# ============================================================
# 3. SAFE FUNCTION TO BUILD ALPHA TABLES
# ============================================================

make_alpha_table <- function(ps_object) {
  
  rich <- phyloseq::estimate_richness(
    ps_object,
    measures = c("Chao1", "Shannon")
  )
  
  rownames(rich) <- sub(
    "^X(?=\\d)",
    "",
    rownames(rich),
    perl = TRUE
  )
  
  meta <- as(sample_data(ps_object), "data.frame")
  
  stopifnot(
    setequal(rownames(rich), rownames(meta))
  )
  
  meta <- meta[rownames(rich), , drop = FALSE]
  
  stopifnot(
    identical(rownames(rich), rownames(meta))
  )
  
  alpha <- cbind(rich, meta)
  
  alpha$Treatment <- factor(
    alpha$Treatment,
    levels = c("NATURAL", "BESE", "BARE", "DEGRADED"),
    labels = c("Natural", "BESE", "BARE", "Degraded")
  )
  
  alpha$Site <- factor(
    alpha$Site,
    levels = c("ACP", "BASS", "TRUG")
  )
  
  alpha$SampleID <- rownames(alpha)
  
  alpha
}

alpha_bacteria <- make_alpha_table(ps_bac_A)
alpha_fungi    <- make_alpha_table(ps_fun_A)

# ============================================================
# 4. FINAL COLOURS
# ============================================================

treatment_cols <- c(
  Natural  = "#9BC1BC",
  BESE     = "#FF5A5F",
  BARE     = "#FFA630",
  Degraded = "#5D576B"
)

# ============================================================
# 5. HELPER TO ADD A BRACKET
# ============================================================

add_bracket <- function(p, x1, x2, y, label, text_y = NULL, line_size = 0.8) {
  
  if (is.null(text_y)) {
    text_y <- y
  }
  
  tick <- 0.03 * y
  
  p +
    annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = line_size) +
    annotate("segment", x = x1, xend = x1, y = y, yend = y - tick, linewidth = line_size) +
    annotate("segment", x = x2, xend = x2, y = y, yend = y - tick, linewidth = line_size) +
    annotate("text", x = mean(c(x1, x2)), y = text_y, label = label, size = 5.5, fontface = "bold")
}

# ============================================================
# 6. BACTERIA PANEL
# ============================================================

ymax_bac   <- max(alpha_bacteria$Chao1, na.rm = TRUE)
yrange_bac <- ymax_bac - min(alpha_bacteria$Chao1, na.rm = TRUE)

plot_bac <- ggplot(
  alpha_bacteria,
  aes(x = Treatment, y = Chao1, fill = Treatment)
) +
  geom_boxplot(
    width = 0.62,
    linewidth = 0.9,
    colour = "black",
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_jitter(
    width = 0.10,
    shape = 21,
    size = 3.0,
    stroke = 0.8,
    colour = "black",
    alpha = 0.9
  ) +
  scale_fill_manual(values = treatment_cols, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25))
  ) +
  labs(
    title = "Bacteria",
    x = NULL,
    y = "Chao1 richness"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black")
  ) +
  annotate(
    "text",
    x = 2.5,
    y = ymax_bac + 0.18 * yrange_bac,
    label = "Natural > mean(BESE, BARE, Degraded) *",
    size = 4.0,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = 2.5,
    y = ymax_bac + 0.12 * yrange_bac,
    label = "p = 0.0248",
    size = 3.6
  )

# ============================================================
# 7. FUNGI PANEL
# ============================================================

ymax_fun   <- max(alpha_fungi$Chao1, na.rm = TRUE)
ymin_fun   <- min(alpha_fungi$Chao1, na.rm = TRUE)
yrange_fun <- ymax_fun - ymin_fun

plot_fun <- ggplot(
  alpha_fungi,
  aes(x = Treatment, y = Chao1, fill = Treatment)
) +
  geom_boxplot(
    width = 0.62,
    linewidth = 0.9,
    colour = "black",
    outlier.shape = NA,
    alpha = 0.9
  ) +
  geom_jitter(
    width = 0.10,
    shape = 21,
    size = 3.0,
    stroke = 0.8,
    colour = "black",
    alpha = 0.9
  ) +
  scale_fill_manual(values = treatment_cols, drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.42))
  ) +
  labs(
    title = "Fungi",
    x = NULL,
    y = "Chao1 richness"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(colour = "black", face = "bold"),
    axis.text.y = element_text(colour = "black")
  )

# lower bracket: Natural vs BESE
plot_fun <- add_bracket(
  plot_fun,
  x1 = 1,
  x2 = 2,
  y = ymax_fun + 0.10 * yrange_fun,
  label = "**",
  text_y = ymax_fun + 0.14 * yrange_fun
)

# upper bracket: Natural vs BARE
plot_fun <- add_bracket(
  plot_fun,
  x1 = 1,
  x2 = 3,
  y = ymax_fun + 0.20 * yrange_fun,
  label = "**",
  text_y = ymax_fun + 0.24 * yrange_fun
)

# move explanatory significance text higher
plot_fun <- plot_fun +
  annotate(
    "text",
    x = 2.15,
    y = ymax_fun + 0.38 * yrange_fun,
    label = "Natural > BESE (p = 0.0056)",
    size = 3.4
  ) +
  annotate(
    "text",
    x = 2.15,
    y = ymax_fun + 0.31 * yrange_fun,
    label = "Natural > BARE (p = 0.0084)",
    size = 3.4
  )

# ============================================================
# 8. COMBINE
# ============================================================

final_plot <- (plot_bac | plot_fun) +
  plot_annotation(tag_levels = "A")

print(final_plot)

# ============================================================
# 9. SAVE FIGURE
# ============================================================

dir.create("figures/H1_alpha", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/H1_alpha/figure_data", recursive = TRUE, showWarnings = FALSE)

ggsave(
  "figures/H1_alpha/H1_bulk_Chao1_richness_final.png",
  final_plot,
  width = 9,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "figures/H1_alpha/H1_bulk_Chao1_richness_final.pdf",
  final_plot,
  width = 9,
  height = 5,
  units = "in"
)

write.csv(
  alpha_bacteria,
  "outputs/H1_alpha/figure_data/H1_bacteria_alpha_plot_data.csv",
  row.names = FALSE
)

write.csv(
  alpha_fungi,
  "outputs/H1_alpha/figure_data/H1_fungi_alpha_plot_data.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "outputs/H1_alpha/figure_data/sessionInfo_H1_alpha_figure.txt"
)

cat("\n========================================\n")
cat("H1 ALPHA FIGURE COMPLETE\n")
cat("========================================\n")



# ============================================================
# H1 — BULK-SOIL BETA DIVERSITY
# PCoA OF BRAY-CURTIS DISSIMILARITY
#
# Main figure:
#   colour = ecological condition
#
# Diagnostic figure:
#   colour = ecological condition
#   shape  = Site
#
# PCoA = descriptive visualization only.
# Statistical inference will be PERMANOVA + PERMDISP.
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(patchwork)
})


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


# Verify correct H1 datasets
stopifnot(
  nsamples(ps_bac_A) == 58,
  ntaxa(ps_bac_A) == 28358,
  nsamples(ps_fun_A) == 46,
  ntaxa(ps_fun_A) == 6495
)


cat("\n========================================\n")
cat("H1 — PCoA INPUT\n")
cat("========================================\n")

cat(
  "Bacteria:",
  nsamples(ps_bac_A),
  "samples |",
  ntaxa(ps_bac_A),
  "detected ASVs\n"
)

cat(
  "Fungi:",
  nsamples(ps_fun_A),
  "samples |",
  ntaxa(ps_fun_A),
  "detected ASVs\n"
)


# ============================================================
# 3. FINAL TREATMENT COLOURS
# ============================================================

treatment_cols <- c(
  Natural  = "#9BC1BC",
  BESE     = "#FF5A5F",
  BARE     = "#FFA630",
  Degraded = "#5D576B"
)


# Diagnostic Site shapes
# All support both fill + black outline

site_shapes <- c(
  ACP  = 21,
  BASS = 22,
  TRUG = 24
)


# ============================================================
# 4. FUNCTION — BRAY-CURTIS + PCoA
# ============================================================

make_pcoa <- function(ps_object, kingdom) {
  
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
  
  # Exact sample order
  comm <- comm[
    sample_names(ps_rel),
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(comm),
      sample_names(ps_rel)
    )
  )
  
  
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
  
  
  # ----------------------------------------------------------
  # Bray-Curtis dissimilarity
  # ----------------------------------------------------------
  
  bray <- vegan::vegdist(
    comm,
    method = "bray"
  )
  
  
  # ----------------------------------------------------------
  # PCoA
  #
  # Lingoes correction is used because Bray-Curtis can produce
  # negative eigenvalues.
  # ----------------------------------------------------------
  
  ord <- vegan::wcmdscale(
    bray,
    k = 2,
    eig = TRUE,
    add = "lingoes"
  )
  
  
  # ----------------------------------------------------------
  # Extract coordinates
  # ----------------------------------------------------------
  
  coords <- as.data.frame(
    ord$points[, 1:2, drop = FALSE]
  )
  
  names(coords) <- c(
    "PCoA1",
    "PCoA2"
  )
  
  coords$SampleID <- rownames(coords)
  
  
  # ----------------------------------------------------------
  # Variance represented by axes
  # Use positive corrected eigenvalues
  # ----------------------------------------------------------
  
  positive_eig <- ord$eig[
    ord$eig > 0
  ]
  
  axis1_percent <-
    100 * ord$eig[1] / sum(positive_eig)
  
  axis2_percent <-
    100 * ord$eig[2] / sum(positive_eig)
  
  
  # ----------------------------------------------------------
  # Attach metadata safely
  # ----------------------------------------------------------
  
  meta$SampleID <- rownames(meta)
  
  scores <- coords %>%
    left_join(
      meta,
      by = "SampleID"
    ) %>%
    mutate(
      Treatment = factor(
        Treatment,
        levels = c(
          "NATURAL",
          "BESE",
          "BARE",
          "DEGRADED"
        ),
        labels = c(
          "Natural",
          "BESE",
          "BARE",
          "Degraded"
        )
      ),
      Site = factor(
        Site,
        levels = c(
          "ACP",
          "BASS",
          "TRUG"
        )
      ),
      Kingdom = kingdom
    )
  
  
  # Critical checks
  stopifnot(
    nrow(scores) == nsamples(ps_object),
    sum(is.na(scores$Treatment)) == 0,
    sum(is.na(scores$Site)) == 0
  )
  
  
  axis_table <- data.frame(
    Kingdom = kingdom,
    Axis = c(
      "PCoA1",
      "PCoA2"
    ),
    Percent = c(
      axis1_percent,
      axis2_percent
    )
  )
  
  
  list(
    ps_relative = ps_rel,
    community_matrix = comm,
    bray = bray,
    ordination = ord,
    scores = scores,
    axis_table = axis_table,
    axis1 = axis1_percent,
    axis2 = axis2_percent
  )
}


# ============================================================
# 5. RUN PCoA
# ============================================================

pcoa_bac <- make_pcoa(
  ps_bac_A,
  "Bacteria"
)

pcoa_fun <- make_pcoa(
  ps_fun_A,
  "Fungi"
)


# ============================================================
# 6. CHECK AXIS VARIANCE
# ============================================================

cat("\n========================================\n")
cat("PCoA VARIANCE EXPLAINED\n")
cat("========================================\n")

print(
  bind_rows(
    pcoa_bac$axis_table,
    pcoa_fun$axis_table
  )
)


# ============================================================
# 7. MAIN PLOT FUNCTION
#    TREATMENT ONLY
# ============================================================

plot_pcoa_treatment <- function(result, title) {
  
  ggplot(
    result$scores,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = Treatment
    )
  ) +
    
    geom_point(
      shape = 21,
      size = 4.0,
      stroke = 0.9,
      colour = "black",
      alpha = 0.90
    ) +
    
    scale_fill_manual(
      values = treatment_cols,
      drop = FALSE
    ) +
    
    coord_equal() +
    
    labs(
      title = title,
      
      x = paste0(
        "PCoA1 (",
        round(result$axis1, 1),
        "%)"
      ),
      
      y = paste0(
        "PCoA2 (",
        round(result$axis2, 1),
        "%)"
      ),
      
      fill = "Ecological condition"
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text = element_text(
        colour = "black"
      ),
      
      legend.title = element_text(
        face = "bold"
      ),
      
      legend.position = "bottom"
    )
}


# ============================================================
# 8. MAIN PCoA PANELS
# ============================================================

plot_pcoa_bac <- plot_pcoa_treatment(
  pcoa_bac,
  "Bacteria"
)

plot_pcoa_fun <- plot_pcoa_treatment(
  pcoa_fun,
  "Fungi"
)


# ============================================================
# 9. COMBINE MAIN FIGURE
# ============================================================

plot_H1_PCoA <- (
  plot_pcoa_bac |
    plot_pcoa_fun
) +
  patchwork::plot_layout(
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

plot_H1_PCoA <- plot_H1_PCoA &
  theme(
    legend.position = "bottom"
  )

print(
  plot_H1_PCoA
)


# ============================================================
# 10. DIAGNOSTIC VERSION
#     TREATMENT COLOUR + SITE SHAPE
# ============================================================

plot_pcoa_diagnostic <- function(result, title) {
  
  ggplot(
    result$scores,
    aes(
      x = PCoA1,
      y = PCoA2,
      fill = Treatment,
      shape = Site
    )
  ) +
    
    geom_point(
      size = 4.0,
      stroke = 0.9,
      colour = "black",
      alpha = 0.90
    ) +
    
    scale_fill_manual(
      values = treatment_cols,
      drop = FALSE
    ) +
    
    scale_shape_manual(
      values = site_shapes
    ) +
    
    coord_equal() +
    
    labs(
      title = title,
      
      x = paste0(
        "PCoA1 (",
        round(result$axis1, 1),
        "%)"
      ),
      
      y = paste0(
        "PCoA2 (",
        round(result$axis2, 1),
        "%)"
      ),
      
      fill = "Ecological condition",
      shape = "Site"
    ) +
    
    theme_classic(
      base_size = 13
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold"
      ),
      
      axis.title = element_text(
        face = "bold"
      ),
      
      axis.text = element_text(
        colour = "black"
      ),
      
      legend.position = "bottom"
    )
}


diagnostic_bac <- plot_pcoa_diagnostic(
  pcoa_bac,
  "Bacteria"
)

diagnostic_fun <- plot_pcoa_diagnostic(
  pcoa_fun,
  "Fungi"
)

plot_H1_PCoA_diagnostic <- (
  diagnostic_bac |
    diagnostic_fun
) +
  patchwork::plot_annotation(
    tag_levels = "A"
  )

print(
  plot_H1_PCoA_diagnostic
)


# ============================================================
# 11. OUTPUT DIRECTORIES
# ============================================================

dir.create(
  "figures/H1_beta/PCoA",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/H1_beta/PCoA",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 12. SAVE MAIN FIGURE
# ============================================================

ggsave(
  "figures/H1_beta/PCoA/H1_bulk_PCoA_treatment.png",
  plot_H1_PCoA,
  width = 10,
  height = 5.3,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "figures/H1_beta/PCoA/H1_bulk_PCoA_treatment.pdf",
  plot_H1_PCoA,
  width = 10,
  height = 5.3,
  units = "in"
)


# ============================================================
# 13. SAVE DIAGNOSTIC SITE-SHAPE VERSION
# ============================================================

ggsave(
  "figures/H1_beta/PCoA/H1_bulk_PCoA_treatment_site_diagnostic.png",
  plot_H1_PCoA_diagnostic,
  width = 10,
  height = 5.3,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ============================================================
# 14. SAVE SCORES
# ============================================================

write.csv(
  pcoa_bac$scores,
  "outputs/H1_beta/PCoA/H1_bacteria_PCoA_scores.csv",
  row.names = FALSE
)

write.csv(
  pcoa_fun$scores,
  "outputs/H1_beta/PCoA/H1_fungi_PCoA_scores.csv",
  row.names = FALSE
)


# ============================================================
# 15. SAVE AXIS INFORMATION
# ============================================================

axis_table <- bind_rows(
  pcoa_bac$axis_table,
  pcoa_fun$axis_table
)

write.csv(
  axis_table,
  "outputs/H1_beta/PCoA/H1_PCoA_axis_variance.csv",
  row.names = FALSE
)


# ============================================================
# 16. SAVE DISTANCE MATRICES
# ============================================================

saveRDS(
  pcoa_bac$bray,
  "outputs/H1_beta/PCoA/H1_bacteria_BrayCurtis.rds"
)

saveRDS(
  pcoa_fun$bray,
  "outputs/H1_beta/PCoA/H1_fungi_BrayCurtis.rds"
)


# ============================================================
# 17. SAVE PCoA OBJECTS
# ============================================================

saveRDS(
  pcoa_bac$ordination,
  "outputs/H1_beta/PCoA/H1_bacteria_PCoA.rds"
)

saveRDS(
  pcoa_fun$ordination,
  "outputs/H1_beta/PCoA/H1_fungi_PCoA.rds"
)


# ============================================================
# 18. SAVE RELATIVE-ABUNDANCE PHYLOSEQ OBJECTS
# ============================================================

saveRDS(
  pcoa_bac$ps_relative,
  "outputs/H1_beta/PCoA/H1_bacteria_relative_abundance_phyloseq.rds"
)

saveRDS(
  pcoa_fun$ps_relative,
  "outputs/H1_beta/PCoA/H1_fungi_relative_abundance_phyloseq.rds"
)


# ============================================================
# 19. SESSION INFO
# ============================================================

capture.output(
  sessionInfo(),
  file =
    "outputs/H1_beta/PCoA/sessionInfo_H1_PCoA.txt"
)


cat("\n========================================\n")
cat("H1 PCoA COMPLETE\n")
cat("========================================\n")


# ============================================================
# H1 — BACTERIA — GENUS
# PLOTS FROM LOCKED ANCOM-BC2 RESULTS
#
# Produces:
# 1. Main figure: directional robust pairwise genera
# 2. Supplementary figure: structural-zero classifications
# 3. CSV tables used for plotting
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
})

# ------------------------------------------------------------
# 1. OUTPUT DIRECTORY
# ------------------------------------------------------------

outdir <- "outputs/H1_ANCOMBC2/bacteria_genus_FINAL"
figdir <- file.path(outdir, "figures")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. COLOURS
# ------------------------------------------------------------

treatment_cols <- c(
  NATURAL   = "#9BC1BC",
  BESE      = "#FF5A5F",
  BARE      = "#FFA630",
  DEGRADED  = "#5D576B"
)

# ------------------------------------------------------------
# 3. BUILD DIRECTIONAL COUNTS FOR MAIN FIGURE
# ------------------------------------------------------------

count_higher <- function(tbl, label) {
  if (is.null(tbl) || nrow(tbl) == 0) return(0L)
  sum(tbl$Higher == label, na.rm = TRUE)
}

directional_counts <- bind_rows(
  
  tibble(
    Contrast = "BESE vs Natural",
    Higher   = c("NATURAL", "BESE"),
    n        = c(
      count_higher(BESE_vs_NATURAL, "NATURAL"),
      count_higher(BESE_vs_NATURAL, "BESE")
    )
  ),
  
  tibble(
    Contrast = "BARE vs Natural",
    Higher   = c("NATURAL", "BARE"),
    n        = c(
      count_higher(BARE_vs_NATURAL, "NATURAL"),
      count_higher(BARE_vs_NATURAL, "BARE")
    )
  ),
  
  tibble(
    Contrast = "Degraded vs Natural",
    Higher   = c("NATURAL", "DEGRADED"),
    n        = c(
      count_higher(DEGRADED_vs_NATURAL, "NATURAL"),
      count_higher(DEGRADED_vs_NATURAL, "DEGRADED")
    )
  )
)

directional_counts$Contrast <- factor(
  directional_counts$Contrast,
  levels = c(
    "BESE vs Natural",
    "BARE vs Natural",
    "Degraded vs Natural"
  )
)

directional_counts$Higher <- factor(
  directional_counts$Higher,
  levels = c("NATURAL", "BESE", "BARE", "DEGRADED")
)

# total robust genera per contrast
directional_totals <- directional_counts %>%
  group_by(Contrast) %>%
  summarise(
    Total = sum(n),
    .groups = "drop"
  )

# positions for segment labels
directional_labels <- directional_counts %>%
  group_by(Contrast) %>%
  mutate(
    ypos = cumsum(n) - 0.5 * n
  ) %>%
  ungroup() %>%
  filter(n > 0)

# ------------------------------------------------------------
# 4. MAIN PLOT — DIRECTIONAL ROBUST GENERA
# ------------------------------------------------------------

plot_directional <- ggplot(
  directional_counts,
  aes(x = Contrast, y = n, fill = Higher)
) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(
    data = directional_labels,
    aes(y = ypos, label = n),
    size = 4
  ) +
  geom_text(
    data = directional_totals,
    aes(x = Contrast, y = Total + 1, label = paste0("Total = ", Total)),
    inherit.aes = FALSE,
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = treatment_cols,
    breaks = c("NATURAL", "BESE", "BARE", "DEGRADED"),
    labels = c("Natural", "BESE", "BARE", "Degraded")
  ) +
  labs(
    title = "H1 — Bulk-soil bacterial genera",
    subtitle = "Most robust genus-level differences were Natural-associated",
    x = NULL,
    y = "Number of robust differentially abundant genera",
    fill = "Higher in"
  ) +
  annotate(
    "text",
    x = 1.2, y = max(directional_totals$Total) + 4.5,
    label = "No robust genus-level differences among BESE, BARE, and Degraded",
    hjust = 0,
    size = 3.8
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

print(plot_directional)

ggsave(
  filename = file.path(figdir, "H1_Bacteria_Genus_directional_pairwise.pdf"),
  plot = plot_directional,
  width = 8.5,
  height = 5.5
)

ggsave(
  filename = file.path(figdir, "H1_Bacteria_Genus_directional_pairwise.png"),
  plot = plot_directional,
  width = 8.5,
  height = 5.5,
  dpi = 400
)

# ------------------------------------------------------------
# 5. STRUCTURAL-ZERO SUMMARY TABLE
# ------------------------------------------------------------

zero_ind <- as.data.frame(fit_A_bac_genus$zero_ind)

zero_cols <- grep(
  "^structural_zero",
  colnames(zero_ind),
  value = TRUE
)

extract_condition_name <- function(x) {
  sub("^structural_zero \\(Condition = ", "", sub("\\)$", "", x))
}

structural_zero_counts <- tibble(
  Condition = sapply(zero_cols, extract_condition_name),
  n = sapply(zero_cols, function(z) sum(zero_ind[[z]], na.rm = TRUE))
)

structural_zero_counts$Condition <- factor(
  structural_zero_counts$Condition,
  levels = c("NATURAL", "BESE", "BARE", "DEGRADED")
)

# ------------------------------------------------------------
# 6. SUPPLEMENTARY PLOT — STRUCTURAL ZEROS
# ------------------------------------------------------------

plot_structural_zero <- ggplot(
  structural_zero_counts,
  aes(x = Condition, y = n, fill = Condition)
) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(
    aes(label = n),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(
    values = treatment_cols,
    labels = c("Natural", "BESE", "BARE", "Degraded")
  ) +
  labs(
    title = "H1 — Bulk-soil bacterial genera",
    subtitle = "Condition-specific structural-zero classifications",
    x = NULL,
    y = "Number of structural-zero classifications"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  ) +
  expand_limits(y = max(structural_zero_counts$n) * 1.15)

print(plot_structural_zero)

ggsave(
  filename = file.path(figdir, "H1_Bacteria_Genus_structural_zero_summary.pdf"),
  plot = plot_structural_zero,
  width = 6.5,
  height = 5
)

ggsave(
  filename = file.path(figdir, "H1_Bacteria_Genus_structural_zero_summary.png"),
  plot = plot_structural_zero,
  width = 6.5,
  height = 5,
  dpi = 400
)

# ------------------------------------------------------------
# 7. SAVE TABLES USED FOR PLOTTING
# ------------------------------------------------------------

write.csv(
  directional_counts,
  file.path(outdir, "H1_Bacteria_Genus_directional_plot_data.csv"),
  row.names = FALSE
)

write.csv(
  directional_totals,
  file.path(outdir, "H1_Bacteria_Genus_directional_plot_totals.csv"),
  row.names = FALSE
)

write.csv(
  structural_zero_counts,
  file.path(outdir, "H1_Bacteria_Genus_structural_zero_plot_data.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. CONSOLE SUMMARY
# ------------------------------------------------------------

cat("\n========================================\n")
cat("PLOTS SAVED\n")
cat("========================================\n")

cat("\nMain directional plot:\n")
print(directional_counts)

cat("\nStructural-zero summary:\n")
print(structural_zero_counts)

cat("\nFigures saved in:\n")
cat(normalizePath(figdir), "\n")


##GENUS###

# ============================================================
# H1 — BACTERIA — GENUS
# ROBUST ANCOM-BC2 GENUS-LEVEL EFFECT-SIZE PLOT
#
# Uses existing:
# BESE_vs_NATURAL
# BARE_vs_NATURAL
# DEGRADED_vs_NATURAL
#
# Positive x = higher in Natural
# Negative x = higher in comparison condition
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

# ------------------------------------------------------------
# 1. OUTPUT DIRECTORY
# ------------------------------------------------------------

outdir <- "outputs/H1_ANCOMBC2/bacteria_genus_FINAL"
figdir <- file.path(outdir, "figures")

dir.create(
  figdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. COLOUR PALETTE
# ------------------------------------------------------------

treatment_cols <- c(
  NATURAL  = "#9BC1BC",
  BESE     = "#FF5A5F",
  BARE     = "#FFA630",
  DEGRADED = "#5D576B"
)

# ------------------------------------------------------------
# 3. COMBINE THE THREE NATURAL CONTRASTS
#
# Original ANCOM-BC2 coefficients are:
# BESE - Natural
# BARE - Natural
# Degraded - Natural
#
# We multiply LFC by -1 so the plotting direction becomes:
#
# positive = Natural higher
# negative = comparison condition higher
# ------------------------------------------------------------

plot_df <- bind_rows(
  
  BESE_vs_NATURAL %>%
    mutate(
      Contrast = "BESE vs Natural",
      Comparison = "BESE",
      Effect = -LFC
    ),
  
  BARE_vs_NATURAL %>%
    mutate(
      Contrast = "BARE vs Natural",
      Comparison = "BARE",
      Effect = -LFC
    ),
  
  DEGRADED_vs_NATURAL %>%
    mutate(
      Contrast = "Degraded vs Natural",
      Comparison = "DEGRADED",
      Effect = -LFC
    )
)

# ------------------------------------------------------------
# 4. STANDARDISE "HIGHER IN" LABELS
# ------------------------------------------------------------

plot_df <- plot_df %>%
  mutate(
    Higher_plot = case_when(
      Higher == "NATURAL"  ~ "NATURAL",
      Higher == "BESE"     ~ "BESE",
      Higher == "BARE"     ~ "BARE",
      Higher == "DEGRADED" ~ "DEGRADED"
    )
  )

# ------------------------------------------------------------
# 5. ORDER CONTRASTS
# ------------------------------------------------------------

plot_df$Contrast <- factor(
  plot_df$Contrast,
  levels = c(
    "BESE vs Natural",
    "BARE vs Natural",
    "Degraded vs Natural"
  )
)

# ------------------------------------------------------------
# 6. ORDER GENERA
#
# Genera occurring in more than one contrast are grouped
# toward the top of the figure.
# ------------------------------------------------------------

genus_order <- plot_df %>%
  group_by(Genus) %>%
  summarise(
    N_contrasts = n(),
    Mean_effect = mean(abs(Effect)),
    .groups = "drop"
  ) %>%
  arrange(
    N_contrasts,
    Mean_effect
  ) %>%
  pull(Genus)

plot_df$Genus <- factor(
  plot_df$Genus,
  levels = genus_order
)

# ------------------------------------------------------------
# 7. DETERMINE SYMMETRICAL X LIMITS
# ------------------------------------------------------------

xmax <- ceiling(
  max(abs(plot_df$Effect), na.rm = TRUE)
)

# ------------------------------------------------------------
# 8. MAIN GENUS-LEVEL FIGURE
# ------------------------------------------------------------

p_genus <- ggplot(
  plot_df,
  aes(
    x = Effect,
    y = Genus,
    fill = Higher_plot
  )
) +
  
  # zero = no directional difference
  geom_vline(
    xintercept = 0,
    linewidth = 0.5,
    linetype = "dashed",
    colour = "grey40"
  ) +
  
  geom_segment(
    aes(
      x = 0,
      xend = Effect,
      y = Genus,
      yend = Genus
    ),
    linewidth = 0.55,
    colour = "grey65"
  ) +
  
  geom_point(
    shape = 21,
    size = 4.2,
    stroke = 0.8,
    colour = "black"
  ) +
  
  facet_wrap(
    ~ Contrast,
    nrow = 1
  ) +
  
  scale_fill_manual(
    values = treatment_cols,
    breaks = c(
      "NATURAL",
      "BESE",
      "BARE",
      "DEGRADED"
    ),
    labels = c(
      "Natural",
      "BESE",
      "BARE",
      "Degraded"
    )
  ) +
  
  scale_x_continuous(
    limits = c(-xmax, xmax),
    breaks = seq(-xmax, xmax, by = 2)
  ) +
  
  labs(
    title = "Bulk-soil bacterial genera",
    subtitle = "Natural-associated genera dominated robust differential abundance",
    x = "Bias-corrected log-fold change\nNatural \u2212 comparison condition",
    y = NULL,
    fill = "Higher estimated\nabundance"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    
    plot.subtitle = element_text(
      size = 11
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    strip.background = element_blank(),
    
    axis.text.y = element_text(
      face = "italic",
      size = 9
    ),
    
    axis.text.x = element_text(
      size = 10
    ),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 10)
    ),
    
    legend.title = element_text(
      face = "bold"
    ),
    
    legend.position = "right",
    
    panel.spacing = unit(
      1.2,
      "lines"
    )
  )

print(p_genus)

# ------------------------------------------------------------
# 9. SAVE
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    figdir,
    "H1_Bacteria_Genus_ANCOMBC2_effect_sizes.pdf"
  ),
  plot = p_genus,
  width = 12,
  height = 8
)

ggsave(
  filename = file.path(
    figdir,
    "H1_Bacteria_Genus_ANCOMBC2_effect_sizes.png"
  ),
  plot = p_genus,
  width = 12,
  height = 8,
  dpi = 400
)

# ------------------------------------------------------------
# 10. SAVE EXACT PLOTTING DATA
# ------------------------------------------------------------

write.csv(
  plot_df,
  file.path(
    outdir,
    "H1_Bacteria_Genus_ANCOMBC2_effect_size_plot_data.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("GENUS EFFECT-SIZE FIGURE COMPLETE\n")
cat("========================================\n")

cat(
  "Unique robust genera shown:",
  dplyr::n_distinct(plot_df$Genus),
  "\n"
)

cat(
  "Robust genus-condition associations:",
  nrow(plot_df),
  "\n"
)

cat(
  "\nPositive LFC = higher in Natural\n",
  "Negative LFC = higher in BESE, BARE or Degraded\n",
  sep = ""
)

cat(
  "\nSaved in:\n",
  normalizePath(figdir),
  "\n"
)



