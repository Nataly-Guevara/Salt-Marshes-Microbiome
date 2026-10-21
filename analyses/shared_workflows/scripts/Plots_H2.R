# ============================================================
# H2 — RHIZOSPHERE REASSEMBLY
# FUNGI — SIGNIFICANT ALPHA DIVERSITY PLOTS
#
# Plots:
#   1. Chao1 richness
#   2. Shannon diversity
#
# Uses:
#   verified historical ITS phyloseq object
#
# Model used for both metrics:
#   metric ~ Site + Condition
#
# Significant pairwise comparisons shown:
#   BESE vs NATURAL
#   BESE vs BARE
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(emmeans)
})

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

# ============================================================
# 1. COLOURS
# ============================================================

treatment_cols <- c(
  NATURAL  = "#9BC1BC",
  BESE     = "#FF5A5F",
  BARE     = "#FFA630"
)

# ============================================================
# 2. LOAD VERIFIED FUNGAL OBJECT
# ============================================================

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)

# ============================================================
# 3. SUBSET H2
# Rhizosphere only: NATURAL, BESE, BARE
# ============================================================

ps_fun_B <- subset_samples(
  ps_fungi,
  SoilPortion == "RZ" &
    Treatment %in% c("NATURAL", "BESE", "BARE")
)

ps_fun_B <- prune_taxa(
  taxa_sums(ps_fun_B) > 0,
  ps_fun_B
)

cat("\n========================================\n")
cat("H2 — FUNGI\n")
cat("========================================\n")
cat("Samples:", nsamples(ps_fun_B), "\n")
cat("ASVs:", ntaxa(ps_fun_B), "\n")
cat("Reads:", sum(sample_sums(ps_fun_B)), "\n")

# ============================================================
# 4. BUILD ALPHA TABLE
# ============================================================

rich <- as.data.frame(
  estimate_richness(
    ps_fun_B,
    measures = c("Chao1", "Shannon")
  )
)

rownames(rich) <- sub(
  "^X(?=\\d)",
  "",
  rownames(rich),
  perl = TRUE
)

meta <- as(
  sample_data(ps_fun_B),
  "data.frame"
)

meta$SampleID <- rownames(meta)

alpha_fun <- data.frame(
  SampleID = rownames(rich),
  Chao1 = rich$Chao1,
  Shannon = rich$Shannon,
  stringsAsFactors = FALSE
) %>%
  left_join(
    meta %>%
      tibble::rownames_to_column("RowID"),
    by = c("SampleID" = "RowID")
  )

alpha_fun$Site <- factor(
  alpha_fun$Site,
  levels = c("ACP", "BASS", "TRUG")
)

alpha_fun$Condition <- factor(
  alpha_fun$Treatment,
  levels = c("NATURAL", "BESE", "BARE")
)

alpha_fun$Condition_label <- factor(
  alpha_fun$Condition,
  levels = c("NATURAL", "BESE", "BARE"),
  labels = c("Natural", "BESE", "BARE")
)

stopifnot(nrow(alpha_fun) == 41)

cat("\nSamples by condition:\n")
print(table(alpha_fun$Condition))

# ============================================================
# 5. FIT FINAL MODELS
# ============================================================

mod_chao <- lm(
  Chao1 ~ Site + Condition,
  data = alpha_fun
)

mod_shannon <- lm(
  Shannon ~ Site + Condition,
  data = alpha_fun
)

cat("\n========================================\n")
cat("MODEL CHECK — CHAO1\n")
cat("========================================\n")
print(anova(mod_chao))

cat("\n========================================\n")
cat("MODEL CHECK — SHANNON\n")
cat("========================================\n")
print(anova(mod_shannon))

# ============================================================
# 6. EMMEANS + PAIRWISE COMPARISONS
# ============================================================

emm_chao <- emmeans(
  mod_chao,
  ~ Condition
)

emm_shannon <- emmeans(
  mod_shannon,
  ~ Condition
)

emm_chao_df <- as.data.frame(
  summary(
    emm_chao,
    infer = TRUE
  )
)

emm_shannon_df <- as.data.frame(
  summary(
    emm_shannon,
    infer = TRUE
  )
)

emm_chao_df$Condition <- factor(
  emm_chao_df$Condition,
  levels = c("NATURAL", "BESE", "BARE")
)

emm_shannon_df$Condition <- factor(
  emm_shannon_df$Condition,
  levels = c("NATURAL", "BESE", "BARE")
)

emm_chao_df$Condition_label <- factor(
  emm_chao_df$Condition,
  levels = c("NATURAL", "BESE", "BARE"),
  labels = c("Natural", "BESE", "BARE")
)

emm_shannon_df$Condition_label <- factor(
  emm_shannon_df$Condition,
  levels = c("NATURAL", "BESE", "BARE"),
  labels = c("Natural", "BESE", "BARE")
)

pairs_chao <- as.data.frame(
  summary(
    pairs(emm_chao, adjust = "holm"),
    infer = TRUE
  )
)

pairs_shannon <- as.data.frame(
  summary(
    pairs(emm_shannon, adjust = "holm"),
    infer = TRUE
  )
)

cat("\n========================================\n")
cat("PAIRWISE — CHAO1\n")
cat("========================================\n")
print(pairs_chao)

cat("\n========================================\n")
cat("PAIRWISE — SHANNON\n")
cat("========================================\n")
print(pairs_shannon)

# ============================================================
# 7. SIGNIFICANCE LABELS
# ============================================================

p_to_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# Current significant comparisons:
# NATURAL - BESE
# BESE - BARE
#
# We standardize x positions to:
# Natural = 1
# BESE    = 2
# BARE    = 3

make_sig_df <- function(metric_values, pairwise_df) {
  
  ymax <- max(metric_values, na.rm = TRUE)
  yrng <- diff(range(metric_values, na.rm = TRUE))
  
  if (yrng == 0) {
    yrng <- ymax * 0.1
  }
  
  data.frame(
    xmin = c(1, 2),
    xmax = c(2, 3),
    y = c(
      ymax + 0.10 * yrng,
      ymax + 0.20 * yrng
    ),
    label = c(
      p_to_stars(
        pairwise_df$p.value[
          pairwise_df$contrast == "NATURAL - BESE"
        ]
      ),
      p_to_stars(
        pairwise_df$p.value[
          pairwise_df$contrast == "BESE - BARE"
        ]
      )
    ),
    p_text = c(
      paste0(
        "p = ",
        formatC(
          pairwise_df$p.value[
            pairwise_df$contrast == "NATURAL - BESE"
          ],
          format = "f",
          digits = 4
        )
      ),
      paste0(
        "p = ",
        formatC(
          pairwise_df$p.value[
            pairwise_df$contrast == "BESE - BARE"
          ],
          format = "f",
          digits = 4
        )
      )
    )
  )
}

sig_chao <- make_sig_df(
  alpha_fun$Chao1,
  pairs_chao
)

sig_shannon <- make_sig_df(
  alpha_fun$Shannon,
  pairs_shannon
)

# ============================================================
# 8. HELPER FUNCTION FOR PLOTS
# ============================================================

make_alpha_plot <- function(
    data,
    emmeans_df,
    sig_df,
    yvar,
    ylab,
    title_text
) {
  
  ymax <- max(data[[yvar]], na.rm = TRUE)
  yrng <- diff(range(data[[yvar]], na.rm = TRUE))
  
  if (yrng == 0) {
    yrng <- ymax * 0.1
  }
  
  ggplot(
    data,
    aes(
      x = Condition_label,
      y = .data[[yvar]],
      fill = Condition
    )
  ) +
    geom_boxplot(
      width = 0.65,
      alpha = 0.45,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.5
    ) +
    geom_jitter(
      width = 0.10,
      size = 2.4,
      alpha = 0.85,
      shape = 21,
      colour = "black",
      stroke = 0.4
    ) +
    geom_point(
      data = emmeans_df,
      aes(
        x = Condition_label,
        y = emmean
      ),
      inherit.aes = FALSE,
      shape = 23,
      size = 4,
      fill = "white",
      colour = "black",
      stroke = 0.9
    ) +
    geom_errorbar(
      data = emmeans_df,
      aes(
        x = Condition_label,
        y = emmean,
        ymin = lower.CL,
        ymax = upper.CL
      ),
      inherit.aes = FALSE,
      width = 0.08,
      linewidth = 0.8,
      colour = "black"
    ) +
    geom_segment(
      data = sig_df,
      aes(
        x = xmin,
        xend = xmax,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.7
    ) +
    geom_segment(
      data = sig_df,
      aes(
        x = xmin,
        xend = xmin,
        y = y - 0.025 * yrng,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.7
    ) +
    geom_segment(
      data = sig_df,
      aes(
        x = xmax,
        xend = xmax,
        y = y - 0.025 * yrng,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.7
    ) +
    geom_text(
      data = sig_df,
      aes(
        x = (xmin + xmax) / 2,
        y = y + 0.04 * yrng,
        label = label
      ),
      inherit.aes = FALSE,
      size = 5.2,
      fontface = "bold"
    ) +
    geom_text(
      data = sig_df,
      aes(
        x = (xmin + xmax) / 2,
        y = y + 0.085 * yrng,
        label = p_text
      ),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    scale_fill_manual(
      values = treatment_cols,
      labels = c(
        NATURAL = "Natural",
        BESE = "BESE",
        BARE = "BARE"
      )
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = ylab
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.x = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    expand_limits(
      y = ymax + 0.32 * yrng
    )
}

# ============================================================
# 9. MAKE PLOTS
# ============================================================

plot_fun_chao <- make_alpha_plot(
  data = alpha_fun,
  emmeans_df = emm_chao_df,
  sig_df = sig_chao,
  yvar = "Chao1",
  ylab = "Chao1 richness",
  title_text = "Fungi"
)

plot_fun_shannon <- make_alpha_plot(
  data = alpha_fun,
  emmeans_df = emm_shannon_df,
  sig_df = sig_shannon,
  yvar = "Shannon",
  ylab = "Shannon diversity",
  title_text = "Fungi"
)

print(plot_fun_chao)
print(plot_fun_shannon)

# ============================================================
# 10. SAVE OUTPUTS
# ============================================================

fig_dir <- "outputs/H2_alpha/figures"
tab_dir <- "outputs/H2_alpha/tables"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(fig_dir, "H2_Fungi_Chao1.png"),
  plot_fun_chao,
  width = 5.4,
  height = 5.3,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "H2_Fungi_Chao1.pdf"),
  plot_fun_chao,
  width = 5.4,
  height = 5.3
)

ggsave(
  file.path(fig_dir, "H2_Fungi_Shannon.png"),
  plot_fun_shannon,
  width = 5.4,
  height = 5.3,
  dpi = 400
)

ggsave(
  file.path(fig_dir, "H2_Fungi_Shannon.pdf"),
  plot_fun_shannon,
  width = 5.4,
  height = 5.3
)

write.csv(
  alpha_fun,
  file.path(tab_dir, "H2_Fungi_alpha_values_for_plot.csv"),
  row.names = FALSE
)

write.csv(
  emm_chao_df,
  file.path(tab_dir, "H2_Fungi_Chao1_emmeans.csv"),
  row.names = FALSE
)

write.csv(
  emm_shannon_df,
  file.path(tab_dir, "H2_Fungi_Shannon_emmeans.csv"),
  row.names = FALSE
)

write.csv(
  pairs_chao,
  file.path(tab_dir, "H2_Fungi_Chao1_pairwise.csv"),
  row.names = FALSE
)

write.csv(
  pairs_shannon,
  file.path(tab_dir, "H2_Fungi_Shannon_pairwise.csv"),
  row.names = FALSE
)

write.csv(
  sig_chao,
  file.path(tab_dir, "H2_Fungi_Chao1_significance_brackets.csv"),
  row.names = FALSE
)

write.csv(
  sig_shannon,
  file.path(tab_dir, "H2_Fungi_Shannon_significance_brackets.csv"),
  row.names = FALSE
)

cat("\n========================================\n")
cat("H2 FUNGI ALPHA PLOTS SAVED\n")
cat("========================================\n")
cat("Figures:", normalizePath(fig_dir), "\n")
cat("Tables:", normalizePath(tab_dir), "\n")

# ============================================================
# H2 — FUNGI ALPHA DIVERSITY
# TIGHT NATURE-STYLE FIGURE
#
# Uses existing objects:
#   alpha_fun
#   pairs_chao
#   pairs_shannon
#
# NO MODELS ARE RERUN
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ============================================================
# 1. COLOURS
# ============================================================

treatment_cols <- c(
  NATURAL = "#9BC1BC",
  BESE    = "#FF5A5F",
  BARE    = "#FFA630"
)

# ============================================================
# 2. CONDITION ORDER
# ============================================================

alpha_fun$Condition_label <- factor(
  alpha_fun$Condition,
  levels = c("NATURAL", "BESE", "BARE"),
  labels = c("Natural", "BESE", "BARE")
)

# ============================================================
# 3. FUNCTION TO CONVERT P TO STARS
# ============================================================

p_to_star <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("ns")
}

# ============================================================
# 4. EXTRACT SIGNIFICANT P VALUES
# ============================================================

p_chao_N_B <- pairs_chao$p.value[
  pairs_chao$contrast == "NATURAL - BESE"
]

p_chao_B_BA <- pairs_chao$p.value[
  pairs_chao$contrast == "BESE - BARE"
]

p_shan_N_B <- pairs_shannon$p.value[
  pairs_shannon$contrast == "NATURAL - BESE"
]

p_shan_B_BA <- pairs_shannon$p.value[
  pairs_shannon$contrast == "BESE - BARE"
]

# ============================================================
# 5. BRACKET POSITIONS — CHAO1
# ============================================================

chao_max <- max(alpha_fun$Chao1, na.rm = TRUE)
chao_range <- diff(range(alpha_fun$Chao1, na.rm = TRUE))

sig_chao <- data.frame(
  xmin = c(1, 2),
  xmax = c(2, 3),
  
  y = c(
    chao_max + 0.08 * chao_range,
    chao_max + 0.20 * chao_range
  ),
  
  label = c(
    p_to_star(p_chao_N_B),
    p_to_star(p_chao_B_BA)
  )
)

# ============================================================
# 6. BRACKET POSITIONS — SHANNON
# ============================================================

shan_max <- max(alpha_fun$Shannon, na.rm = TRUE)
shan_range <- diff(range(alpha_fun$Shannon, na.rm = TRUE))

sig_shannon <- data.frame(
  xmin = c(1, 2),
  xmax = c(2, 3),
  
  y = c(
    shan_max + 0.08 * shan_range,
    shan_max + 0.20 * shan_range
  ),
  
  label = c(
    p_to_star(p_shan_N_B),
    p_to_star(p_shan_B_BA)
  )
)

# ============================================================
# 7. SIGNIFICANCE-BRACKET FUNCTION
# ============================================================

add_brackets <- function(p, sig, yrange) {
  
  p +
    geom_segment(
      data = sig,
      aes(
        x = xmin,
        xend = xmax,
        y = y,
        yend = y
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    
    geom_segment(
      data = sig,
      aes(
        x = xmin,
        xend = xmin,
        y = y,
        yend = y - 0.025 * yrange
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    
    geom_segment(
      data = sig,
      aes(
        x = xmax,
        xend = xmax,
        y = y,
        yend = y - 0.025 * yrange
      ),
      inherit.aes = FALSE,
      linewidth = 0.55
    ) +
    
    geom_text(
      data = sig,
      aes(
        x = (xmin + xmax) / 2,
        y = y + 0.035 * yrange,
        label = label
      ),
      inherit.aes = FALSE,
      size = 4.1
    )
}

# ============================================================
# 8. COMMON THEME
# ============================================================

nature_theme <- theme_classic(base_size = 10.5) +
  
  theme(
    axis.title = element_text(
      colour = "black",
      size = 11
    ),
    
    axis.text = element_text(
      colour = "black",
      size = 10
    ),
    
    axis.text.x = element_text(
      face = "plain"
    ),
    
    axis.line = element_line(
      linewidth = 0.55,
      colour = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    
    legend.position = "none",
    
    plot.margin = margin(
      7, 7, 5, 5
    )
  )

# ============================================================
# 9. CHAO1
# ============================================================

p_chao <- ggplot(
  alpha_fun,
  aes(
    x = Condition_label,
    y = Chao1,
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
    aes(fill = Condition),
    width = 0.09,
    height = 0,
    shape = 21,
    size = 2.5,
    stroke = 0.45,
    colour = "black",
    alpha = 0.9
  ) +
  
  scale_fill_manual(
    values = treatment_cols
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.03, 0.28)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Chao1 richness"
  ) +
  
  nature_theme

p_chao <- add_brackets(
  p_chao,
  sig_chao,
  chao_range
)

# ============================================================
# 10. SHANNON
# ============================================================

p_shannon <- ggplot(
  alpha_fun,
  aes(
    x = Condition_label,
    y = Shannon,
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
    aes(fill = Condition),
    width = 0.09,
    height = 0,
    shape = 21,
    size = 2.5,
    stroke = 0.45,
    colour = "black",
    alpha = 0.9
  ) +
  
  scale_fill_manual(
    values = treatment_cols
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(0.03, 0.28)
    )
  ) +
  
  labs(
    x = NULL,
    y = "Shannon diversity"
  ) +
  
  nature_theme

p_shannon <- add_brackets(
  p_shannon,
  sig_shannon,
  shan_range
)

# ============================================================
# 11. COMBINE AS FIGURE
# ============================================================

fig_H2_fungi_alpha <- (
  p_chao +
    p_shannon
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

print(fig_H2_fungi_alpha)

# ============================================================
# 12. SAVE
# ============================================================

fig_dir <- "outputs/H2_alpha/figures"

dir.create(
  fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  file.path(
    fig_dir,
    "H2_Fungi_alpha_Nature_style.pdf"
  ),
  fig_H2_fungi_alpha,
  width = 6.6,
  height = 3.6
)

ggsave(
  file.path(
    fig_dir,
    "H2_Fungi_alpha_Nature_style.png"
  ),
  fig_H2_fungi_alpha,
  width = 6.6,
  height = 3.6,
  dpi = 600
)

cat("\n========================================\n")
cat("NATURE-STYLE FIGURE SAVED\n")
cat("========================================\n")
cat(normalizePath(fig_dir), "\n")

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

# ============================================================
# H2 — DISPLAY MISSING BETA-DIVERSITY RESULTS
# No analyses are rerun
# ============================================================

cat("\n========================================\n")
cat("PAIRWISE PERMANOVA + PERMDISP — FULL\n")
cat("========================================\n")

print(
  pairwise_beta %>%
    dplyr::select(
      Kingdom,
      Contrast,
      n,
      PERMANOVA_F,
      PERMANOVA_R2,
      PERMANOVA_p,
      PERMANOVA_p_Holm,
      PERMDISP_F,
      PERMDISP_p,
      PERMDISP_p_Holm
    ),
  n = Inf,
  width = Inf
)

cat("\n========================================\n")
cat("PCoA VARIANCE\n")
cat("========================================\n")

cat(
  "Bacteria: PCoA1 =",
  round(pcoa_bac$variance["PCoA1"], 1),
  "%; PCoA2 =",
  round(pcoa_bac$variance["PCoA2"], 1),
  "%\n"
)

cat(
  "Fungi: PCoA1 =",
  round(pcoa_fun$variance["PCoA1"], 1),
  "%; PCoA2 =",
  round(pcoa_fun$variance["PCoA2"], 1),
  "%\n"
)



# ============================================================
# H2 — PROKARYOTIC GENUS ANCOM-BC2 PLOT
# NATURAL TO THE RIGHT
#
# Assumes this object already exists in memory:
#   H2_Bacteria_Genus
#
# Uses:
#   H2_Bacteria_Genus$BESE_vs_NATURAL
#   H2_Bacteria_Genus$BARE_vs_NATURAL
#   H2_Bacteria_Genus$BARE_vs_BESE
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. OUTPUT DIRECTORY
# ------------------------------------------------------------

outdir <- "outputs/H2_ANCOMBC2/figures"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. HELPER TO SAFELY EXTRACT TABLES
# ------------------------------------------------------------

safe_extract <- function(x, contrast_label) {
  
  if (is.null(x) || nrow(x) == 0) {
    
    return(
      data.frame(
        taxon = character(0),
        Contrast = character(0),
        LFC = numeric(0),
        q = numeric(0),
        Higher = character(0),
        Genus = character(0),
        Prevalence = numeric(0),
        Total_reads = numeric(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  out <- as.data.frame(x, stringsAsFactors = FALSE)
  
  if (!"Contrast" %in% colnames(out)) {
    out$Contrast <- contrast_label
  }
  
  out
}

# ------------------------------------------------------------
# 3. EXTRACT ROBUST CONTRAST TABLES
# ------------------------------------------------------------

bese_nat <- safe_extract(
  H2_Bacteria_Genus$BESE_vs_NATURAL,
  "BESE - NATURAL"
)

bare_nat <- safe_extract(
  H2_Bacteria_Genus$BARE_vs_NATURAL,
  "BARE - NATURAL"
)

bare_bese <- safe_extract(
  H2_Bacteria_Genus$BARE_vs_BESE,
  "BARE - BESE"
)

plot_dat <- bind_rows(
  bese_nat,
  bare_nat,
  bare_bese
)

# ------------------------------------------------------------
# 4. IF ALL ARE EMPTY, STOP
# ------------------------------------------------------------

if (nrow(plot_dat) == 0) {
  stop("No robust H2 bacterial genus contrasts were found to plot.")
}

# ------------------------------------------------------------
# 5. CLEAN LABELS
# ------------------------------------------------------------

plot_dat$Higher <- toupper(plot_dat$Higher)
plot_dat$Contrast <- as.character(plot_dat$Contrast)
plot_dat$Genus <- as.character(plot_dat$Genus)

# archaeal genera detected in H2
archaea_genera <- c(
  "Halobellus",
  "Halobaculum",
  "Halapricum",
  "Salinigranum",
  "Haloarchaeobius",
  "Halorussus"
)

plot_dat <- plot_dat %>%
  mutate(
    Genus_label = ifelse(
      Genus %in% archaea_genera,
      paste0(Genus, " (A)"),
      Genus
    )
  )

# ------------------------------------------------------------
# 6. PANEL LABELS
# ------------------------------------------------------------

plot_dat <- plot_dat %>%
  mutate(
    Panel = case_when(
      Contrast == "BESE - NATURAL" ~ "BESE vs Natural",
      Contrast == "BARE - NATURAL" ~ "BARE vs Natural",
      Contrast == "BARE - BESE"    ~ "BARE vs BESE",
      TRUE ~ Contrast
    )
  )

# ------------------------------------------------------------
# 7. FLIP EFFECTS SO NATURAL IS TO THE RIGHT
#
# For BESE - NATURAL and BARE - NATURAL:
#   positive = higher in Natural
#   negative = higher in BESE/BARE
#
# For BARE - BESE:
#   positive = higher in BESE
#   negative = higher in BARE
# ------------------------------------------------------------

plot_dat <- plot_dat %>%
  mutate(
    Plot_LFC = case_when(
      Contrast %in% c("BESE - NATURAL", "BARE - NATURAL") ~ -LFC,
      Contrast == "BARE - BESE" & Higher == "BESE" ~ abs(LFC),
      Contrast == "BARE - BESE" & Higher == "BARE" ~ -abs(LFC),
      TRUE ~ LFC
    ),
    
    Higher_plot = case_when(
      Contrast %in% c("BESE - NATURAL", "BARE - NATURAL") & Higher == "NATURAL" ~ "Natural",
      Contrast %in% c("BESE - NATURAL", "BARE - NATURAL") & Higher == "BESE"    ~ "BESE",
      Contrast %in% c("BESE - NATURAL", "BARE - NATURAL") & Higher == "BARE"    ~ "BARE",
      Contrast == "BARE - BESE" & Higher == "BESE" ~ "BESE",
      Contrast == "BARE - BESE" & Higher == "BARE" ~ "BARE",
      TRUE ~ Higher
    )
  )

# ------------------------------------------------------------
# 8. ADD EMPTY PANEL MESSAGE IF BARE vs BESE HAS NO ROBUST GENERA
# ------------------------------------------------------------

if (!any(plot_dat$Panel == "BARE vs BESE")) {
  
  empty_row <- data.frame(
    taxon = NA_character_,
    Contrast = "BARE - BESE",
    LFC = NA_real_,
    q = NA_real_,
    Higher = NA_character_,
    Genus = NA_character_,
    Prevalence = NA_real_,
    Total_reads = NA_real_,
    Genus_label = "No robust genera",
    Panel = "BARE vs BESE",
    Plot_LFC = 0,
    Higher_plot = NA_character_,
    stringsAsFactors = FALSE
  )
  
  plot_dat <- bind_rows(plot_dat, empty_row)
}

# ------------------------------------------------------------
# 9. ORDER WITHIN PANEL
# ------------------------------------------------------------

plot_dat <- plot_dat %>%
  mutate(
    Placeholder = Genus_label == "No robust genera"
  )

plot_order <- plot_dat %>%
  filter(!Placeholder) %>%
  arrange(Panel, Plot_LFC) %>%
  mutate(
    Genus_panel = paste0(Genus_label, "__", Panel)
  )

if (nrow(plot_order) > 0) {
  ordered_levels <- unique(plot_order$Genus_panel)
} else {
  ordered_levels <- character(0)
}

plot_dat <- plot_dat %>%
  mutate(
    Genus_panel = ifelse(
      Placeholder,
      paste0("No robust genera__", Panel),
      paste0(Genus_label, "__", Panel)
    ),
    Genus_panel = factor(
      Genus_panel,
      levels = c(
        ordered_levels,
        "No robust genera__BARE vs BESE"
      )
    ),
    Panel = factor(
      Panel,
      levels = c(
        "BESE vs Natural",
        "BARE vs Natural",
        "BARE vs BESE"
      )
    )
  )

# ------------------------------------------------------------
# 10. COLOURS
# ------------------------------------------------------------

cols_higher <- c(
  "Natural" = "#9BC1BC",
  "BESE"    = "#FF5A5F",
  "BARE"    = "#FFA630"
)

# ------------------------------------------------------------
# 11. X LIMITS
# ------------------------------------------------------------

xmax <- max(abs(plot_dat$Plot_LFC), na.rm = TRUE)
xlim_use <- max(3, ceiling(xmax * 1.10))

# ------------------------------------------------------------
# 12. PLOT
# ------------------------------------------------------------

p_H2_ancom <- ggplot() +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6,
    colour = "grey45"
  ) +
  
  geom_segment(
    data = plot_dat %>% filter(!Placeholder),
    aes(
      x = 0,
      xend = Plot_LFC,
      y = Genus_panel,
      yend = Genus_panel
    ),
    linewidth = 1,
    colour = "grey65"
  ) +
  
  geom_point(
    data = plot_dat %>% filter(!Placeholder),
    aes(
      x = Plot_LFC,
      y = Genus_panel,
      fill = Higher_plot
    ),
    shape = 21,
    size = 5.5,
    stroke = 1,
    colour = "black"
  ) +
  
  geom_text(
    data = plot_dat %>% filter(Placeholder),
    aes(
      x = 0,
      y = Genus_panel,
      label = "No robust genera"
    ),
    colour = "grey40",
    size = 7
  ) +
  
  facet_wrap(
    ~ Panel,
    nrow = 1,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = cols_higher,
    name = "Higher estimated abundance"
  ) +
  
  scale_y_discrete(
    labels = function(x) sub("__.*$", "", x)
  ) +
  
  coord_cartesian(
    xlim = c(-xlim_use, xlim_use)
  ) +
  
  labs(
    title = "Rhizosphere prokaryotic genera",
    subtitle = "Robust ANCOM-BC2 contrasts after accounting for Site",
    x = "Bias-corrected log-fold change",
    y = NULL,
    caption = "(A) = archaeal genus; for contrasts against Natural, positive values indicate higher estimated abundance in Natural"
  ) +
  
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 24),
    plot.subtitle = element_text(size = 17),
    axis.title.x = element_text(face = "bold", size = 17),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14, face = "italic"),
    strip.text = element_text(face = "bold", size = 18),
    legend.title = element_text(face = "bold", size = 17),
    legend.text = element_text(size = 14),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

print(p_H2_ancom)

# ------------------------------------------------------------
# 13. SAVE FIGURE
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    outdir,
    "H2_rhizosphere_prokaryotic_genera_ANCOMBC2_Natural_right.png"
  ),
  plot = p_H2_ancom,
  width = 16,
  height = 8,
  dpi = 400
)

ggsave(
  filename = file.path(
    outdir,
    "H2_rhizosphere_prokaryotic_genera_ANCOMBC2_Natural_right.pdf"
  ),
  plot = p_H2_ancom,
  width = 16,
  height = 8
)

cat("\n========================================\n")
cat("H2 ANCOM PLOT SAVED\n")
cat("========================================\n")
cat(
  "PNG:",
  file.path(outdir, "H2_rhizosphere_prokaryotic_genera_ANCOMBC2_Natural_right.png"),
  "\n"
)
cat(
  "PDF:",
  file.path(outdir, "H2_rhizosphere_prokaryotic_genera_ANCOMBC2_Natural_right.pdf"),
  "\n"
)


# ============================================================
# H2 FUNGI
# DO LOW-PREVALENCE ASVs EXPLAIN HIGHER BESE RICHNESS?
#
# PRIMARY QUESTION:
# Does BESE contain more low-prevalence fungal ASVs than
# Natural or BARE after accounting for Site and library depth?
#
# Definitions:
#
# Primary low-prevalence:
#   ASVs detected in <10% of H2 samples
#   = prevalence 1-4 of 41 samples
#
# Sensitivity low-prevalence:
#   ASVs detected in 2-4 samples
#   = excludes H2 singletons
#
# Higher-prevalence:
#   ASVs detected in >=5 samples
#
# Model:
#   Negative-binomial GLM
#   richness ~ Site + Condition + log10(LibrarySize)
#
# NO RAREFACTION
# ============================================================


suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(MASS)
  library(car)
  library(emmeans)
})


# ============================================================
# 1. WORKING DIRECTORY
# ============================================================

# Original personal working directory removed. Run from the Chapter 1 analysis workspace root.
# Optionally set CHAPTER1_WORKDIR to that workspace path before running.
if (nzchar(Sys.getenv("CHAPTER1_WORKDIR"))) setwd(Sys.getenv("CHAPTER1_WORKDIR"))

outdir <- "outputs/H2_fungi_low_prevalence"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. LOAD HISTORICAL FUNGI OBJECT
# ============================================================

ps_fungi <- readRDS(
  "phyloseq_ITS_corrected_metadata copy.rds"
)


# ============================================================
# 3. H2 SUBSET
# Rhizosphere: Natural, BESE, BARE
# ============================================================

ps_B_fun <- subset_samples(
  ps_fungi,
  SoilPortion == "RZ" &
    Treatment %in% c(
      "NATURAL",
      "BESE",
      "BARE"
    )
)

ps_B_fun <- prune_taxa(
  taxa_sums(ps_B_fun) > 0,
  ps_B_fun
)


cat("\n========================================\n")
cat("H2 FUNGI INPUT\n")
cat("========================================\n")

cat(
  "Samples:",
  nsamples(ps_B_fun),
  "\n"
)

cat(
  "Detected ASVs:",
  ntaxa(ps_B_fun),
  "\n"
)

cat(
  "Reads:",
  sum(sample_sums(ps_B_fun)),
  "\n"
)


# ============================================================
# 4. RAW COUNT MATRIX
# samples x ASVs
# ============================================================

counts <- as(
  otu_table(ps_B_fun),
  "matrix"
)

if (taxa_are_rows(ps_B_fun)) {
  counts <- t(counts)
}

counts <- counts[
  sample_names(ps_B_fun),
  taxa_names(ps_B_fun),
  drop = FALSE
]


# ============================================================
# 5. METADATA
# ============================================================

meta <- as.data.frame(
  sample_data(ps_B_fun),
  stringsAsFactors = FALSE
)

meta <- meta[
  rownames(counts),
  ,
  drop = FALSE
]


meta$SampleID <- rownames(meta)

meta$Site <- factor(
  meta$Site,
  levels = c(
    "ACP",
    "BASS",
    "TRUG"
  )
)

meta$Condition <- factor(
  meta$Treatment,
  levels = c(
    "NATURAL",
    "BESE",
    "BARE"
  )
)


cat("\n========================================\n")
cat("SAMPLE DESIGN\n")
cat("========================================\n")

print(
  table(
    meta$Site,
    meta$Condition
  )
)


# ============================================================
# 6. ASV PREVALENCE ACROSS H2
# ============================================================

prevalence <- colSums(
  counts > 0
)

n_samples <- nrow(counts)

threshold_10pct <- ceiling(
  0.10 * n_samples
)


cat("\n========================================\n")
cat("PREVALENCE DEFINITION\n")
cat("========================================\n")

cat(
  "H2 samples:",
  n_samples,
  "\n"
)

cat(
  "10% prevalence threshold:",
  threshold_10pct,
  "samples\n"
)


# ============================================================
# 7. DEFINE PREVALENCE CLASSES
# ============================================================

lowprev_1to4 <- names(prevalence)[
  prevalence < threshold_10pct
]

lowprev_2to4 <- names(prevalence)[
  prevalence >= 2 &
    prevalence < threshold_10pct
]

higherprev <- names(prevalence)[
  prevalence >= threshold_10pct
]


cat("\n========================================\n")
cat("ASV PREVALENCE CLASSES\n")
cat("========================================\n")

cat(
  "Low prevalence (1-4 samples):",
  length(lowprev_1to4),
  "\n"
)

cat(
  "Low prevalence sensitivity (2-4 samples):",
  length(lowprev_2to4),
  "\n"
)

cat(
  "Higher prevalence (>=5 samples):",
  length(higherprev),
  "\n"
)

cat(
  "Check total:",
  length(lowprev_1to4) +
    length(higherprev),
  "\n"
)


# ============================================================
# 8. SAMPLE-LEVEL RICHNESS WITHIN EACH PREVALENCE CLASS
# ============================================================

count_detected <- function(taxa_vector) {
  
  if (length(taxa_vector) == 0) {
    return(
      rep(
        0,
        nrow(counts)
      )
    )
  }
  
  rowSums(
    counts[
      ,
      taxa_vector,
      drop = FALSE
    ] > 0
  )
}


sample_dat <- tibble(
  SampleID = rownames(counts),
  
  Site =
    meta$Site,
  
  Condition =
    meta$Condition,
  
  LibrarySize =
    rowSums(counts),
  
  Total_ASV_richness =
    rowSums(counts > 0),
  
  LowPrev_1to4 =
    count_detected(
      lowprev_1to4
    ),
  
  LowPrev_2to4 =
    count_detected(
      lowprev_2to4
    ),
  
  HigherPrev_5plus =
    count_detected(
      higherprev
    )
)


# ============================================================
# 9. DESCRIPTIVE SUMMARY
# ============================================================

descriptive <- sample_dat %>%
  
  group_by(
    Condition
  ) %>%
  
  summarise(
    n = n(),
    
    Mean_library =
      mean(LibrarySize),
    
    Median_library =
      median(LibrarySize),
    
    Mean_total_richness =
      mean(Total_ASV_richness),
    
    Mean_lowprev_1to4 =
      mean(LowPrev_1to4),
    
    SD_lowprev_1to4 =
      sd(LowPrev_1to4),
    
    Median_lowprev_1to4 =
      median(LowPrev_1to4),
    
    Mean_lowprev_2to4 =
      mean(LowPrev_2to4),
    
    Mean_higherprev_5plus =
      mean(HigherPrev_5plus),
    
    .groups = "drop"
  )


cat("\n========================================\n")
cat("DESCRIPTIVE SUMMARY\n")
cat("========================================\n")

print(
  descriptive,
  n = Inf,
  width = Inf
)


# ============================================================
# 10. MODEL FUNCTION
#
# Negative-binomial count model:
#
#   richness ~ Site + Condition + log10(LibrarySize)
#
# Library size is explicitly adjusted because rare-taxon
# detection is sequencing-depth sensitive.
# ============================================================

run_nb_model <- function(
    response,
    label
) {
  
  formula_use <- as.formula(
    paste0(
      response,
      " ~ Site + Condition + log10(LibrarySize)"
    )
  )
  
  
  fit <- MASS::glm.nb(
    formula_use,
    data = sample_dat
  )
  
  
  # ----------------------------------------------------------
  # Type-II likelihood-ratio tests
  # ----------------------------------------------------------
  
  anova_tab <- as.data.frame(
    car::Anova(
      fit,
      type = 2,
      test.statistic = "LR"
    )
  ) %>%
    rownames_to_column(
      "Effect"
    ) %>%
    as_tibble()
  
  
  # ----------------------------------------------------------
  # Estimated marginal means
  # Response scale
  # ----------------------------------------------------------
  
  emm <- emmeans::emmeans(
    fit,
    ~ Condition,
    type = "response"
  )
  
  
  emm_tab <- summary(
    emm,
    infer = c(
      TRUE,
      TRUE
    ),
    type = "response"
  ) %>%
    as.data.frame() %>%
    as_tibble()
  
  
  # ----------------------------------------------------------
  # Pairwise comparisons
  # Holm corrected
  # ----------------------------------------------------------
  
  pairs_tab <- summary(
    pairs(
      emm,
      adjust = "holm"
    ),
    infer = c(
      TRUE,
      TRUE
    ),
    type = "response"
  ) %>%
    as.data.frame() %>%
    as_tibble()
  
  
  # ----------------------------------------------------------
  # Pearson dispersion diagnostic
  # ----------------------------------------------------------
  
  pearson_dispersion <-
    sum(
      residuals(
        fit,
        type = "pearson"
      )^2
    ) /
    df.residual(fit)
  
  
  diagnostics <- tibble(
    Analysis = label,
    
    Theta =
      fit$theta,
    
    AIC =
      AIC(fit),
    
    Pearson_dispersion =
      pearson_dispersion,
    
    Converged =
      fit$converged
  )
  
  
  cat("\n\n========================================\n")
  cat(label, "\n")
  cat("========================================\n")
  
  
  cat("\nTYPE-II TESTS\n")
  
  print(
    anova_tab,
    n = Inf,
    width = Inf
  )
  
  
  cat("\nESTIMATED MARGINAL MEANS\n")
  
  print(
    emm_tab,
    n = Inf,
    width = Inf
  )
  
  
  cat("\nHOLM-ADJUSTED PAIRWISE COMPARISONS\n")
  
  print(
    pairs_tab,
    n = Inf,
    width = Inf
  )
  
  
  cat("\nMODEL DIAGNOSTICS\n")
  
  print(
    diagnostics,
    n = Inf,
    width = Inf
  )
  
  
  list(
    fit = fit,
    anova = anova_tab,
    emmeans = emm_tab,
    pairs = pairs_tab,
    diagnostics = diagnostics
  )
}


# ============================================================
# 11. PRIMARY TEST
# LOW-PREVALENCE ASVs = 1-4 SAMPLES
# ============================================================

rare_primary <- run_nb_model(
  response = "LowPrev_1to4",
  label = "LOW-PREVALENCE ASVs — 1 to 4 samples"
)


# ============================================================
# 12. SENSITIVITY
# EXCLUDE H2 SINGLETONS
# ============================================================

rare_sensitivity <- run_nb_model(
  response = "LowPrev_2to4",
  label = "LOW-PREVALENCE ASVs — 2 to 4 samples"
)


# ============================================================
# 13. HIGHER-PREVALENCE ASVs
# ============================================================

higher_prevalence <- run_nb_model(
  response = "HigherPrev_5plus",
  label = "HIGHER-PREVALENCE ASVs — >=5 samples"
)


# ============================================================
# 14. TAXONOMY + PREVALENCE TABLE
# ============================================================

tax <- as.data.frame(
  tax_table(ps_B_fun),
  stringsAsFactors = FALSE
)

tax <- tax[
  colnames(counts),
  ,
  drop = FALSE
]


taxa_prevalence <- tibble(
  ASV =
    colnames(counts),
  
  H2_prevalence =
    as.integer(
      prevalence[colnames(counts)]
    ),
  
  Prevalence_percent =
    100 *
    prevalence[colnames(counts)] /
    n_samples,
  
  Prevalence_class =
    case_when(
      prevalence[colnames(counts)] <
        threshold_10pct ~
        "Low prevalence (<10%)",
      
      TRUE ~
        "Higher prevalence (>=10%)"
    ),
  
  Total_reads_H2 =
    as.numeric(
      colSums(counts)
    )
) %>%
  
  bind_cols(
    as_tibble(tax)
  )


# ============================================================
# 15. CONDITION-SPECIFIC OCCUPANCY
# ============================================================

condition_prevalence <- lapply(
  levels(meta$Condition),
  
  function(cond) {
    
    ids <- rownames(meta)[
      meta$Condition == cond
    ]
    
    tibble(
      ASV =
        colnames(counts),
      
      Condition =
        cond,
      
      Samples_in_condition =
        length(ids),
      
      Samples_detected =
        as.integer(
          colSums(
            counts[
              ids,
              ,
              drop = FALSE
            ] > 0
          )
        ),
      
      Condition_prevalence_percent =
        100 *
        Samples_detected /
        Samples_in_condition
    )
  }
) %>%
  bind_rows()


condition_prevalence_wide <-
  condition_prevalence %>%
  
  dplyr::select(
    ASV,
    Condition,
    Samples_detected
  ) %>%
  
  tidyr::pivot_wider(
    names_from = Condition,
    values_from = Samples_detected,
    names_prefix = "Detected_"
  )


taxa_prevalence <- taxa_prevalence %>%
  dplyr::left_join(
    condition_prevalence_wide,
    by = "ASV"
  )

# ============================================================
# 16. LOW-PREVALENCE TAXA — CONDITION ASSOCIATION DESCRIPTIVE
#
# Counts taxa detected exclusively within one ecological
# condition. DESCRIPTIVE ONLY.
# ============================================================

lowprev_condition_summary <- taxa_prevalence %>%
  
  filter(
    H2_prevalence <
      threshold_10pct
  ) %>%
  
  mutate(
    Conditions_present =
      (Detected_NATURAL > 0) +
      (Detected_BESE > 0) +
      (Detected_BARE > 0),
    
    Exclusive_condition =
      case_when(
        
        Conditions_present == 1 &
          Detected_NATURAL > 0 ~
          "Natural only",
        
        Conditions_present == 1 &
          Detected_BESE > 0 ~
          "BESE only",
        
        Conditions_present == 1 &
          Detected_BARE > 0 ~
          "BARE only",
        
        TRUE ~
          "Shared"
      )
  ) %>%
  
  count(
    Exclusive_condition,
    name = "Low_prevalence_ASVs"
  )


cat("\n========================================\n")
cat("LOW-PREVALENCE TAXA — DESCRIPTIVE OCCUPANCY\n")
cat("========================================\n")

print(
  lowprev_condition_summary,
  n = Inf,
  width = Inf
)


# ============================================================
# 17. SAVE EVERYTHING
# ============================================================

write.csv(
  sample_dat,
  file.path(
    outdir,
    "H2_fungi_sample_low_prevalence_richness.csv"
  ),
  row.names = FALSE
)


write.csv(
  descriptive,
  file.path(
    outdir,
    "H2_fungi_low_prevalence_descriptive_summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  taxa_prevalence,
  file.path(
    outdir,
    "H2_fungi_ASV_prevalence_and_taxonomy.csv"
  ),
  row.names = FALSE
)


write.csv(
  lowprev_condition_summary,
  file.path(
    outdir,
    "H2_fungi_low_prevalence_condition_occupancy.csv"
  ),
  row.names = FALSE
)


# Primary model
write.csv(
  rare_primary$anova,
  file.path(
    outdir,
    "H2_fungi_lowprev_1to4_ANOVA.csv"
  ),
  row.names = FALSE
)

write.csv(
  rare_primary$emmeans,
  file.path(
    outdir,
    "H2_fungi_lowprev_1to4_EMMeans.csv"
  ),
  row.names = FALSE
)

write.csv(
  rare_primary$pairs,
  file.path(
    outdir,
    "H2_fungi_lowprev_1to4_pairwise_Holm.csv"
  ),
  row.names = FALSE
)


# Sensitivity
write.csv(
  rare_sensitivity$anova,
  file.path(
    outdir,
    "H2_fungi_lowprev_2to4_ANOVA.csv"
  ),
  row.names = FALSE
)

write.csv(
  rare_sensitivity$pairs,
  file.path(
    outdir,
    "H2_fungi_lowprev_2to4_pairwise_Holm.csv"
  ),
  row.names = FALSE
)


# Higher prevalence
write.csv(
  higher_prevalence$anova,
  file.path(
    outdir,
    "H2_fungi_higherprev_5plus_ANOVA.csv"
  ),
  row.names = FALSE
)

write.csv(
  higher_prevalence$pairs,
  file.path(
    outdir,
    "H2_fungi_higherprev_5plus_pairwise_Holm.csv"
  ),
  row.names = FALSE
)


saveRDS(
  list(
    primary_low_prevalence =
      rare_primary,
    
    sensitivity_low_prevalence =
      rare_sensitivity,
    
    higher_prevalence =
      higher_prevalence
  ),
  
  file.path(
    outdir,
    "H2_fungi_low_prevalence_models.rds"
  )
)


capture.output(
  sessionInfo(),
  
  file = file.path(
    outdir,
    "sessionInfo_H2_fungi_low_prevalence.txt"
  )
)


cat("\n========================================\n")
cat("H2 LOW-PREVALENCE ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "Results saved in:\n",
  normalizePath(outdir),
  "\n"
)

p_paired_prok_condition <- ggplot(
  pcoa_prok_paired,
  aes(x = PCoA1, y = PCoA2)
) +
  geom_line(
    aes(group = PlotID),
    colour = "grey75",
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      fill = Ecological_condition,
      shape = Soil_habitat
    ),
    colour = "black",
    size = 3.2,
    stroke = 0.7
  ) +
  facet_wrap(
    ~Ecological_condition,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Natural" = "#A9C9C5",
      "BESE"    = "#FF6B6B",
      "BARE"    = "#F4A340"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "Bulk soil" = 21,
      "Rhizosphere" = 24
    )
  ) +
  labs(
    x = p_pcoa_C_prok$labels$x,
    y = p_pcoa_C_prok$labels$y,
    shape = "Soil habitat"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.position = "bottom",
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )

p_paired_prok_condition



##FUNGAL B DIVERSITY
# Bray–Curtis distances for the paired fungal dataset
bray_fun_pair <- phyloseq::distance(
  ps_fun_pair,
  method = "bray"
)

# PCoA using the same approach as the prokaryotic ordination
pcoa_C_fun_pair <- vegan::wcmdscale(
  bray_fun_pair,
  eig = TRUE,
  add = TRUE
)

# Extract coordinates
pcoa_fun_pair_dat <- data.frame(
  PCoA1 = pcoa_C_fun_pair$points[, 1],
  PCoA2 = pcoa_C_fun_pair$points[, 2],
  SampleID = rownames(pcoa_C_fun_pair$points),
  row.names = NULL
)

# Extract metadata
fungal_pair_metadata <- data.frame(
  phyloseq::sample_data(ps_fun_pair)
)

fungal_pair_metadata$SampleID <- rownames(
  fungal_pair_metadata
)

# Join coordinates and metadata
pcoa_fun_pair_dat <- merge(
  pcoa_fun_pair_dat,
  fungal_pair_metadata[
    ,
    c(
      "SampleID",
      "Site",
      "Ecological_condition",
      "Soil_habitat",
      "PlotID"
    )
  ],
  by = "SampleID"
)

# Set plotting order
pcoa_fun_pair_dat$Ecological_condition <- factor(
  pcoa_fun_pair_dat$Ecological_condition,
  levels = c("Natural", "BESE", "BARE")
)

pcoa_fun_pair_dat$Soil_habitat <- factor(
  pcoa_fun_pair_dat$Soil_habitat,
  levels = c("Bulk soil", "Rhizosphere")
)

# Percentage of variation represented by each axis
fungal_positive_eigenvalues <- pcoa_C_fun_pair$eig[
  pcoa_C_fun_pair$eig > 0
]

fungal_axis_1 <- round(
  100 * pcoa_C_fun_pair$eig[1] /
    sum(fungal_positive_eigenvalues),
  1
)

fungal_axis_2 <- round(
  100 * pcoa_C_fun_pair$eig[2] /
    sum(fungal_positive_eigenvalues),
  1
)

# Paired fungal PCoA
p_paired_fun_condition <- ggplot(
  pcoa_fun_pair_dat,
  aes(x = PCoA1, y = PCoA2)
) +
  geom_line(
    aes(group = PlotID),
    colour = "grey75",
    linewidth = 0.6
  ) +
  geom_point(
    aes(
      fill = Ecological_condition,
      shape = Soil_habitat
    ),
    colour = "black",
    size = 3.2,
    stroke = 0.7
  ) +
  facet_wrap(
    ~Ecological_condition,
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Natural" = "#A9C9C5",
      "BESE"    = "#FF6B6B",
      "BARE"    = "#F4A340"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "Bulk soil" = 21,
      "Rhizosphere" = 24
    )
  ) +
  labs(
    x = paste0("PCoA1 (", fungal_axis_1, "%)"),
    y = paste0("PCoA2 (", fungal_axis_2, "%)"),
    shape = "Soil habitat"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.position = "bottom",
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black")
  )

p_paired_fun_condition




