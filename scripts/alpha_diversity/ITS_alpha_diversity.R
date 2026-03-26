## ============================================================
## ITS ALPHA DIVERSITY — SUPPLEMENTARY SCRIPT (LMM + Site)
## Models:
##   A) StatusGroup × Compartment (NATURAL vs MANAGED) within BS/RZ
##   B) Restoration trajectory in RZ (NATURAL vs BESE vs BARE)
##   C) Compartment effect (BS vs RZ) across ALL treatments (incl DEGRADED)
##
## Output:
##   - Figures: SUPP_Figures_ITS_SITE_LMM/*.png
##   - Tables : SUPP_ITS_Alpha_SITE_LMM.xlsx
## ============================================================

rm(list = ls())
graphics.off()
set.seed(1)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(openxlsx)
})

## -------------------------
## 1) Paths (FIXED)
## -------------------------
base_dir <- "~/Library/CloudStorage/OneDrive-DeakinUniversity/Desktop/ALL DESKTOP/Phyloseq_ITS"

## CHANGE IF NEEDED:
phy_file <- file.path(base_dir, "phyloseq_ITS_fixed.rds")

cat("Reading:", phy_file, "\n")
stopifnot(file.exists(phy_file))

out_fig_dir <- file.path(base_dir, "SUPP_Figures_ITS_SITE_LMM")
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

out_xlsx <- file.path(base_dir, "SUPP_ITS_Alpha_SITE_LMM.xlsx")

## -------------------------
## 2) Settings
## -------------------------
metrics   <- c("Chao1","Shannon","Simpson")
min_depth <- 5000   # ITS often lower depth than 16S; adjust if you want

## -------------------------
## 3) Helper functions
## -------------------------

# Fix SampleID format from estimate_richness() if R added leading "X"
fix_alpha_ids <- function(x) {
  x <- as.character(x)
  sub("^X(?=\\d)", "", x, perl = TRUE)
}

to_factor_clean <- function(x){
  factor(toupper(trimws(as.character(x))))
}

# Standardise SoilPortion -> Compartment (BS/RZ)
std_compartment <- function(x){
  sp <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    sp %in% c("RZ","RHIZOSPHERE") ~ "RZ",
    sp %in% c("BS","BULK","BULK_SOIL","BULKSOIL","BULK SOIL") ~ "BS",
    TRUE ~ NA_character_
  )
}

# Standardise Treatment values (keeps NATURAL/BESE/BARE/DEGRADED)
std_treatment <- function(x){
  tr <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    tr %in% c("NAT","NATURAL") ~ "NATURAL",
    tr %in% c("DEGREST","DEGRADED") ~ "DEGRADED",
    tr %in% c("BESE") ~ "BESE",
    tr %in% c("BARE") ~ "BARE",
    TRUE ~ tr
  )
}

# StatusGroup = NATURAL vs MANAGED (everything not NATURAL -> MANAGED)
make_statusgroup <- function(tr){
  tr <- toupper(trimws(as.character(tr)))
  ifelse(tr == "NATURAL", "NATURAL", "MANAGED")
}

# Safe ANOVA table from lmerTest::anova()
safe_anova_df <- function(mod, model_tag, metric){
  df <- as.data.frame(anova(mod)) %>% rownames_to_column("Term")
  nms <- names(df)
  
  if("F value" %in% nms)    df <- dplyr::rename(df, F.value = `F value`)
  if("Pr(>F)" %in% nms)     df <- dplyr::rename(df, p.value = `Pr(>F)`)
  if("Pr(>Chisq)" %in% nms) df <- dplyr::rename(df, p.value = `Pr(>Chisq)`)
  
  df %>%
    mutate(Model = model_tag, Metric = metric) %>%
    select(any_of(c("Model","Metric","Term","NumDF","DenDF","Df","Chisq","F.value","p.value")))
}

# mean ± SD + n
summarise_metric <- function(df, group_vars, metric, model_tag){
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n = n(),
      mean = mean(.data[[metric]], na.rm = TRUE),
      sd   = sd(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Model = model_tag, Metric = metric) %>%
    relocate(Model, Metric)
}

# pairs() -> ggpubr bracket table with asterisks
pairs_to_ggpubr <- function(pairs_df, y_base, step_frac = 0.06){
  if(nrow(pairs_df) == 0) return(pairs_df)
  
  out <- pairs_df %>%
    mutate(
      contrast = as.character(contrast),
      group1 = str_trim(str_split_fixed(contrast, "-", 2)[,1]),
      group2 = str_trim(str_split_fixed(contrast, "-", 2)[,2]),
      p = p.value
    ) %>%
    mutate(
      p.adj.signif = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE ~ "ns"
      )
    )
  
  out$y.position <- y_base + (seq_len(nrow(out)) - 1) * (step_frac * y_base)
  out
}

## -------------------------
## 4) Load phyloseq + filter by depth
## -------------------------
phy <- readRDS(phy_file)

phy_filt <- prune_samples(sample_sums(phy) >= min_depth, phy)

cat("Samples before:", nsamples(phy), "\n")
cat("Samples after :", nsamples(phy_filt), "\n")
cat("Min reads after filter:", min(sample_sums(phy_filt)), "\n")

## -------------------------
## 5) Build alpha diversity table
## -------------------------
get_alpha <- function(phy_obj){
  
  md <- as(sample_data(phy_obj), "data.frame")
  
  # ensure SampleID exists
  if(!"SampleID" %in% names(md)) md$SampleID <- rownames(md)
  md$SampleID <- as.character(md$SampleID)
  
  # standardize key columns (if they exist)
  if(!"Site" %in% names(md)) stop("sample_data must include a 'Site' column.")
  if(!"SoilPortion" %in% names(md)) stop("sample_data must include a 'SoilPortion' column.")
  if(!"Treatment" %in% names(md)) stop("sample_data must include a 'Treatment' column.")
  
  alpha <- estimate_richness(phy_obj, measures = metrics) %>%
    rownames_to_column("AlphaID_raw") %>%
    mutate(SampleID = fix_alpha_ids(AlphaID_raw)) %>%
    select(SampleID, all_of(metrics))
  
  out <- alpha %>%
    left_join(
      md %>% select(SampleID, Site, SoilPortion, Treatment),
      by = "SampleID"
    ) %>%
    mutate(
      Site = to_factor_clean(Site),
      Compartment = factor(std_compartment(SoilPortion), levels = c("BS","RZ")),
      Treatment = factor(std_treatment(Treatment), levels = c("NATURAL","BESE","BARE","DEGRADED")),
      StatusGroup = factor(make_statusgroup(Treatment), levels = c("NATURAL","MANAGED"))
    ) %>%
    filter(!is.na(Site), !is.na(Compartment), !is.na(Treatment))
  
  out
}

alpha_main <- get_alpha(phy_filt)

cat("\nCounts (StatusGroup x Compartment):\n")
print(table(alpha_main$StatusGroup, alpha_main$Compartment, useNA="ifany"))

cat("\nCounts (Treatment x Compartment):\n")
print(table(alpha_main$Treatment, alpha_main$Compartment, useNA="ifany"))

## -------------------------
## 6) Define the three LMMs (Site random effect)
## -------------------------

# Model A: StatusGroup × Compartment
fit_A <- function(df, metric){
  dfA <- df %>% filter(Compartment %in% c("BS","RZ")) %>% droplevels()
  lmer(as.formula(paste0(metric, " ~ StatusGroup * Compartment + (1|Site)")), data = dfA)
}

# Model B: RZ only, NATURAL vs BESE vs BARE
fit_B <- function(df, metric){
  dfB <- df %>% filter(Compartment == "RZ", Treatment %in% c("NATURAL","BESE","BARE")) %>% droplevels()
  lmer(as.formula(paste0(metric, " ~ Treatment + (1|Site)")), data = dfB)
}

# Model C: Compartment only across all treatments (incl DEGRADED)
fit_C <- function(df, metric){
  dfC <- df %>% filter(Compartment %in% c("BS","RZ")) %>% droplevels()
  lmer(as.formula(paste0(metric, " ~ Compartment + (1|Site)")), data = dfC)
}

## -------------------------
## 7) Plotters (with brackets + asterisks)
## -------------------------

plot_A <- function(df, metric, mod, tag){
  dfA <- df %>% filter(Compartment %in% c("BS","RZ")) %>% droplevels()
  
  emA <- emmeans(mod, ~ StatusGroup | Compartment)
  prA <- as.data.frame(pairs(emA, adjust = "tukey")) %>%
    rename(p.value = p.value)
  
  y_base <- max(dfA[[metric]], na.rm = TRUE) * 1.10
  ggpr <- pairs_to_ggpubr(prA, y_base = y_base, step_frac = 0.08)
  
  # IMPORTANT for faceting: carry Compartment into bracket table
  if(nrow(ggpr) > 0 && "Compartment" %in% names(prA)){
    ggpr$Compartment <- prA$Compartment
  }
  
  p <- ggplot(dfA, aes(x = StatusGroup, y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    facet_wrap(~ Compartment, nrow = 1) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("ITS Model A (", tag, "): ", metric),
      x = "Status group", y = metric
    )
  
  if(nrow(ggpr) > 0 && "Compartment" %in% names(ggpr)){
    p <- p + stat_pvalue_manual(
      ggpr,
      label = "p.adj.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01
    )
  }
  
  p
}

plot_B <- function(df, metric, mod, tag){
  dfB <- df %>% filter(Compartment == "RZ", Treatment %in% c("NATURAL","BESE","BARE")) %>% droplevels()
  
  emB <- emmeans(mod, ~ Treatment)
  prB <- as.data.frame(pairs(emB, adjust = "tukey")) %>%
    rename(p.value = p.value)
  
  y_base <- max(dfB[[metric]], na.rm = TRUE) * 1.12
  ggpr <- pairs_to_ggpubr(prB, y_base = y_base, step_frac = 0.06)
  
  p <- ggplot(dfB, aes(x = Treatment, y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("ITS Model B (", tag, "): ", metric, " (RZ only)"),
      x = "Treatment", y = metric
    )
  
  if(nrow(ggpr) > 0){
    p <- p + stat_pvalue_manual(
      ggpr,
      label = "p.adj.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01
    )
  }
  
  p
}

plot_C <- function(df, metric, mod, tag){
  dfC <- df %>% filter(Compartment %in% c("BS","RZ")) %>% droplevels()
  
  emC <- emmeans(mod, ~ Compartment)
  prC <- as.data.frame(pairs(emC, adjust = "none")) %>%
    rename(p.value = p.value)
  
  y_base <- max(dfC[[metric]], na.rm = TRUE) * 1.10
  ggpr <- pairs_to_ggpubr(prC, y_base = y_base, step_frac = 0.08)
  
  p <- ggplot(dfC, aes(x = Compartment, y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("ITS Model C (", tag, "): ", metric, " (BS vs RZ)"),
      x = "Compartment", y = metric
    )
  
  if(nrow(ggpr) > 0){
    p <- p + stat_pvalue_manual(
      ggpr,
      label = "p.adj.signif",
      xmin = "group1", xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01
    )
  }
  
  p
}

## -------------------------
## 8) Run all models, save plots + tables
## -------------------------
ALL_ANOVA <- list()
ALL_PAIRS <- list()
ALL_SUMM  <- list()

for(metric in metrics){
  
  ## =====================
  ## MODEL A
  ## =====================
  modA <- fit_A(alpha_main, metric)
  anA  <- safe_anova_df(modA, "A", metric)
  
  emA <- emmeans(modA, ~ StatusGroup | Compartment)
  prA <- as.data.frame(pairs(emA, adjust = "tukey")) %>%
    rename(p.value = p.value) %>%
    mutate(Model = "A", Metric = metric, ContrastSet = "StatusGroup|Compartment")
  
  sumA <- summarise_metric(
    alpha_main,
    c("StatusGroup","Compartment"), metric, "A"
  )
  
  pA <- plot_A(alpha_main, metric, modA, "MAIN")
  print(pA)
  ggsave(file.path(out_fig_dir, paste0("SUPP_ITS_ModelA_MAIN_", metric, "_LMM.png")),
         pA, width = 8, height = 4.8, dpi = 300)
  
  ## =====================
  ## MODEL B
  ## =====================
  modB <- fit_B(alpha_main, metric)
  anB  <- safe_anova_df(modB, "B", metric)
  
  emB <- emmeans(modB, ~ Treatment)
  prB <- as.data.frame(pairs(emB, adjust = "tukey")) %>%
    rename(p.value = p.value) %>%
    mutate(Model = "B", Metric = metric, ContrastSet = "Treatment_RZ")
  
  sumB <- summarise_metric(
    alpha_main %>% filter(Compartment == "RZ", Treatment %in% c("NATURAL","BESE","BARE")),
    c("Treatment"), metric, "B"
  )
  
  pB <- plot_B(alpha_main, metric, modB, "MAIN")
  print(pB)
  ggsave(file.path(out_fig_dir, paste0("SUPP_ITS_ModelB_MAIN_", metric, "_LMM.png")),
         pB, width = 6.2, height = 5.0, dpi = 300)
  
  ## =====================
  ## MODEL C
  ## =====================
  modC <- fit_C(alpha_main, metric)
  anC  <- safe_anova_df(modC, "C", metric)
  
  emC <- emmeans(modC, ~ Compartment)
  prC <- as.data.frame(pairs(emC, adjust = "none")) %>%
    rename(p.value = p.value) %>%
    mutate(Model = "C", Metric = metric, ContrastSet = "Compartment")
  
  sumC <- summarise_metric(
    alpha_main,
    c("Compartment"), metric, "C"
  )
  
  pC <- plot_C(alpha_main, metric, modC, "MAIN")
  print(pC)
  ggsave(file.path(out_fig_dir, paste0("SUPP_ITS_ModelC_MAIN_", metric, "_LMM.png")),
         pC, width = 6.2, height = 5.0, dpi = 300)
  
  ## Collect
  ALL_ANOVA[[metric]] <- bind_rows(anA, anB, anC)
  ALL_PAIRS[[metric]] <- bind_rows(prA, prB, prC)
  ALL_SUMM[[metric]]  <- bind_rows(sumA, sumB, sumC)
}

ANOVA_ALL <- bind_rows(ALL_ANOVA)
PAIRS_ALL <- bind_rows(ALL_PAIRS)
SUMM_ALL  <- bind_rows(ALL_SUMM)

## -------------------------
## 9) Export to Excel
## -------------------------
wb <- createWorkbook()

add_sheet <- function(wb, name, df){
  clean <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", name)
  clean <- substr(clean, 1, 31)
  if(clean %in% names(wb)){
    base <- substr(clean, 1, 28)
    i <- 1
    while(paste0(base,"_",i) %in% names(wb)) i <- i + 1
    clean <- paste0(base,"_",i)
  }
  addWorksheet(wb, clean)
  writeData(wb, clean, df)
}

add_sheet(wb, "SUMMARIES_ALL", SUMM_ALL)
add_sheet(wb, "ANOVA_ALL",     ANOVA_ALL)
add_sheet(wb, "PAIRWISE_ALL",  PAIRS_ALL)

for(metric in metrics){
  add_sheet(wb, paste0("SUMM_", metric),  SUMM_ALL  %>% filter(Metric == metric))
  add_sheet(wb, paste0("ANOVA_", metric), ANOVA_ALL %>% filter(Metric == metric))
  add_sheet(wb, paste0("PAIR_", metric),  PAIRS_ALL %>% filter(Metric == metric))
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\nDONE.\n")
cat("Plots saved in: ", normalizePath(out_fig_dir), "\n")
cat("Tables saved as: ", normalizePath(out_xlsx), "\n")

