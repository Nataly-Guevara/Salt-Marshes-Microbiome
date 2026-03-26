## ============================================================
## 16S ALPHA DIVERSITY — FINAL SCRIPT (LMM + Site)
## Models:
##   A) Status × Compartment (NATURAL vs MANAGED) within BS/RZ
##   B) Restoration trajectory in RZ (NAT vs BESE vs BARE)
##   C) Compartment effect (BS vs RZ) across ALL treatments (incl DEGRADED)
##
## What this script does:
##   1) Loads a phyloseq object
##   2) Filters samples by min read depth
##   3) Computes alpha diversity (Chao1, Shannon, Simpson)
##   4) Fits LMMs with Site as a random effect: (1|Site)
##   5) Runs emmeans + pairwise tests
##   6) Generates boxplots WITH brackets + asterisks
##   7) Prints plots to RStudio and saves PNGs
##   8) Exports summary tables + ANOVA + pairwise results to Excel
##
## Output:
##   - Figures: SUPP_Figures_16S_SITE_LMM/*.png
##   - Tables : SUPP_16S_Alpha_SITE_LMM.xlsx
## ============================================================

rm(list = ls())
graphics.off()
set.seed(1)

## -------------------------
## 0) Libraries
## -------------------------
suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(lme4)
  library(lmerTest)   # provides p-values in anova() for lmer
  library(emmeans)
  library(openxlsx)
})

## -------------------------
## 1) Paths (FIXED)
## -------------------------
base_dir <- "~/Library/CloudStorage/OneDrive-DeakinUniversity/Desktop/ALL DESKTOP/DeSeq2_16S and ITS"
phy_file <- file.path(base_dir, "phyloseq_16S_clean_recodeNATURAL.rds")

cat("Reading:", phy_file, "\n")
stopifnot(file.exists(phy_file))

out_fig_dir <- "SUPP_Figures_16S_SITE_LV"
dir.create(out_fig_dir, showWarnings = FALSE, recursive = TRUE)

out_xlsx <- "SUPP_16S_Alpha_SITE_LV.xlsx"

## -------------------------
## 2) Settings
## -------------------------
metrics   <- c("Chao1","Shannon","Simpson")
min_depth <- 10000

## -------------------------
## 3) Helper functions
## -------------------------

## (A) Fix SampleID format from estimate_richness() if R added leading "X"
fix_alpha_ids <- function(x) {
  x <- as.character(x)
  sub("^X(?=\\d)", "", x, perl = TRUE)
}

## (B) Recode EcologicalStatus so DEGREST etc becomes MANAGED
recode_status_16S <- function(x){
  x <- toupper(trimws(as.character(x)))
  dplyr::case_when(
    x %in% c("NATURAL","NAT") ~ "NATURAL",
    x %in% c("DEGREST","MANAGED","DEGRADED","RESTORED") ~ "MANAGED",
    TRUE ~ x
  )
}

## (C) Clean factor helper
to_factor_clean <- function(x){
  factor(toupper(trimws(as.character(x))))
}

## (D) Make a safe ANOVA table from lmerTest::anova()
##     Different anova outputs can have different column names;
##     we standardize key columns.
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

## (E) Summary table: mean ± SD + n
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

## (F) Convert emmeans pairs() table into a ggpubr bracket table
##     Adds group1/group2 and asterisk codes.
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
  
  ## Stacked y positions so brackets don't overlap
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
## 5) Build alpha diversity table (MAIN)
## -------------------------
get_alpha <- function(phy_obj){
  meta <- as(sample_data(phy_obj), "data.frame") %>%
    mutate(
      SampleID = as.character(SampleID %||% rownames(.)),
      Site = as.character(Site),
      SoilPortion = as.character(SoilPortion),
      Treatment = as.character(Treatment),
      EcologicalStatus = as.character(EcologicalStatus)
    )
  
  alpha <- estimate_richness(phy_obj, measures = metrics) %>%
    rownames_to_column("AlphaID_raw") %>%
    mutate(SampleID = fix_alpha_ids(AlphaID_raw)) %>%
    select(SampleID, all_of(metrics))
  
  out <- alpha %>%
    left_join(meta %>% select(SampleID, Site, SoilPortion, Treatment, EcologicalStatus), by = "SampleID") %>%
    mutate(
      ## Site
      Site = to_factor_clean(Site),
      
      ## Compartment (BS vs RZ)
      SoilPortion = toupper(trimws(SoilPortion)),
      SoilPortion = case_when(
        SoilPortion %in% c("RZ","RHIZOSPHERE") ~ "RZ",
        SoilPortion %in% c("BS","BULK","BULK_SOIL","BULKSOIL","BULK SOIL") ~ "BS",
        TRUE ~ SoilPortion
      ),
      SoilPortion = factor(SoilPortion, levels = c("BS","RZ")),
      
      ## Treatment (kept for Model B only)
      Treatment = toupper(trimws(Treatment)),
      Treatment = case_when(Treatment == "NATURAL" ~ "NAT", TRUE ~ Treatment),
      Treatment = factor(Treatment, levels = c("NAT","BESE","BARE","DEGRADED")),
      
      ## EcologicalStatus -> NATURAL vs MANAGED
      EcologicalStatus = recode_status_16S(EcologicalStatus),
      EcologicalStatus = factor(EcologicalStatus, levels = c("NATURAL","MANAGED"))
    ) %>%
    filter(!is.na(Site), !is.na(SoilPortion))
  
  out
}

alpha_main <- get_alpha(phy_filt)

cat("\nMAIN counts (Status x Portion):\n")
print(table(alpha_main$EcologicalStatus, alpha_main$SoilPortion, useNA="ifany"))

cat("\nMAIN counts (Treatment x Portion):\n")
print(table(alpha_main$Treatment, alpha_main$SoilPortion, useNA="ifany"))

## -------------------------
## 6) Define the three LMMs (Site as random effect)
## -------------------------

## Model A:
## Question: Does diversity differ by EcologicalStatus, Compartment, and their interaction,
## while accounting for Site-level clustering?
fit_A <- function(df, metric){
  dfA <- df %>%
    filter(!is.na(EcologicalStatus), SoilPortion %in% c("BS","RZ")) %>%
    droplevels()
  
  lmer(as.formula(paste0(metric, " ~ EcologicalStatus * SoilPortion + (1|Site)")), data = dfA)
}

## Model B:
## Question: In rhizosphere only, do treatments differ (NAT vs BESE vs BARE),
## controlling for Site?
fit_B <- function(df, metric){
  dfB <- df %>%
    filter(SoilPortion == "RZ", Treatment %in% c("NAT","BESE","BARE")) %>%
    droplevels()
  
  lmer(as.formula(paste0(metric, " ~ Treatment + (1|Site)")), data = dfB)
}

## Model C:
## Question: Overall compartment effect (BS vs RZ) across ALL treatments,
## including DEGRADED (which only contributes BS). NO Treatment term by design.
fit_C <- function(df, metric){
  dfC <- df %>%
    filter(SoilPortion %in% c("BS","RZ")) %>%
    droplevels()
  
  lmer(as.formula(paste0(metric, " ~ SoilPortion + (1|Site)")), data = dfC)
}

## -------------------------
## 7) Plotters (with brackets + asterisks)
## -------------------------

## Plot A: Facet by compartment, compare Status within each compartment
plot_A <- function(df, metric, mod, tag){
  dfA <- df %>% filter(!is.na(EcologicalStatus), SoilPortion %in% c("BS","RZ")) %>% droplevels()
  
  ## emmeans: Status within each compartment
  emA <- emmeans(mod, ~ EcologicalStatus | SoilPortion)
  prA <- as.data.frame(pairs(emA, adjust = "tukey")) %>%
    rename(p.value = p.value)
  
  ## bracket table, keep SoilPortion column for facet placement
  y_base <- max(dfA[[metric]], na.rm = TRUE) * 1.10
  ggpr <- pairs_to_ggpubr(prA, y_base = y_base, step_frac = 0.08)
  
  p <- ggplot(dfA, aes(x = EcologicalStatus, y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    facet_wrap(~ SoilPortion, nrow = 1) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("16S Model A (", tag, "): ", metric),
      x = "Ecological status", y = metric
    )
  
  if(nrow(ggpr) > 0 && "SoilPortion" %in% colnames(prA)){
    ggpr$SoilPortion <- prA$SoilPortion
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

## Plot B: RZ only, pairwise among NAT/BESE/BARE
plot_B <- function(df, metric, mod, tag){
  dfB <- df %>% filter(SoilPortion == "RZ", Treatment %in% c("NAT","BESE","BARE")) %>% droplevels()
  
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
      title = paste0("16S Model B (", tag, "): ", metric, " (RZ only)"),
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

## Plot C: BS vs RZ only (across all treatments), one comparison
plot_C <- function(df, metric, mod, tag){
  dfC <- df %>% filter(SoilPortion %in% c("BS","RZ")) %>% droplevels()
  
  emC <- emmeans(mod, ~ SoilPortion)
  prC <- as.data.frame(pairs(emC, adjust = "none")) %>%
    rename(p.value = p.value)
  
  y_base <- max(dfC[[metric]], na.rm = TRUE) * 1.10
  ggpr <- pairs_to_ggpubr(prC, y_base = y_base, step_frac = 0.08)
  
  p <- ggplot(dfC, aes(x = SoilPortion, y = .data[[metric]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.6) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("16S Model C (", tag, "): ", metric, " (BS vs RZ; Site random effect)"),
      x = "Soil compartment", y = metric
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
## 8) Run ALL models, save plots + tables
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
  
  emA <- emmeans(modA, ~ EcologicalStatus | SoilPortion)
  prA <- as.data.frame(pairs(emA, adjust = "tukey")) %>%
    rename(p.value = p.value) %>%
    mutate(Model = "A", Metric = metric, ContrastSet = "Status|Portion")
  
  sumA <- summarise_metric(
    alpha_main %>% filter(!is.na(EcologicalStatus), SoilPortion %in% c("BS","RZ")),
    c("EcologicalStatus","SoilPortion"), metric, "A"
  )
  
  pA <- plot_A(alpha_main, metric, modA, "MAIN")
  print(pA)  # IMPORTANT: ensures plot shows in RStudio
  ggsave(file.path(out_fig_dir, paste0("SUPP_16S_ModelA_MAIN_", metric, "_LMM.png")),
         pA, width = 8, height = 4.8, dpi = 300)
  
  ## =====================
  ## MODEL B
  ## =====================
  modB <- fit_B(alpha_main, metric)
  anB  <- safe_anova_df(modB, "B", metric)
  
  emB <- emmeans(modB, ~ Treatment)
  prB <- as.data.frame(pairs(emB, adjust = "tukey")) %>%
    rename(p.value = p.value) %>%
    mutate(Model = "B", Metric = metric, ContrastSet = "Treatment")
  
  sumB <- summarise_metric(
    alpha_main %>% filter(SoilPortion == "RZ", Treatment %in% c("NAT","BESE","BARE")),
    c("Treatment"), metric, "B"
  )
  
  pB <- plot_B(alpha_main, metric, modB, "MAIN")
  print(pB)
  ggsave(file.path(out_fig_dir, paste0("SUPP_16S_ModelB_MAIN_", metric, "_LMM.png")),
         pB, width = 6.2, height = 5.0, dpi = 300)
  
  ## =====================
  ## MODEL C  (includes DEGRADED; no Treatment term)
  ## =====================
  modC <- fit_C(alpha_main, metric)
  anC  <- safe_anova_df(modC, "C", metric)
  
  emC <- emmeans(modC, ~ SoilPortion)
  prC <- as.data.frame(pairs(emC, adjust = "none")) %>%
    rename(p.value = p.value) %>%
    mutate(Model = "C", Metric = metric, ContrastSet = "Portion")
  
  sumC <- summarise_metric(
    alpha_main %>% filter(SoilPortion %in% c("BS","RZ")),
    c("SoilPortion"), metric, "C"
  )
  
  pC <- plot_C(alpha_main, metric, modC, "MAIN")
  print(pC)
  ggsave(file.path(out_fig_dir, paste0("SUPP_16S_ModelC_MAIN_", metric, "_LMM.png")),
         pC, width = 6.2, height = 5.0, dpi = 300)
  
  ## =====================
  ## Collect across models
  ## =====================
  ALL_ANOVA[[metric]] <- bind_rows(anA, anB, anC)
  ALL_PAIRS[[metric]] <- bind_rows(prA, prB, prC)
  ALL_SUMM[[metric]]  <- bind_rows(sumA, sumB, sumC)
}

ANOVA_ALL <- bind_rows(ALL_ANOVA)
PAIRS_ALL <- bind_rows(ALL_PAIRS)
SUMM_ALL  <- bind_rows(ALL_SUMM)

## -------------------------
## 9) Export to Excel (clean sheet names)
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

## Optional: per-metric sheets
for(metric in metrics){
  add_sheet(wb, paste0("SUMM_", metric),  SUMM_ALL  %>% filter(Metric == metric))
  add_sheet(wb, paste0("ANOVA_", metric), ANOVA_ALL %>% filter(Metric == metric))
  add_sheet(wb, paste0("PAIR_", metric),  PAIRS_ALL %>% filter(Metric == metric))
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("\nDONE.\n")
cat("Plots saved in: ", normalizePath(out_fig_dir), "\n")
cat("Tables saved as: ", normalizePath(out_xlsx), "\n")
