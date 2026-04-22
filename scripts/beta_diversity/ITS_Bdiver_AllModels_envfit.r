## =========================================================
## ITS Model A beta diversity (FINAL CLEAN)
## NATURAL vs MANAGED x Compartment (BS/RZ), conditioned on Site
## Bray (relative abundance) + db-RDA + PERMANOVA + dispersion + envfit (FDR)
## Blocked permutations by Site
## Outputs: plot + CSV tables for Supplementary
## =========================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(permute)
  library(grid)
})

set.seed(1)

## ---------- USER SETTINGS ----------
in_rds  <- "~/Library/CloudStorage/OneDrive-DeakinUniversity/Desktop/ALL DESKTOP/Phyloseq_ITS/phyloseq_ITS_fixed.rds"
out_dir <- "supp_beta_ITS_ModelA"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Helpers ----------
relab_comm <- function(ps){
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  mat <- as(otu_table(ps_rel), "matrix")
  if (taxa_are_rows(ps_rel)) mat <- t(mat)   # samples x taxa
  mat
}

to_num_if_possible <- function(x){
  if(is.numeric(x) || is.integer(x)) return(as.numeric(x))
  if(is.factor(x)) x <- as.character(x)
  if(is.character(x)){
    x2 <- suppressWarnings(as.numeric(x))
    if(mean(!is.na(x2)) >= 0.7) return(x2)
  }
  x
}

fmt_stats <- function(Fval, pval, adjR2){
  p_txt <- ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))
  paste0("F = ", sprintf("%.2f", Fval),
         ", p = ", p_txt,
         ", Adj R² = ", sprintf("%.3f", adjR2))
}

make_env_numeric <- function(ps, meta, drop_like){
  md <- as(sample_data(ps), "data.frame")
  
  # ensure SampleID exists (yours does)
  if(!"SampleID" %in% names(md)) md$SampleID <- rownames(md)
  
  # drop non-env columns
  md2 <- md[, !tolower(names(md)) %in% tolower(drop_like), drop = FALSE]
  
  md3 <- as.data.frame(lapply(md2, to_num_if_possible))
  md3$SampleID <- md$SampleID
  
  env_numeric <- md3 %>%
    dplyr::select(where(is.numeric)) %>%
    mutate(SampleID = md$SampleID)
  
  # align to meta order via SampleID
  env_numeric <- env_numeric %>%
    dplyr::semi_join(meta %>% dplyr::select(SampleID), by = "SampleID") %>%
    dplyr::right_join(meta %>% dplyr::select(SampleID), by = "SampleID") %>%
    tibble::column_to_rownames("SampleID")
  
  # drop all-NA columns
  keep <- colSums(!is.na(env_numeric)) > 0
  env_numeric[, keep, drop = FALSE]
}

## ---------- Load ----------
ps <- readRDS(in_rds)

## =========================================
## 1) SUBSET + recode Compartment + StatusGroup
## =========================================
ps_A <- subset_samples(ps, !is.na(Site) & !is.na(SoilPortion))
ps_A <- prune_samples(sample_sums(ps_A) > 0, ps_A)
ps_A <- prune_taxa(taxa_sums(ps_A) > 0, ps_A)

meta_A <- as(sample_data(ps_A), "data.frame")
if(!"SampleID" %in% colnames(meta_A)) meta_A$SampleID <- rownames(meta_A)

# Compartment recode
sp_raw <- toupper(trimws(as.character(meta_A$SoilPortion)))
meta_A$Compartment <- dplyr::case_when(
  sp_raw %in% c("RZ","RHIZOSPHERE") ~ "RZ",
  sp_raw %in% c("BS","BULK","BULK_SOIL","BULKSOIL","BULK SOIL") ~ "BS",
  TRUE ~ NA_character_
)
meta_A$Compartment <- factor(meta_A$Compartment, levels = c("BS","RZ"))

# keep only BS/RZ
keep_ids <- meta_A$SampleID[!is.na(meta_A$Compartment)]
ps_A <- prune_samples(keep_ids, ps_A)
ps_A <- prune_taxa(taxa_sums(ps_A) > 0, ps_A)

# refresh meta
meta_A <- as(sample_data(ps_A), "data.frame")
if(!"SampleID" %in% colnames(meta_A)) meta_A$SampleID <- rownames(meta_A)

sp_raw <- toupper(trimws(as.character(meta_A$SoilPortion)))
meta_A$Compartment <- dplyr::case_when(
  sp_raw %in% c("RZ","RHIZOSPHERE") ~ "RZ",
  sp_raw %in% c("BS","BULK","BULK_SOIL","BULKSOIL","BULK SOIL") ~ "BS",
  TRUE ~ NA_character_
)
meta_A$Compartment <- factor(meta_A$Compartment, levels = c("BS","RZ"))

# StatusGroup
if("EcologicalStatus" %in% colnames(meta_A)){
  es <- toupper(trimws(as.character(meta_A$EcologicalStatus)))
  meta_A$StatusGroup <- ifelse(es == "NATURAL", "NATURAL", "MANAGED")
} else if("Treatment" %in% colnames(meta_A)){
  tr <- toupper(trimws(as.character(meta_A$Treatment)))
  meta_A$StatusGroup <- ifelse(tr == "NATURAL", "NATURAL", "MANAGED")
} else {
  stop("No EcologicalStatus or Treatment found in sample_data().")
}
meta_A$StatusGroup <- factor(meta_A$StatusGroup, levels = c("NATURAL","MANAGED"))
meta_A$Site <- factor(meta_A$Site)

cat("\nSample counts by StatusGroup x Compartment x Site:\n")
print(table(meta_A$StatusGroup, meta_A$Compartment, meta_A$Site))

## =========================================
## 2) COMMUNITY MATRIX + align meta
## =========================================
comm_A <- relab_comm(ps_A)

# Ensure rownames(comm_A) are SampleIDs
# phyloseq uses sample_names(ps_A) as rownames of comm_A
stopifnot(identical(rownames(comm_A), sample_names(ps_A)))

# Align meta to comm order by SampleID
meta_A <- meta_A %>%
  dplyr::mutate(SampleID = as.character(SampleID)) %>%
  dplyr::slice(match(rownames(comm_A), SampleID))

stopifnot(identical(meta_A$SampleID, rownames(comm_A)))

## =========================================
## 3) Blocked permutations by Site
## =========================================
ctrl_A <- how(nperm = 999, blocks = meta_A$Site)

## =========================================
## 4) db-RDA (capscale) + ANOVAs
## =========================================
cap_A <- capscale(comm_A ~ StatusGroup * Compartment + Condition(Site),
                  distance = "bray", data = meta_A)

a_all   <- anova(cap_A, permutations = ctrl_A)
a_terms <- anova(cap_A, by = "terms",  permutations = ctrl_A)
a_marg  <- anova(cap_A, by = "margin", permutations = ctrl_A)
a_axis  <- anova(cap_A, by = "axis",   permutations = ctrl_A)

adjR2 <- RsquareAdj(cap_A)$adj.r.squared

write.csv(as.data.frame(a_all),   file.path(out_dir, "ITS_ModelA_capscale_ANOVA_overall.csv"))
write.csv(as.data.frame(a_terms), file.path(out_dir, "ITS_ModelA_capscale_ANOVA_terms_sequential.csv"))
write.csv(as.data.frame(a_marg),  file.path(out_dir, "ITS_ModelA_capscale_ANOVA_terms_marginal.csv"))
write.csv(as.data.frame(a_axis),  file.path(out_dir, "ITS_ModelA_capscale_ANOVA_axis.csv"))
write.csv(data.frame(AdjR2 = adjR2), file.path(out_dir, "ITS_ModelA_capscale_AdjR2.csv"), row.names = FALSE)

term_int <- "StatusGroup:Compartment"
if(term_int %in% rownames(a_marg)){
  subtitle_txt <- paste0("Interaction ", fmt_stats(a_marg[term_int,"F"], a_marg[term_int,"Pr(>F)"], adjR2))
} else {
  subtitle_txt <- fmt_stats(a_all[1,"F"], a_all[1,"Pr(>F)"], adjR2)
}

## =========================================
## 5) PERMANOVA + dispersion
## =========================================
dist_A <- vegdist(comm_A, method = "bray")

permA <- adonis2(dist_A ~ StatusGroup * Compartment, data = meta_A,
                 permutations = 999, strata = meta_A$Site, by = "margin")
write.csv(as.data.frame(permA), file.path(out_dir, "ITS_ModelA_PERMANOVA_adonis2_margin.csv"))

disp_Status <- betadisper(dist_A, meta_A$StatusGroup)
disp_Comp   <- betadisper(dist_A, meta_A$Compartment)
disp_Int    <- betadisper(dist_A, interaction(meta_A$StatusGroup, meta_A$Compartment))

disp_tab <- bind_rows(
  tibble(Test = "betadisper_StatusGroup",
         F = permutest(disp_Status, permutations = ctrl_A)$tab[1,"F"],
         p = permutest(disp_Status, permutations = ctrl_A)$tab[1,"Pr(>F)"]),
  tibble(Test = "betadisper_Compartment",
         F = permutest(disp_Comp, permutations = ctrl_A)$tab[1,"F"],
         p = permutest(disp_Comp, permutations = ctrl_A)$tab[1,"Pr(>F)"]),
  tibble(Test = "betadisper_StatusGroup_x_Compartment",
         F = permutest(disp_Int, permutations = ctrl_A)$tab[1,"F"],
         p = permutest(disp_Int, permutations = ctrl_A)$tab[1,"Pr(>F)"])
)
write.csv(disp_tab, file.path(out_dir, "ITS_ModelA_Dispersion_betadisper_permutest_blocked.csv"), row.names = FALSE)

## =========================================
## 5b) Variance partitioning (varpart)
## =========================================
X_site   <- data.frame(Site = meta_A$Site)
X_status <- data.frame(StatusGroup = meta_A$StatusGroup)
X_comp   <- data.frame(Compartment = meta_A$Compartment)

vp_A <- varpart(dist_A, X_site, X_status, X_comp)

cat("\n=== ITS Model A varpart ===\n")
print(vp_A)

vp_tab_A <- as.data.frame(vp_A$part$indfract)
vp_tab_A$Fraction <- rownames(vp_tab_A)
rownames(vp_tab_A) <- NULL

# grab adjusted column robustly
adj_col <- grep("Adj", names(vp_tab_A), value = TRUE)[1]
vp_tab_A <- vp_tab_A %>%
  dplyr::rename(AdjR2 = !!adj_col) %>%
  dplyr::mutate(Percent = 100 * AdjR2)

labels <- c(
  "[a]"="Unique Site",
  "[b]"="Unique StatusGroup",
  "[c]"="Unique Compartment",
  "[d]"="Shared Site∩Status",
  "[e]"="Shared Site∩Compartment",
  "[f]"="Shared Status∩Compartment",
  "[g]"="Shared Site∩Status∩Compartment"
)

vp_tab_A$Meaning <- dplyr::recode(vp_tab_A$Fraction, !!!labels)

write.csv(vp_tab_A,
          file.path(out_dir, "ITS_ModelA_VariancePartition_varpart.csv"),
          row.names = FALSE)

cat("\nVariance partitioning (AdjR2 and %):\n")
print(vp_tab_A)


## =========================================
## 6) envfit (ALL numeric env vars) + FDR
## =========================================
drop_like <- c("SampleID","sampleid","Treatment","EcologicalStatus","StatusGroup",
               "SoilPortion","Compartment","Site","Block","Plot")

env_numeric <- make_env_numeric(ps_A, meta_A, drop_like)

set.seed(123)
fit_all <- envfit(cap_A, env_numeric, permutations = ctrl_A, na.rm = TRUE)

envfit_table_all <- data.frame(
  var = rownames(fit_all$vectors$arrows),
  CAP1 = fit_all$vectors$arrows[,1],
  CAP2 = fit_all$vectors$arrows[,2],
  r2  = fit_all$vectors$r,
  p   = fit_all$vectors$pvals,
  stringsAsFactors = FALSE
) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p, desc(r2))

write.csv(envfit_table_all, file.path(out_dir, "ITS_ModelA_envfit_ALL_numeric_FDR.csv"), row.names = FALSE)

df_vec <- envfit_table_all %>% filter(p_adj < 0.05)

## =========================================
## 7) Plot (95% ellipses + % axis + envfit arrows)
## =========================================
scr <- vegan::scores(cap_A, display = "sites", choices = 1:2)
df_sites <- tibble::as_tibble(scr, rownames = "SampleID") %>%
  left_join(meta_A %>% dplyr::select(SampleID, StatusGroup, Compartment),
            by = "SampleID")

stopifnot(!anyDuplicated(df_sites$SampleID))

eig_constrained <- cap_A$CCA$eig
cap1_pct <- 100 * eig_constrained[1] / sum(eig_constrained)
cap2_pct <- 100 * eig_constrained[2] / sum(eig_constrained)

x_lab <- paste0("CAP1 (", round(cap1_pct, 1), "%)")
y_lab <- paste0("CAP2 (", round(cap2_pct, 1), "%)")

arrow_mult <- 0.9 * max(
  diff(range(df_sites$CAP1, na.rm = TRUE)),
  diff(range(df_sites$CAP2, na.rm = TRUE))
)

col_map   <- c("NATURAL" = "green3", "MANAGED" = "purple3")
shape_map <- c("BS" = 16, "RZ" = 17)

pA <- ggplot(df_sites, aes(CAP1, CAP2, color = StatusGroup, shape = Compartment)) +
  geom_point(size = 2.2, alpha = 0.85) +
  stat_ellipse(aes(group = interaction(StatusGroup, Compartment)),
               type = "t", level = 0.95, linewidth = 0.8) +
  scale_color_manual(values = col_map) +
  scale_shape_manual(values = shape_map) +
  labs(title = "ITS Model A: db-RDA + envfit",
       subtitle = subtitle_txt,
       x = x_lab, y = y_lab) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

if(nrow(df_vec) > 0){
  pA <- pA +
    geom_segment(
      data = df_vec,
      aes(x = 0, y = 0, xend = CAP1 * arrow_mult, yend = CAP2 * arrow_mult),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.02, "npc")),
      linewidth = 0.7,
      color = "black"
    ) +
    geom_text(
      data = df_vec,
      aes(x = CAP1 * arrow_mult * 1.07, y = CAP2 * arrow_mult * 1.07, label = var),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    )
}

ggsave(file.path(out_dir, "ITS_ModelA_dbRDA_envfit_ellipses.png"), pA, width = 6.6, height = 5.0, dpi = 600)
ggsave(file.path(out_dir, "ITS_ModelA_dbRDA_envfit_ellipses.pdf"), pA, width = 6.6, height = 5.0, device = "pdf")

cat("\nDONE. Outputs saved to: ", normalizePath(out_dir), "\n")


print(pA)



## =========================================================
## ITS Model B beta diversity (FINAL CLEAN)
## RZ only: NATURAL vs BESE vs BARE (exclude DEGRADED)
## Conditioned on Site + blocked permutations by Site
## Bray (relative abundance) + db-RDA + PERMANOVA + dispersion + envfit (FDR)
## Outputs: plot + CSV tables for Supplementary
## =========================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(permute)
  library(grid)
})

set.seed(1)

## ---------- USER SETTINGS ----------
in_rds  <- "~/Library/CloudStorage/OneDrive-DeakinUniversity/Desktop/ALL DESKTOP/Phyloseq_ITS/phyloseq_ITS_fixed.rds"
out_dir <- "supp_beta_ITS_ModelB"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Helpers ----------
relab_comm <- function(ps){
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  mat <- as(otu_table(ps_rel), "matrix")
  if (taxa_are_rows(ps_rel)) mat <- t(mat)   # samples x taxa
  mat
}

to_num_if_possible <- function(x){
  if(is.numeric(x) || is.integer(x)) return(as.numeric(x))
  if(is.factor(x)) x <- as.character(x)
  if(is.character(x)){
    x2 <- suppressWarnings(as.numeric(x))
    if(mean(!is.na(x2)) >= 0.7) return(x2)
  }
  x
}

fmt_stats <- function(Fval, pval, adjR2){
  p_txt <- ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))
  paste0("F = ", sprintf("%.2f", Fval),
         ", p = ", p_txt,
         ", Adj R² = ", sprintf("%.3f", adjR2))
}

make_env_numeric <- function(ps, meta, drop_like){
  md <- as(sample_data(ps), "data.frame")
  if(!"SampleID" %in% names(md)) md$SampleID <- rownames(md)
  
  md2 <- md[, !tolower(names(md)) %in% tolower(drop_like), drop = FALSE]
  md3 <- as.data.frame(lapply(md2, to_num_if_possible))
  md3$SampleID <- md$SampleID
  
  env_numeric <- md3 %>%
    dplyr::select(where(is.numeric)) %>%
    mutate(SampleID = md$SampleID)
  
  env_numeric <- env_numeric %>%
    dplyr::semi_join(meta %>% dplyr::select(SampleID), by = "SampleID") %>%
    dplyr::right_join(meta %>% dplyr::select(SampleID), by = "SampleID") %>%
    tibble::column_to_rownames("SampleID")
  
  keep <- colSums(!is.na(env_numeric)) > 0
  env_numeric[, keep, drop = FALSE]
}

## ---------- Load ----------
ps <- readRDS(in_rds)

## =========================================
## 1) SUBSET: RZ only + Treatment3 (NATURAL/BESE/BARE) + Site
## =========================================
ps_B <- subset_samples(ps,
                       !is.na(Site) &
                         SoilPortion == "RZ" &
                         !is.na(Treatment) &
                         Treatment %in% c("NATURAL","BESE","BARE"))
ps_B <- prune_samples(sample_sums(ps_B) > 0, ps_B)
ps_B <- prune_taxa(taxa_sums(ps_B) > 0, ps_B)

meta_B <- as(sample_data(ps_B), "data.frame")
if(!"SampleID" %in% colnames(meta_B)) meta_B$SampleID <- rownames(meta_B)

meta_B$Site <- factor(meta_B$Site)
meta_B$Treatment3 <- factor(as.character(meta_B$Treatment),
                            levels = c("NATURAL","BESE","BARE"))

# sanity checks
stopifnot(all(levels(meta_B$Treatment3) == c("NATURAL","BESE","BARE")))
stopifnot(all(unique(as.character(meta_B$SoilPortion)) == "RZ"))
stopifnot(length(unique(meta_B$Site)) > 1)

cat("\nSample counts by Treatment3 x Site:\n")
print(table(meta_B$Treatment3, meta_B$Site))

## =========================================
## 2) COMMUNITY MATRIX + align meta
## =========================================
comm_B <- relab_comm(ps_B)

stopifnot(identical(rownames(comm_B), sample_names(ps_B)))

meta_B <- meta_B %>%
  dplyr::mutate(SampleID = as.character(SampleID)) %>%
  dplyr::slice(match(rownames(comm_B), SampleID))

stopifnot(identical(meta_B$SampleID, rownames(comm_B)))

## =========================================
## 3) Blocked permutations by Site
## =========================================
ctrl_B <- how(nperm = 999, blocks = meta_B$Site)

## =========================================
## 4) db-RDA (capscale) + ANOVAs
##    Model B: Treatment3 + Condition(Site)
## =========================================
cap_B <- capscale(comm_B ~ Treatment3 + Condition(Site),
                  distance = "bray", data = meta_B)

a_all   <- anova(cap_B, permutations = ctrl_B)
a_terms <- anova(cap_B, by = "terms",  permutations = ctrl_B)
a_marg  <- anova(cap_B, by = "margin", permutations = ctrl_B)
a_axis  <- anova(cap_B, by = "axis",   permutations = ctrl_B)

adjR2 <- RsquareAdj(cap_B)$adj.r.squared

write.csv(as.data.frame(a_all),   file.path(out_dir, "ITS_ModelB_capscale_ANOVA_overall.csv"))
write.csv(as.data.frame(a_terms), file.path(out_dir, "ITS_ModelB_capscale_ANOVA_terms_sequential.csv"))
write.csv(as.data.frame(a_marg),  file.path(out_dir, "ITS_ModelB_capscale_ANOVA_terms_marginal.csv"))
write.csv(as.data.frame(a_axis),  file.path(out_dir, "ITS_ModelB_capscale_ANOVA_axis.csv"))
write.csv(data.frame(AdjR2 = adjR2), file.path(out_dir, "ITS_ModelB_capscale_AdjR2.csv"), row.names = FALSE)

# subtitle: overall model (and marginal term if present)
if("Treatment3" %in% rownames(a_marg)){
  subtitle_txt <- paste0("Model ", fmt_stats(a_all[1,"F"], a_all[1,"Pr(>F)"], adjR2),
                         " | Treatment ", fmt_stats(a_marg["Treatment3","F"], a_marg["Treatment3","Pr(>F)"], adjR2))
} else {
  subtitle_txt <- paste0("Model ", fmt_stats(a_all[1,"F"], a_all[1,"Pr(>F)"], adjR2))
}

## =========================================
## 5) PERMANOVA + dispersion
## =========================================
dist_B <- vegdist(comm_B, method = "bray")

permB <- adonis2(dist_B ~ Treatment3, data = meta_B,
                 permutations = 999, strata = meta_B$Site, by = "margin")
write.csv(as.data.frame(permB), file.path(out_dir, "ITS_ModelB_PERMANOVA_adonis2_margin.csv"))

# dispersion across treatments
disp_Tr <- betadisper(dist_B, meta_B$Treatment3)

disp_tab <- tibble(
  Test = "betadisper_Treatment3",
  F = permutest(disp_Tr, permutations = ctrl_B)$tab[1,"F"],
  p = permutest(disp_Tr, permutations = ctrl_B)$tab[1,"Pr(>F)"]
)
write.csv(disp_tab, file.path(out_dir, "ITS_ModelB_Dispersion_betadisper_permutest_blocked.csv"), row.names = FALSE)

## =========================================
## 6) envfit (ALL numeric env vars) + FDR
## =========================================
drop_like <- c("SampleID","sampleid","Treatment","Treatment3","EcologicalStatus",
               "SoilPortion","Site","Block","Plot")

env_numeric <- make_env_numeric(ps_B, meta_B, drop_like)

set.seed(123)
fit_all <- envfit(cap_B, env_numeric, permutations = ctrl_B, na.rm = TRUE)

envfit_table_all <- data.frame(
  var = rownames(fit_all$vectors$arrows),
  CAP1 = fit_all$vectors$arrows[,1],
  CAP2 = fit_all$vectors$arrows[,2],
  r2  = fit_all$vectors$r,
  p   = fit_all$vectors$pvals,
  stringsAsFactors = FALSE
) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p, desc(r2))

write.csv(envfit_table_all, file.path(out_dir, "ITS_ModelB_envfit_ALL_numeric_FDR.csv"), row.names = FALSE)

df_vec <- envfit_table_all %>% filter(p_adj < 0.05)

## =========================================
## 7) Plot (95% ellipses + % axis + envfit arrows)
## =========================================
scr <- vegan::scores(cap_B, display = "sites", choices = 1:2)
df_sites <- tibble::as_tibble(scr, rownames = "SampleID") %>%
  left_join(meta_B %>% dplyr::select(SampleID, Treatment3),
            by = "SampleID")

stopifnot(!anyDuplicated(df_sites$SampleID))

eig_constrained <- cap_B$CCA$eig
cap1_pct <- 100 * eig_constrained[1] / sum(eig_constrained)
cap2_pct <- 100 * eig_constrained[2] / sum(eig_constrained)

x_lab <- paste0("CAP1 (", round(cap1_pct, 1), "%)")
y_lab <- paste0("CAP2 (", round(cap2_pct, 1), "%)")

arrow_mult <- 0.9 * max(
  diff(range(df_sites$CAP1, na.rm = TRUE)),
  diff(range(df_sites$CAP2, na.rm = TRUE))
)

col_map <- c("NATURAL"="green3", "BESE"="dodgerblue3", "BARE"="orange3")

pB <- ggplot(df_sites, aes(CAP1, CAP2, color = Treatment3)) +
  geom_point(size = 2.2, alpha = 0.85) +
  stat_ellipse(aes(group = Treatment3), type = "t", level = 0.95, linewidth = 0.8) +
  scale_color_manual(values = col_map) +
  labs(title = "ITS Model B (RZ): db-RDA + envfit",
       subtitle = subtitle_txt,
       x = x_lab, y = y_lab,
       color = "Treatment") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

if(nrow(df_vec) > 0){
  pB <- pB +
    geom_segment(
      data = df_vec,
      aes(x = 0, y = 0, xend = CAP1 * arrow_mult, yend = CAP2 * arrow_mult),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.02, "npc")),
      linewidth = 0.7,
      color = "black"
    ) +
    geom_text(
      data = df_vec,
      aes(x = CAP1 * arrow_mult * 1.07, y = CAP2 * arrow_mult * 1.07, label = var),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    )
}

ggsave(file.path(out_dir, "ITS_ModelB_dbRDA_envfit_ellipses.png"), pB, width = 6.6, height = 5.0, dpi = 600)
ggsave(file.path(out_dir, "ITS_ModelB_dbRDA_envfit_ellipses.pdf"), pB, width = 6.6, height = 5.0, device = "pdf")

cat("\nDONE. Outputs saved to: ", normalizePath(out_dir), "\n")

print(pB) 


## =========================================================
## ITS Model C beta diversity (FINAL CLEAN)
## Compare Compartment: RZ vs BS (all treatments)
## Blocked permutations by Site (Site = design factor)
## Bray (relative abundance) + db-RDA + PERMANOVA + dispersion + envfit (FDR)
## Outputs: plot + CSV tables for Supplementary
## =========================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(permute)
  library(grid)
})

set.seed(1)

## ---------- USER SETTINGS ----------
in_rds  <- "~/Library/CloudStorage/OneDrive-DeakinUniversity/Desktop/ALL DESKTOP/Phyloseq_ITS/phyloseq_ITS_fixed.rds"
out_dir <- "supp_beta_ITS_ModelC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- Helpers ----------
relab_comm <- function(ps){
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  mat <- as(otu_table(ps_rel), "matrix")
  if (taxa_are_rows(ps_rel)) mat <- t(mat)   # samples x taxa
  mat
}

to_num_if_possible <- function(x){
  if(is.numeric(x) || is.integer(x)) return(as.numeric(x))
  if(is.factor(x)) x <- as.character(x)
  if(is.character(x)){
    x2 <- suppressWarnings(as.numeric(x))
    if(mean(!is.na(x2)) >= 0.7) return(x2)
  }
  x
}

fmt_stats <- function(Fval, pval, adjR2){
  p_txt <- ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))
  paste0("F = ", sprintf("%.2f", Fval),
         ", p = ", p_txt,
         ", Adj R² = ", sprintf("%.3f", adjR2))
}

make_env_numeric <- function(ps, meta, drop_like){
  md <- as(sample_data(ps), "data.frame")
  if(!"SampleID" %in% names(md)) md$SampleID <- rownames(md)
  
  md2 <- md[, !tolower(names(md)) %in% tolower(drop_like), drop = FALSE]
  md3 <- as.data.frame(lapply(md2, to_num_if_possible))
  md3$SampleID <- md$SampleID
  
  env_numeric <- md3 %>%
    dplyr::select(where(is.numeric)) %>%
    mutate(SampleID = md$SampleID)
  
  env_numeric <- env_numeric %>%
    dplyr::semi_join(meta %>% dplyr::select(SampleID), by = "SampleID") %>%
    dplyr::right_join(meta %>% dplyr::select(SampleID), by = "SampleID") %>%
    tibble::column_to_rownames("SampleID")
  
  keep <- colSums(!is.na(env_numeric)) > 0
  env_numeric[, keep, drop = FALSE]
}

## ---------- Load ----------
ps <- readRDS(in_rds)

## =========================================
## 1) SUBSET: BS + RZ only + Site (all treatments)
## =========================================
ps_C <- subset_samples(ps,
                       !is.na(Site) & !is.na(SoilPortion) &
                         SoilPortion %in% c("BS","RZ","BULK","BULK_SOIL","BULKSOIL","BULK SOIL","RHIZOSPHERE"))
ps_C <- prune_samples(sample_sums(ps_C) > 0, ps_C)
ps_C <- prune_taxa(taxa_sums(ps_C) > 0, ps_C)

meta_C <- as(sample_data(ps_C), "data.frame")
if(!"SampleID" %in% colnames(meta_C)) meta_C$SampleID <- rownames(meta_C)

# Recode Compartment (BS vs RZ)
sp_raw <- toupper(trimws(as.character(meta_C$SoilPortion)))
meta_C$Compartment <- dplyr::case_when(
  sp_raw %in% c("RZ","RHIZOSPHERE") ~ "RZ",
  sp_raw %in% c("BS","BULK","BULK_SOIL","BULKSOIL","BULK SOIL") ~ "BS",
  TRUE ~ NA_character_
)
meta_C$Compartment <- factor(meta_C$Compartment, levels = c("BS","RZ"))
meta_C$Site <- factor(meta_C$Site)

# Keep only BS/RZ after recode
keep_ids <- meta_C$SampleID[!is.na(meta_C$Compartment)]
ps_C <- prune_samples(keep_ids, ps_C)
ps_C <- prune_taxa(taxa_sums(ps_C) > 0, ps_C)

# Refresh meta after pruning
meta_C <- as(sample_data(ps_C), "data.frame")
if(!"SampleID" %in% colnames(meta_C)) meta_C$SampleID <- rownames(meta_C)

sp_raw <- toupper(trimws(as.character(meta_C$SoilPortion)))
meta_C$Compartment <- dplyr::case_when(
  sp_raw %in% c("RZ","RHIZOSPHERE") ~ "RZ",
  sp_raw %in% c("BS","BULK","BULK_SOIL","BULKSOIL","BULK SOIL") ~ "BS",
  TRUE ~ NA_character_
)
meta_C$Compartment <- factor(meta_C$Compartment, levels = c("BS","RZ"))
meta_C$Site <- factor(meta_C$Site)

stopifnot(length(unique(meta_C$Site)) > 1)

cat("\nSample counts by Compartment x Site:\n")
print(table(meta_C$Compartment, meta_C$Site))

## =========================================
## 2) COMMUNITY MATRIX + align meta
## =========================================
comm_C <- relab_comm(ps_C)
stopifnot(identical(rownames(comm_C), sample_names(ps_C)))

meta_C <- meta_C %>%
  dplyr::mutate(SampleID = as.character(SampleID)) %>%
  dplyr::slice(match(rownames(comm_C), SampleID))

stopifnot(identical(meta_C$SampleID, rownames(comm_C)))

## =========================================
## 3) Blocked permutations by Site
## =========================================
ctrl_C <- how(nperm = 999, blocks = meta_C$Site)

## =========================================
## 4) db-RDA (capscale) + ANOVAs (blocked)
##    Model C: Compartment + Condition(Site)
## =========================================
cap_C <- capscale(comm_C ~ Compartment + Condition(Site),
                  distance = "bray", data = meta_C)

a_all   <- anova(cap_C, permutations = ctrl_C)
a_terms <- anova(cap_C, by = "terms",  permutations = ctrl_C)
a_marg  <- anova(cap_C, by = "margin", permutations = ctrl_C)
a_axis  <- anova(cap_C, by = "axis",   permutations = ctrl_C)

adjR2 <- RsquareAdj(cap_C)$adj.r.squared

write.csv(as.data.frame(a_all),   file.path(out_dir, "ITS_ModelC_capscale_ANOVA_overall.csv"))
write.csv(as.data.frame(a_terms), file.path(out_dir, "ITS_ModelC_capscale_ANOVA_terms_sequential.csv"))
write.csv(as.data.frame(a_marg),  file.path(out_dir, "ITS_ModelC_capscale_ANOVA_terms_marginal.csv"))
write.csv(as.data.frame(a_axis),  file.path(out_dir, "ITS_ModelC_capscale_ANOVA_axis.csv"))
write.csv(data.frame(AdjR2 = adjR2), file.path(out_dir, "ITS_ModelC_capscale_AdjR2.csv"), row.names = FALSE)

# Subtitle (use marginal Compartment test if available)
if("Compartment" %in% rownames(a_marg)){
  subtitle_txt <- paste0("Compartment ", fmt_stats(a_marg["Compartment","F"], a_marg["Compartment","Pr(>F)"], adjR2))
} else {
  subtitle_txt <- paste0("Model ", fmt_stats(a_all[1,"F"], a_all[1,"Pr(>F)"], adjR2))
}

## =========================================
## 5) PERMANOVA + dispersion (Site-blocked)
## =========================================
dist_C <- vegdist(comm_C, method = "bray")

permC <- adonis2(dist_C ~ Compartment, data = meta_C,
                 permutations = 999, strata = meta_C$Site, by = "margin")
write.csv(as.data.frame(permC), file.path(out_dir, "ITS_ModelC_PERMANOVA_adonis2_margin.csv"))

disp_C <- betadisper(dist_C, meta_C$Compartment)
disp_tab <- tibble(
  Test = "betadisper_Compartment",
  F = permutest(disp_C, permutations = ctrl_C)$tab[1,"F"],
  p = permutest(disp_C, permutations = ctrl_C)$tab[1,"Pr(>F)"]
)
write.csv(disp_tab, file.path(out_dir, "ITS_ModelC_Dispersion_betadisper_permutest_blocked.csv"), row.names = FALSE)

## =========================================
## 6) envfit (ALL numeric env vars) + FDR (blocked)
## =========================================
drop_like <- c("SampleID","sampleid","Treatment","EcologicalStatus",
               "SoilPortion","Compartment","Site","Block","Plot")

env_numeric <- make_env_numeric(ps_C, meta_C, drop_like)

set.seed(123)
fit_all <- envfit(cap_C, env_numeric, permutations = ctrl_C, na.rm = TRUE)

# IMPORTANT: do NOT assume CAP2 exists (rank-1 models often give CAP1 + MDS1)
vec_mat <- as.data.frame(scores(fit_all, display = "vectors"))
vec_mat$var <- rownames(vec_mat)

# take first two axes returned (whatever their names are)
ax_names <- colnames(vec_mat)[1:2]
colnames(vec_mat)[1:2] <- c("AX1","AX2")

envfit_table_all <- data.frame(
  var = vec_mat$var,
  AX1 = vec_mat$AX1,
  AX2 = vec_mat$AX2,
  r2  = fit_all$vectors$r,
  p   = fit_all$vectors$pvals,
  stringsAsFactors = FALSE
) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(p, desc(r2))

write.csv(envfit_table_all, file.path(out_dir, "ITS_ModelC_envfit_ALL_numeric_FDR.csv"), row.names = FALSE)

df_vec <- envfit_table_all %>% filter(p_adj < 0.05)


## =========================================
## 7) Plot (95% ellipses + % axis + envfit arrows)
## Robust to CAP1+MDS1 (rank-1 constrained models)
## =========================================

scr <- vegan::scores(cap_C, display = "sites", choices = 1:2)

# Make a clean sites table with guaranteed numeric axes
df_sites <- as.data.frame(scr) %>%
  tibble::rownames_to_column("SampleID")

# Identify axis columns from scores() (usually CAP1 and MDS1)
axis_cols <- setdiff(colnames(df_sites), "SampleID")
stopifnot(length(axis_cols) >= 2)

# Standardize names + force numeric
df_sites <- df_sites %>%
  mutate(
    AX1 = as.numeric(.data[[axis_cols[1]]]),
    AX2 = as.numeric(.data[[axis_cols[2]]])
  ) %>%
  select(SampleID, AX1, AX2)

# Join grouping variable AFTER axes are locked in
df_sites <- df_sites %>%
  left_join(meta_C %>% dplyr::select(SampleID, Compartment), by = "SampleID")

stopifnot(!anyDuplicated(df_sites$SampleID))
stopifnot(is.numeric(df_sites$AX1), is.numeric(df_sites$AX2))

# Axis labels (only CAP1 % is interpretable as constrained)
eig_constrained <- cap_C$CCA$eig
cap1_pct <- 100 * eig_constrained[1] / sum(eig_constrained)

x_lab <- paste0(axis_cols[1], " (", round(cap1_pct, 1), "% constrained)")
y_lab <- paste0(axis_cols[2], " (unconstrained)")

arrow_mult <- 0.9 * max(
  diff(range(df_sites$AX1, na.rm = TRUE)),
  diff(range(df_sites$AX2, na.rm = TRUE))
)

col_map <- c("BS"="grey30", "RZ"="firebrick3")

pC <- ggplot(df_sites, aes(AX1, AX2, color = Compartment)) +
  geom_point(size = 2.2, alpha = 0.85) +
  stat_ellipse(aes(group = Compartment), type = "t", level = 0.95, linewidth = 0.8) +
  scale_color_manual(values = col_map) +
  labs(title = "ITS Model C: RZ vs BS (db-RDA + envfit)",
       subtitle = subtitle_txt,
       x = x_lab, y = y_lab,
       color = "Compartment") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# df_vec must use AX1/AX2 now (make sure your envfit table uses these names)
if(nrow(df_vec) > 0){
  pC <- pC +
    geom_segment(
      data = df_vec,
      aes(x = 0, y = 0, xend = AX1 * arrow_mult, yend = AX2 * arrow_mult),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.02, "npc")),
      linewidth = 0.7,
      color = "black"
    ) +
    geom_text(
      data = df_vec,
      aes(x = AX1 * arrow_mult * 1.07, y = AX2 * arrow_mult * 1.07, label = var),
      inherit.aes = FALSE,
      size = 3,
      color = "black"
    )
}

print(pC)

ggsave(file.path(out_dir, "ITS_ModelC_dbRDA_envfit_ellipses.png"), pC, width = 6.6, height = 5.0, dpi = 600)
ggsave(file.path(out_dir, "ITS_ModelC_dbRDA_envfit_ellipses.pdf"), pC, width = 6.6, height = 5.0, device = "pdf")
