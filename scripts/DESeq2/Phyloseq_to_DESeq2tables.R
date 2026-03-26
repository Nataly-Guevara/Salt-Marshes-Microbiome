# ============================================================
# 0) Working directory + packages
# ============================================================

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/DEAKIN/Supplementary/DESeq2_main/DESeq2")

suppressPackageStartupMessages({
  library(phyloseq)
  library(DESeq2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
})

# VERY IMPORTANT: ensure treatment contrasts
options(contrasts = c("contr.treatment", "contr.poly"))

# ============================================================
# 1) Paths
# ============================================================

project_dir <- getwd()

in16 <- file.path(project_dir, "inputs", "phyloseq", "phyloseq_16S_clean.rds")
inIT <- file.path(project_dir, "inputs", "phyloseq", "phyloseq_ITS_final.rds")

outA <- file.path(project_dir, "outputs", "deseq2_tables", "ModelA")
outB <- file.path(project_dir, "outputs", "deseq2_tables", "ModelB")
outC <- file.path(project_dir, "outputs", "deseq2_tables", "ModelC")

dir.create(outA, recursive = TRUE, showWarnings = FALSE)
dir.create(outB, recursive = TRUE, showWarnings = FALSE)
dir.create(outC, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(in16), file.exists(inIT))

# ============================================================
# 2) Helper functions
# ============================================================

clean_tax <- function(ps) {
  taxa_names(ps) <- make.names(taxa_names(ps), unique = TRUE)
  ps
}

tax_df_from_ps <- function(ps_ord) {
  as.data.frame(tax_table(ps_ord)) %>%
    rownames_to_column("ASV")
}

annotate_res <- function(res, tax_df, contrast_label, g1, g2) {
  as.data.frame(res) %>%
    rownames_to_column("ASV") %>%
    left_join(tax_df, by = "ASV") %>%
    mutate(
      Contrast = contrast_label,
      Group1 = g1,
      Group2 = g2,
      EnrichedIn = case_when(
        is.na(log2FoldChange) ~ NA_character_,
        log2FoldChange > 0 ~ g1,
        log2FoldChange < 0 ~ g2,
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(padj)
}

# ============================================================
# 3) Load phyloseq objects
# ============================================================

ps16 <- readRDS(in16)
psit <- readRDS(inIT)

# Standardize common factors
sample_data(ps16)$SoilPortion <- factor(sample_data(ps16)$SoilPortion, levels = c("BS","RZ"))
sample_data(psit)$SoilPortion <- factor(sample_data(psit)$SoilPortion, levels = c("BS","RZ"))

sample_data(ps16)$Site <- factor(sample_data(ps16)$Site)
sample_data(psit)$Site <- factor(sample_data(psit)$Site)

sample_data(ps16)$EcologicalStatus <- factor(sample_data(ps16)$EcologicalStatus,
                                             levels = c("DEGREST","NATURAL"))
sample_data(psit)$EcologicalStatus <- factor(sample_data(psit)$EcologicalStatus,
                                             levels = c("DEGREST","NATURAL"))

# ============================================================
# ======================= MODEL A =============================
# ============================================================

ps_ord <- tax_glom(ps16, taxrank = "Order")
ps_ord <- clean_tax(ps_ord)
ps_ord <- filter_taxa(ps_ord, function(x) sum(x > 5) > (0.2 * length(x)), prune = TRUE)

sample_data(ps_ord)$SoilPortion <- factor(sample_data(ps_ord)$SoilPortion, levels = c("BS","RZ"))

tax_df <- tax_df_from_ps(ps_ord)

dds <- phyloseq_to_deseq2(ps_ord, ~ EcologicalStatus * SoilPortion + Site)
dds <- DESeq(dds, fitType = "parametric")

res_RZvsBS_DEGREST <- results(dds, name = "SoilPortion_RZ_vs_BS")
df1 <- annotate_res(res_RZvsBS_DEGREST, tax_df,
                    "RZ vs BS within DEGREST","RZ","BS")

res_RZvsBS_NATURAL <- results(dds,
                              contrast = list(c("SoilPortion_RZ_vs_BS",
                                                "EcologicalStatusNATURAL.SoilPortionRZ")))
df2 <- annotate_res(res_RZvsBS_NATURAL, tax_df,
                    "RZ vs BS within NATURAL","RZ","BS")

res_NATvsDEG_BS <- results(dds,
                           name = "EcologicalStatus_NATURAL_vs_DEGREST")
df3 <- annotate_res(res_NATvsDEG_BS, tax_df,
                    "NATURAL vs DEGREST within BS",
                    "NATURAL","DEGREST")

res_NATvsDEG_RZ <- results(dds,
                           contrast = list(c("EcologicalStatus_NATURAL_vs_DEGREST",
                                             "EcologicalStatusNATURAL.SoilPortionRZ")))
df4 <- annotate_res(res_NATvsDEG_RZ, tax_df,
                    "NATURAL vs DEGREST within RZ",
                    "NATURAL","DEGREST")

all_contrasts_A <- bind_rows(df1,df2,df3,df4)

write_csv(all_contrasts_A,
          file.path(outA,"DESeq2_16S_ModelA_all_contrasts.csv"))

# ============================================================
# ======================= MODEL A (ITS) ======================
# ============================================================

ps_ord_it <- tax_glom(psit, taxrank = "Order")
ps_ord_it <- clean_tax(ps_ord_it)
ps_ord_it <- filter_taxa(ps_ord_it,
                         function(x) sum(x > 5) > (0.2 * length(x)),
                         prune = TRUE)

sample_data(ps_ord_it)$SoilPortion <-
  factor(sample_data(ps_ord_it)$SoilPortion,
         levels = c("BS","RZ"))

sample_data(ps_ord_it)$EcologicalStatus <-
  factor(sample_data(ps_ord_it)$EcologicalStatus,
         levels = c("DEGREST","NATURAL"))

sample_data(ps_ord_it)$Site <- factor(sample_data(ps_ord_it)$Site)

tax_df_it <- tax_df_from_ps(ps_ord_it)

dds_it <- phyloseq_to_deseq2(ps_ord_it,
                             ~ EcologicalStatus * SoilPortion + Site)

dds_it <- DESeq(dds_it, fitType = "parametric")

# 1) RZ vs BS within DEGREST
res1 <- results(dds_it, name = "SoilPortion_RZ_vs_BS")
df1 <- annotate_res(res1, tax_df_it,
                    "RZ vs BS within DEGREST","RZ","BS")

# 2) RZ vs BS within NATURAL
res2 <- results(dds_it,
                contrast = list(c("SoilPortion_RZ_vs_BS",
                                  "EcologicalStatusNATURAL.SoilPortionRZ")))
df2 <- annotate_res(res2, tax_df_it,
                    "RZ vs BS within NATURAL","RZ","BS")

# 3) NATURAL vs DEGREST within BS
res3 <- results(dds_it,
                name = "EcologicalStatus_NATURAL_vs_DEGREST")
df3 <- annotate_res(res3, tax_df_it,
                    "NATURAL vs DEGREST within BS",
                    "NATURAL","DEGREST")

# 4) NATURAL vs DEGREST within RZ
res4 <- results(dds_it,
                contrast = list(c("EcologicalStatus_NATURAL_vs_DEGREST",
                                  "EcologicalStatusNATURAL.SoilPortionRZ")))
df4 <- annotate_res(res4, tax_df_it,
                    "NATURAL vs DEGREST within RZ",
                    "NATURAL","DEGREST")

all_contrasts_A_ITS <- bind_rows(df1,df2,df3,df4)

write_csv(all_contrasts_A_ITS,
          file.path(outA,
                    "DESeq2_ITS_ModelA_all_contrasts.csv"))

# ============================================================
# ======================= MODEL B =============================
# ============================================================

run_modelB <- function(ps, label, out_dir) {
  
  # Special ITS fix (NATURAL stored as "NA")
  if (label == "ITS") {
    md <- as.data.frame(sample_data(ps))
    md$Treatment <- as.character(md$Treatment)
    md$Treatment[md$Treatment == "NA"] <- "NATURAL"
    md$Treatment <- factor(md$Treatment,
                           levels = c("NATURAL","BESE","BARE"))
    sample_data(ps)$Treatment <- md$Treatment
    ps <- prune_samples(!is.na(sample_data(ps)$Treatment), ps)
    ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  }
  
  ps_rz <- subset_samples(ps, SoilPortion == "RZ")
  ps_rz <- prune_samples(!is.na(sample_data(ps_rz)$Treatment), ps_rz)
  
  ps_ord <- tax_glom(ps_rz, taxrank = "Order")
  ps_ord <- clean_tax(ps_ord)
  ps_ord <- filter_taxa(ps_ord,
                        function(x) sum(x > 5) > (0.2 * length(x)),
                        prune = TRUE)
  
  sample_data(ps_ord)$Treatment <-
    factor(sample_data(ps_ord)$Treatment,
           levels = c("NATURAL","BESE","BARE"))
  sample_data(ps_ord)$Treatment <-
    relevel(sample_data(ps_ord)$Treatment, ref = "NATURAL")
  
  sample_data(ps_ord)$Site <- factor(sample_data(ps_ord)$Site)
  
  tax_df <- tax_df_from_ps(ps_ord)
  
  dds <- phyloseq_to_deseq2(ps_ord, ~ Site + Treatment)
  dds <- DESeq(dds, fitType = "parametric")
  
  res1 <- results(dds,
                  contrast = c("Treatment","BESE","NATURAL"))
  df1 <- annotate_res(res1, tax_df,
                      "BESE vs NATURAL",
                      "BESE","NATURAL")
  
  res2 <- results(dds,
                  contrast = c("Treatment","BARE","NATURAL"))
  df2 <- annotate_res(res2, tax_df,
                      "BARE vs NATURAL",
                      "BARE","NATURAL")
  
  res3 <- results(dds,
                  contrast = c("Treatment","BESE","BARE"))
  df3 <- annotate_res(res3, tax_df,
                      "BESE vs BARE",
                      "BESE","BARE")
  
  all_contrasts <- bind_rows(df1,df2,df3)
  
  write_csv(all_contrasts,
            file.path(out_dir,
                      paste0("DESeq2_",label,
                             "_ModelB_all_contrasts.csv")))
}

run_modelB(ps16,"16S",outB)
run_modelB(psit,"ITS",outB)

# ============================================================
# ======================= MODEL C =============================
# ============================================================

run_modelC <- function(ps, label, out_dir) {
  
  ps_ord <- tax_glom(ps, taxrank="Order")
  ps_ord <- clean_tax(ps_ord)
  
  sample_data(ps_ord)$SoilPortion <-
    factor(sample_data(ps_ord)$SoilPortion,
           levels = c("BS","RZ"))
  
  sample_data(ps_ord)$Site <- factor(sample_data(ps_ord)$Site)
  
  dds <- phyloseq_to_deseq2(ps_ord, ~ SoilPortion + Site)
  dds <- DESeq(dds)
  
  tax_df <- tax_df_from_ps(ps_ord)
  
  res <- results(dds,
                 contrast = c("SoilPortion","RZ","BS"))
  
  res_df <- annotate_res(res, tax_df,
                         "RZ vs BS",
                         "RZ","BS")
  
  write_csv(res_df,
            file.path(out_dir,
                      paste0("DESeq2_",label,
                             "_ModelC_SoilPortion.csv")))
}

run_modelC(ps16,"16S",outC)
run_modelC(psit,"ITS",outC)

message("All models completed successfully.")
