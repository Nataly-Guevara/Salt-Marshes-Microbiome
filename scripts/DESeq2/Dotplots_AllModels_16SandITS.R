# ============================================================
# DESeq2 Dotplots — Models A/B/C — 16S and ITS
#
# Purpose:
#   Generate dotplots of Top 10 differentially abundant Orders
#   (ranked by adjusted p-value) for DESeq2 outputs exported to CSV.
#
# Models:
#   Model A design:
#     EcologicalStatus * SoilPortion + Site
#
#   Model B design:
#     Treatment effects (RZ only) — contrasts in all_contrasts export
#
#   Model C design:
#     SoilPortion (RZ vs BS) — contrasts in SoilPortion export
#
# Input (relative to project dir):
#   inputs/
#     DESeq2_16S_ModelA_all_contrasts.csv
#     DESeq2_ITS_ModelA_all_contrasts.csv
#     DESeq2_16S_ModelB_all_contrasts.csv
#     DESeq2_ITS_ModelB_all_contrasts.csv
#     DESeq2_16S_ModelC_SoilPortion.csv
#     DESeq2_ITS_ModelC_SoilPortion.csv
#
# Output (auto-created):
#   outputs/figures/deseq2_dotplots_ModelA/ (and B/C)
#   outputs/tables/deseq2_ModelA/           (and B/C)
# ============================================================


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/DEAKIN/Supplementary/DESeq2_main/Dotplots")


install.packages("purrr")
library(purrr)
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(writexl)
})

# ---------------------------
# 0) Project structure
# ---------------------------
project_dir <- getwd()

input_dir  <- file.path(project_dir, "inputs")
output_dir <- file.path(project_dir, "outputs")

stopifnot(dir.exists(input_dir))

# ---------------------------
# 1) Helpers
# ---------------------------
sig_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

std_contrast <- function(x) {
  x %>%
    as.character() %>%
    str_squish() %>%
    str_replace_all("VS\\.", "vs") %>%
    str_replace_all("\\bVS\\b", "vs") %>%
    str_replace_all("\\s+vs\\s+", " vs ")
}

# Extract g1/g2 from "A vs B ..." (ignores trailing text like "within BS")
get_groups <- function(contrast) {
  main <- str_split(contrast, " within | \\(| — | - ", simplify = TRUE)[,1]
  parts <- str_split(main, " vs ", simplify = TRUE)
  tibble(
    g1 = str_squish(parts[,1]),
    g2 = str_squish(parts[,2])
  )
}

# Robust reader: harmonize log2FC column name and required fields
read_clean <- function(file, target, model) {
  df <- read_csv(file, show_col_types = FALSE)
  
  # Harmonize log2FC column
  if (!("log2FC" %in% names(df)) && ("log2FoldChange" %in% names(df))) {
    df <- df %>% mutate(log2FC = log2FoldChange)
  }
  # Some exports might already have log2FC and not log2FoldChange
  stopifnot(all(c("baseMean", "padj", "log2FC", "Order") %in% names(df)))
  
  # Contrast column:
  # - Models A/B usually have Contrast
  # - Model C sometimes may omit; if missing, create a default label
  if (!("Contrast" %in% names(df))) {
    df <- df %>% mutate(Contrast = "RZ vs BS")
  }
  
  df %>%
    mutate(
      Target = target,
      Model = model,
      Contrast = std_contrast(Contrast),
      Order = ifelse(is.na(Order) | Order == "" | Order == "NA", "Unclassified", as.character(Order))
    )
}

# Order-proxy summary (aggregated across features within each Order)
make_order_proxy <- function(df) {
  df %>%
    filter(!is.na(padj)) %>%
    group_by(Target, Model, Contrast, Order) %>%
    summarise(
      n_features   = n(),
      baseMean_sum = sum(baseMean, na.rm = TRUE),
      log2FC_w     = weighted.mean(log2FC, w = baseMean, na.rm = TRUE),
      padj_min     = min(padj, na.rm = TRUE),
      pvalue_min   = if ("pvalue" %in% names(df)) min(pvalue, na.rm = TRUE) else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      stars = sig_stars(padj_min),
      label = sprintf("%.2f%s", log2FC_w, stars)
    )
}

# Dotplot: Top 10 Orders by padj_min
plot_dotplot <- function(order_proxy, target, model, contrast, fig_dir) {
  
  dfp <- order_proxy %>%
    filter(Target == target, Model == model, Contrast == contrast) %>%
    arrange(padj_min) %>%
    slice_head(n = 10)
  
  if (nrow(dfp) == 0) {
    warning("No data for: ", target, " | ", model, " | ", contrast)
    return(NULL)
  }
  
  groups <- get_groups(contrast)
  g1 <- groups$g1[1]
  g2 <- groups$g2[1]
  
  dfp <- dfp %>%
    mutate(
      EnrichedIn = ifelse(log2FC_w > 0, g1, g2),
      Order = fct_reorder(Order, log2FC_w)
    )
  
  p <- ggplot(dfp, aes(x = log2FC_w, y = Order)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point(aes(size = baseMean_sum, color = EnrichedIn), alpha = 0.85) +
    geom_text(aes(label = label), nudge_x = 0.10, size = 3.5, check_overlap = TRUE) +
    scale_size_continuous(name = "baseMean (sum)") +
    labs(
      title = paste0(target, " — ", model, " — ", contrast),
      subtitle = "Top 10 Orders by adjusted p-value",
      x = "log2 Fold Change",
      y = "Order",
      color = "Enriched in"
    ) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  file_name <- paste0(
    "Dotplot_", target, "_", model, "_",
    str_replace_all(contrast, "[^A-Za-z0-9]+", "_"),
    ".png"
  )
  
  ggsave(file.path(fig_dir, file_name), p, width = 9, height = 6, dpi = 300)
  p
}

# ---------------------------
# 2) File map (your exact filenames)
# ---------------------------
files <- tibble::tribble(
  ~Model,  ~Target, ~File,
  "ModelA","16S",   "DESeq2_16S_ModelA_all_contrasts.csv",
  "ModelA","ITS",   "DESeq2_ITS_ModelA_all_contrasts.csv",
  "ModelB","16S",   "DESeq2_16S_ModelB_all_contrasts.csv",
  "ModelB","ITS",   "DESeq2_ITS_ModelB_all_contrasts.csv",
  "ModelC","16S",   "DESeq2_16S_ModelC_SoilPortion.csv",
  "ModelC","ITS",   "DESeq2_ITS_ModelC_SoilPortion.csv"
) %>%
  mutate(Path = file.path(input_dir, File))

# Confirm all inputs exist
missing <- files %>% filter(!file.exists(Path))
if (nrow(missing) > 0) {
  stop("Missing input files in inputs/: \n", paste(missing$File, collapse = "\n"))
}

# ---------------------------
# 3) Read all results
# ---------------------------
res_all <- purrr::pmap_dfr(
  list(files$Path, files$Target, files$Model),
  ~ read_clean(..1, target = ..2, model = ..3)
)

# ---------------------------
# 4) Build Order-proxy summary
# ---------------------------
order_proxy <- make_order_proxy(res_all)

# ---------------------------
# 5) Loop through models: export tables + plots
# ---------------------------
models <- c("ModelA", "ModelB", "ModelC")

for (m in models) {
  
  # output folders per model
  fig_dir <- file.path(output_dir, "figures", paste0("deseq2_dotplots_", m))
  tab_dir <- file.path(output_dir, "tables",  paste0("deseq2_", m))
  
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
  
  # subset for this model
  op_m <- order_proxy %>% filter(Model == m)
  
  # export tables
  write_csv(op_m, file.path(tab_dir, paste0("OrderProxy_", m, "_ALL.csv")))
  write_csv(op_m %>% filter(padj_min < 0.05),
            file.path(tab_dir, paste0("OrderProxy_", m, "_SIG_padj<0.05.csv")))
  
  # Excel supplementary per model
  writexl::write_xlsx(
    list(
      "OrderProxy_ALL" = op_m,
      "OrderProxy_SIG_padj<0.05" = op_m %>% filter(padj_min < 0.05)
    ),
    path = file.path(tab_dir, paste0("Supplement_Table_DESeq2_", m, "_16S_ITS.xlsx"))
  )
  
  # plots per target and per contrast found
  for (tgt in c("16S", "ITS")) {
    
    contrasts <- op_m %>%
      filter(Target == tgt) %>%
      pull(Contrast) %>%
      unique()
    
    for (ct in contrasts) {
      p <- plot_dotplot(order_proxy, target = tgt, model = m, contrast = ct, fig_dir = fig_dir)
      if (!is.null(p)) print(p)
    }
  }
  
  message("DONE: ", m, " | Figures -> ", fig_dir, " | Tables -> ", tab_dir)
}

message("\nALL MODELS COMPLETE.\n")

