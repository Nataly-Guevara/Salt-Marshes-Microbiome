# ============================================================
# MICROBIAL OVERVIEW — PHYLUM-LEVEL EXPLORATION
# Objects used:
#   ps_prok = final prokaryotic phyloseq object
#   ps_fun  = final fungal phyloseq object
#
# Aim:
#   Explore dominant phyla across:
#     1. Overall
#     2. Bulk soil
#     3. Rhizosphere
#
# No treatment comparison yet.
# No graph yet.
# ============================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)

# ============================================================
# 1. Output folder
# ============================================================

dir.create(
  "microbial_overview_outputs",
  showWarnings = FALSE
)

# ============================================================
# 2. Check objects exist
# ============================================================

if (!exists("ps_prok")) {
  stop("Object ps_prok was not found in the R session.")
}

if (!exists("ps_fun")) {
  stop("Object ps_fun was not found in the R session.")
}

# ============================================================
# 3. Helper functions
# ============================================================

clean_taxon_label <- function(x, rank_label = "Phylum") {
  
  x <- as.character(x)
  
  x <- case_when(
    is.na(x) ~ paste0("Unclassified_", rank_label),
    x == "" ~ paste0("Unclassified_", rank_label),
    x %in% c("NA", "Unclassified", "unclassified", "unidentified", "Unidentified") ~ paste0("Unclassified_", rank_label),
    TRUE ~ x
  )
  
  x %>%
    str_replace("^p__", "") %>%
    str_replace("^k__", "") %>%
    str_replace("^d__", "") %>%
    str_replace("^Phylum_", "") %>%
    str_replace("^Unclassified_", "Unclassified ")
}


standardise_soil_habitat <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  case_when(
    x %in% c("BS", "BULK", "BULK_SOIL") ~ "Bulk soil",
    x %in% c("RZ", "RHIZOSPHERE") ~ "Rhizosphere",
    TRUE ~ x
  )
}


make_phylum_summary <- function(ps_object,
                                community_label,
                                top_n = 10) {
  
  # -----------------------------
  # Basic object audit
  # -----------------------------
  
  cat("\n========================================\n")
  cat("OBJECT AUDIT:", community_label, "\n")
  cat("========================================\n")
  
  cat("Samples:", phyloseq::nsamples(ps_object), "\n")
  cat("ASVs:", phyloseq::ntaxa(ps_object), "\n")
  cat("Reads:", sum(phyloseq::sample_sums(ps_object)), "\n")
  
  # -----------------------------
  # Metadata
  # -----------------------------
  
  meta <- data.frame(
    phyloseq::sample_data(ps_object)
  )
  
  meta$SampleID <- rownames(meta)
  
  if (!"SoilPortion" %in% names(meta)) {
    stop(
      "SoilPortion column not found. Available sample_data columns are:\n",
      paste(names(meta), collapse = ", ")
    )
  }
  
  meta_clean <- meta %>%
    transmute(
      SampleID = SampleID,
      Soil_habitat = standardise_soil_habitat(SoilPortion)
    )
  
  cat("\nSoil habitat counts:\n")
  print(
    meta_clean %>%
      count(Soil_habitat)
  )
  
  # -----------------------------
  # Taxonomy
  # -----------------------------
  
  tax <- data.frame(
    phyloseq::tax_table(ps_object)
  )
  
  tax$ASV <- rownames(tax)
  
  if (!"Phylum" %in% names(tax)) {
    stop(
      "Phylum column not found. Available taxonomy columns are:\n",
      paste(names(tax), collapse = ", ")
    )
  }
  
  tax_clean <- tax %>%
    transmute(
      ASV = ASV,
      Phylum = clean_taxon_label(Phylum, rank_label = "Phylum")
    )
  
  # -----------------------------
  # Convert to sample-level relative abundance
  # -----------------------------
  
  ps_rel <- phyloseq::transform_sample_counts(
    ps_object,
    function(x) x / sum(x)
  )
  
  otu <- as(
    phyloseq::otu_table(ps_rel),
    "matrix"
  )
  
  if (phyloseq::taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  
  otu_long <- as.data.frame(otu) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(
      cols = -SampleID,
      names_to = "ASV",
      values_to = "Relative_abundance"
    )
  
  # -----------------------------
  # Collapse ASVs to phylum per sample
  # -----------------------------
  
  sample_phylum <- otu_long %>%
    left_join(
      tax_clean,
      by = "ASV"
    ) %>%
    left_join(
      meta_clean,
      by = "SampleID"
    ) %>%
    group_by(
      SampleID,
      Soil_habitat,
      Phylum
    ) %>%
    summarise(
      Relative_abundance = sum(Relative_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # Add Overall grouping
  # -----------------------------
  
  sample_phylum_with_overall <- bind_rows(
    sample_phylum,
    sample_phylum %>%
      mutate(
        Soil_habitat = "Overall"
      )
  )
  
  # -----------------------------
  # Summarise mean RA and prevalence
  # -----------------------------
  
  phylum_summary_all <- sample_phylum_with_overall %>%
    group_by(
      Community = community_label,
      Soil_habitat,
      Phylum
    ) %>%
    summarise(
      mean_relative_abundance = mean(Relative_abundance, na.rm = TRUE),
      median_relative_abundance = median(Relative_abundance, na.rm = TRUE),
      prevalence = mean(Relative_abundance > 0, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    mutate(
      mean_relative_abundance_percent = mean_relative_abundance * 100,
      median_relative_abundance_percent = median_relative_abundance * 100,
      prevalence_percent = prevalence * 100
    )
  
  # -----------------------------
  # Choose top phyla by Overall mean abundance
  # -----------------------------
  
  top_phyla <- phylum_summary_all %>%
    filter(
      Soil_habitat == "Overall"
    ) %>%
    arrange(
      desc(mean_relative_abundance_percent)
    ) %>%
    slice_head(
      n = top_n
    ) %>%
    pull(Phylum)
  
  # -----------------------------
  # Collapse non-top phyla into Other phyla
  # -----------------------------
  
  sample_phylum_top <- sample_phylum_with_overall %>%
    mutate(
      Phylum_plot = if_else(
        Phylum %in% top_phyla,
        Phylum,
        "Other phyla"
      )
    ) %>%
    group_by(
      Community = community_label,
      SampleID,
      Soil_habitat,
      Phylum_plot
    ) %>%
    summarise(
      Relative_abundance = sum(Relative_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  phylum_summary_top <- sample_phylum_top %>%
    group_by(
      Community,
      Soil_habitat,
      Phylum = Phylum_plot
    ) %>%
    summarise(
      mean_relative_abundance = mean(Relative_abundance, na.rm = TRUE),
      median_relative_abundance = median(Relative_abundance, na.rm = TRUE),
      prevalence = mean(Relative_abundance > 0, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    ) %>%
    mutate(
      mean_relative_abundance_percent = mean_relative_abundance * 100,
      median_relative_abundance_percent = median_relative_abundance * 100,
      prevalence_percent = prevalence * 100
    )
  
  # -----------------------------
  # Reporting table: Overall only
  # -----------------------------
  
  reporting_table <- phylum_summary_top %>%
    filter(
      Soil_habitat == "Overall"
    ) %>%
    arrange(
      desc(mean_relative_abundance_percent)
    ) %>%
    mutate(
      rank = row_number(),
      reporting_value = paste0(
        Phylum,
        " = ",
        sprintf("%.2f", mean_relative_abundance_percent),
        "%"
      )
    ) %>%
    select(
      Community,
      rank,
      Phylum,
      mean_relative_abundance_percent,
      median_relative_abundance_percent,
      prevalence_percent,
      reporting_value
    )
  
  cat("\nTop phyla by Overall mean relative abundance:\n")
  print(
    reporting_table,
    n = Inf
  )
  
  # -----------------------------
  # Check bar totals
  # -----------------------------
  
  bar_totals <- phylum_summary_top %>%
    group_by(
      Community,
      Soil_habitat
    ) %>%
    summarise(
      total_percent = sum(mean_relative_abundance_percent),
      .groups = "drop"
    )
  
  cat("\nBar totals after collapsing to top phyla + Other:\n")
  print(bar_totals)
  
  return(
    list(
      sample_phylum_all = sample_phylum,
      phylum_summary_all = phylum_summary_all,
      phylum_summary_top = phylum_summary_top,
      reporting_table = reporting_table,
      bar_totals = bar_totals,
      top_phyla = top_phyla
    )
  )
}


# ============================================================
# 4. Run summaries
# ============================================================

prok_phylum_overview <- make_phylum_summary(
  ps_object = ps_prok,
  community_label = "Prokaryotic",
  top_n = 10
)

fungal_phylum_overview <- make_phylum_summary(
  ps_object = ps_fun,
  community_label = "Fungal",
  top_n = 10
)


# ============================================================
# 5. Combine outputs
# ============================================================

phylum_summary_all <- bind_rows(
  prok_phylum_overview$phylum_summary_all,
  fungal_phylum_overview$phylum_summary_all
)

phylum_summary_top <- bind_rows(
  prok_phylum_overview$phylum_summary_top,
  fungal_phylum_overview$phylum_summary_top
)

phylum_reporting_values <- bind_rows(
  prok_phylum_overview$reporting_table,
  fungal_phylum_overview$reporting_table
)

phylum_bar_totals <- bind_rows(
  prok_phylum_overview$bar_totals,
  fungal_phylum_overview$bar_totals
)


# ============================================================
# 6. Save CSV outputs
# ============================================================

write_csv(
  phylum_summary_all,
  "microbial_overview_outputs/phylum_summary_all_taxa_full.csv"
)

write_csv(
  phylum_summary_top,
  "microbial_overview_outputs/phylum_summary_top10_plus_other_for_plotting.csv"
)

write_csv(
  phylum_reporting_values,
  "microbial_overview_outputs/phylum_reporting_values_for_results.csv"
)

write_csv(
  phylum_bar_totals,
  "microbial_overview_outputs/phylum_bar_totals_check.csv"
)


# ============================================================
# 7. Save Excel workbook if available
# ============================================================

if (requireNamespace("writexl", quietly = TRUE)) {
  
  writexl::write_xlsx(
    list(
      README = tibble(
        Note = c(
          "Phylum-level overview generated from final manuscript phyloseq objects: ps_prok and ps_fun.",
          "Relative abundance was calculated within each sample first, then averaged across samples.",
          "No treatment-level comparison is included in this overview.",
          "Top phyla were selected by Overall mean relative abundance within each community.",
          "Non-top phyla were collapsed into Other phyla.",
          "Use phylum_summary_top10_plus_other_for_plotting for the bar plot.",
          "Use phylum_reporting_values_for_results for main-text descriptive reporting."
        )
      ),
      
      Reporting_values = phylum_reporting_values,
      Top10_plus_other = phylum_summary_top,
      All_phyla = phylum_summary_all,
      Bar_total_check = phylum_bar_totals
    ),
    path = "microbial_overview_outputs/Supplementary_Table_S3_phylum_overview_MASTER.xlsx"
  )
}


# ============================================================
# 8. Save complete R object
# ============================================================

saveRDS(
  list(
    prok_phylum_overview = prok_phylum_overview,
    fungal_phylum_overview = fungal_phylum_overview,
    phylum_summary_all = phylum_summary_all,
    phylum_summary_top = phylum_summary_top,
    phylum_reporting_values = phylum_reporting_values,
    phylum_bar_totals = phylum_bar_totals
  ),
  "microbial_overview_outputs/phylum_overview_complete.rds"
)


# ============================================================
# 9. Final print for audit
# ============================================================

cat("\n========================================\n")
cat("FINAL PHYLUM REPORTING VALUES\n")
cat("========================================\n")

print(
  phylum_reporting_values,
  n = Inf
)

cat("\n========================================\n")
cat("BAR TOTALS CHECK\n")
cat("========================================\n")

print(
  phylum_bar_totals
)

cat("\nSaved outputs in folder:\n")
cat("microbial_overview_outputs/\n")


# ============================================================
# EXPLORE PHYLUM → CLASS → ORDER COMPOSITION
# Objects:
#   ps_prok
#   ps_fun
#
# No graph created.
# Output:
#   Mean relative abundance of each Order
#   Percent of each Order within its Phylum
#   Cumulative percent within Phylum
# ============================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)

dir.create(
  "microbial_overview_outputs",
  showWarnings = FALSE
)

# ============================================================
# Helper functions
# ============================================================

clean_tax_label <- function(x, rank_name) {
  
  x <- as.character(x)
  
  x <- case_when(
    is.na(x) ~ paste0("Unclassified_", rank_name),
    x == "" ~ paste0("Unclassified_", rank_name),
    x %in% c("NA", "na", "Unclassified", "unclassified", "unidentified", "Unidentified") ~ paste0("Unclassified_", rank_name),
    TRUE ~ x
  )
  
  x %>%
    str_replace("^d__", "") %>%
    str_replace("^k__", "") %>%
    str_replace("^p__", "") %>%
    str_replace("^c__", "") %>%
    str_replace("^o__", "") %>%
    str_replace("^f__", "") %>%
    str_replace("^g__", "") %>%
    str_replace("^s__", "")
}


find_rank_col <- function(tax_names, possible_names) {
  
  hit <- intersect(
    possible_names,
    tax_names
  )
  
  if (length(hit) == 0) {
    stop(
      "Could not find any of these taxonomy columns: ",
      paste(possible_names, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(tax_names, collapse = ", ")
    )
  }
  
  hit[1]
}


make_phylum_order_summary <- function(ps_object,
                                      community_label,
                                      phylum_cutoff_percent = 1,
                                      order_cutoff_percent = 0.5,
                                      max_orders_per_phylum = 12) {
  
  cat("\n========================================\n")
  cat("OBJECT:", community_label, "\n")
  cat("========================================\n")
  
  cat("Samples:", phyloseq::nsamples(ps_object), "\n")
  cat("ASVs:", phyloseq::ntaxa(ps_object), "\n")
  cat("Reads:", sum(phyloseq::sample_sums(ps_object)), "\n")
  
  # -----------------------------
  # Taxonomy
  # -----------------------------
  
  tax <- data.frame(
    phyloseq::tax_table(ps_object)
  )
  
  tax$ASV <- rownames(tax)
  
  phylum_col <- find_rank_col(
    names(tax),
    c("Phylum", "phylum")
  )
  
  class_col <- find_rank_col(
    names(tax),
    c("Class", "class")
  )
  
  order_col <- find_rank_col(
    names(tax),
    c("Order", "order")
  )
  
  tax_clean <- tax %>%
    transmute(
      ASV = ASV,
      Phylum = clean_tax_label(.data[[phylum_col]], "Phylum"),
      Class = clean_tax_label(.data[[class_col]], "Class"),
      Order = clean_tax_label(.data[[order_col]], "Order")
    )
  
  # -----------------------------
  # Relative abundance
  # -----------------------------
  
  ps_rel <- phyloseq::transform_sample_counts(
    ps_object,
    function(x) x / sum(x)
  )
  
  otu <- as(
    phyloseq::otu_table(ps_rel),
    "matrix"
  )
  
  if (phyloseq::taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  
  otu_long <- as.data.frame(otu) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(
      cols = -SampleID,
      names_to = "ASV",
      values_to = "Relative_abundance"
    )
  
  # -----------------------------
  # Collapse ASVs to Phylum-Class-Order per sample
  # -----------------------------
  
  sample_order <- otu_long %>%
    left_join(
      tax_clean,
      by = "ASV"
    ) %>%
    group_by(
      SampleID,
      Phylum,
      Class,
      Order
    ) %>%
    summarise(
      Relative_abundance = sum(Relative_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  # -----------------------------
  # Overall order summary
  # -----------------------------
  
  order_summary <- sample_order %>%
    group_by(
      Community = community_label,
      Phylum,
      Class,
      Order
    ) %>%
    summarise(
      Mean_RA = mean(Relative_abundance, na.rm = TRUE) * 100,
      Median_RA = median(Relative_abundance, na.rm = TRUE) * 100,
      Prevalence_percent = mean(Relative_abundance > 0, na.rm = TRUE) * 100,
      n_samples = n(),
      .groups = "drop"
    )
  
  phylum_summary <- order_summary %>%
    group_by(
      Community,
      Phylum
    ) %>%
    summarise(
      Phylum_Mean_RA = sum(Mean_RA),
      .groups = "drop"
    ) %>%
    arrange(
      desc(Phylum_Mean_RA)
    )
  
  retained_phyla <- phylum_summary %>%
    filter(
      Phylum_Mean_RA >= phylum_cutoff_percent
    ) %>%
    pull(Phylum)
  
  order_summary_within_phylum <- order_summary %>%
    left_join(
      phylum_summary,
      by = c("Community", "Phylum")
    ) %>%
    mutate(
      Percent_within_phylum =
        100 * Mean_RA / Phylum_Mean_RA
    ) %>%
    group_by(
      Community,
      Phylum
    ) %>%
    arrange(
      desc(Mean_RA),
      .by_group = TRUE
    ) %>%
    mutate(
      Order_rank_within_phylum = row_number(),
      Cumulative_within_phylum = cumsum(Percent_within_phylum),
      Status = case_when(
        !Phylum %in% retained_phyla ~ "Collapsed into Other phyla",
        Mean_RA >= order_cutoff_percent ~ "Retained",
        Order_rank_within_phylum <= max_orders_per_phylum ~ "Retained",
        TRUE ~ "Collapsed into Other orders"
      )
    ) %>%
    ungroup()
  
  # -----------------------------
  # Collapse small orders for future plotting
  # -----------------------------
  
  order_plot_table <- order_summary_within_phylum %>%
    filter(
      Phylum %in% retained_phyla
    ) %>%
    mutate(
      Order_plot = if_else(
        Status == "Retained",
        Order,
        "Other orders"
      )
    ) %>%
    group_by(
      Community,
      Phylum,
      Class,
      Order_plot
    ) %>%
    summarise(
      Mean_RA = sum(Mean_RA, na.rm = TRUE),
      Median_RA = sum(Median_RA, na.rm = TRUE),
      Prevalence_percent = max(Prevalence_percent, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      phylum_summary,
      by = c("Community", "Phylum")
    ) %>%
    mutate(
      Percent_within_phylum =
        100 * Mean_RA / Phylum_Mean_RA
    ) %>%
    group_by(
      Community,
      Phylum
    ) %>%
    arrange(
      desc(Mean_RA),
      .by_group = TRUE
    ) %>%
    mutate(
      Cumulative_within_phylum = cumsum(Percent_within_phylum)
    ) %>%
    ungroup()
  
  # -----------------------------
  # Print Phylum → Class → Order
  # -----------------------------
  
  cat("\n========================================\n")
  cat("PHYLUM SUMMARY:", community_label, "\n")
  cat("========================================\n")
  
  print(
    phylum_summary,
    n = Inf
  )
  
  cat("\n========================================\n")
  cat("PHYLUM → CLASS → ORDER:", community_label, "\n")
  cat("Only phyla >= ", phylum_cutoff_percent, "% mean RA printed\n", sep = "")
  cat("Orders retained if >= ", order_cutoff_percent, "% mean RA or among top ",
      max_orders_per_phylum, " within phylum\n", sep = "")
  cat("========================================\n")
  
  for (phy in retained_phyla) {
    
    cat("\n----------------------------------------\n")
    cat(phy, "\n")
    cat("----------------------------------------\n")
    
    tmp <- order_summary_within_phylum %>%
      filter(
        Phylum == phy
      ) %>%
      arrange(
        desc(Mean_RA)
      ) %>%
      select(
        Class,
        Order,
        Mean_RA,
        Percent_within_phylum,
        Cumulative_within_phylum,
        Prevalence_percent,
        Status
      )
    
    print(
      tmp,
      n = Inf
    )
  }
  
  list(
    phylum_summary = phylum_summary,
    order_summary_full = order_summary_within_phylum,
    order_plot_table = order_plot_table,
    retained_phyla = retained_phyla
  )
}


# ============================================================
# Run prokaryotic and fungal summaries
# ============================================================

prok_phylum_order <- make_phylum_order_summary(
  ps_object = ps_prok,
  community_label = "Prokaryotic",
  phylum_cutoff_percent = 1,
  order_cutoff_percent = 0.5,
  max_orders_per_phylum = 12
)

fungal_phylum_order <- make_phylum_order_summary(
  ps_object = ps_fun,
  community_label = "Fungal",
  phylum_cutoff_percent = 1,
  order_cutoff_percent = 0.5,
  max_orders_per_phylum = 12
)


# ============================================================
# Combine outputs
# ============================================================

phylum_order_summary_full <- bind_rows(
  prok_phylum_order$order_summary_full,
  fungal_phylum_order$order_summary_full
)

phylum_order_plot_table <- bind_rows(
  prok_phylum_order$order_plot_table,
  fungal_phylum_order$order_plot_table
)

phylum_summary_combined <- bind_rows(
  prok_phylum_order$phylum_summary,
  fungal_phylum_order$phylum_summary
)


# ============================================================
# Save outputs
# ============================================================

write_csv(
  phylum_summary_combined,
  "microbial_overview_outputs/phylum_summary_for_order_exploration.csv"
)

write_csv(
  phylum_order_summary_full,
  "microbial_overview_outputs/phylum_class_order_full_exploration.csv"
)

write_csv(
  phylum_order_plot_table,
  "microbial_overview_outputs/phylum_class_order_retained_for_possible_plot.csv"
)

if (requireNamespace("writexl", quietly = TRUE)) {
  
  writexl::write_xlsx(
    list(
      README = tibble(
        Note = c(
          "Generated from final manuscript objects ps_prok and ps_fun.",
          "Relative abundance was calculated within each sample first, then averaged across samples.",
          "This is an exploration of Phylum → Class → Order composition.",
          "No treatment grouping and no statistical inference are included.",
          "Phyla with overall mean relative abundance >= 1% are printed as dominant phyla.",
          "Orders are retained if mean relative abundance >= 0.5% or among the top 12 within a retained phylum."
        )
      ),
      Phylum_summary = phylum_summary_combined,
      Full_phylum_class_order = phylum_order_summary_full,
      Retained_for_possible_plot = phylum_order_plot_table
    ),
    path = "microbial_overview_outputs/Supplementary_Table_Sx_phylum_class_order_exploration_MASTER.xlsx"
  )
}

saveRDS(
  list(
    prok_phylum_order = prok_phylum_order,
    fungal_phylum_order = fungal_phylum_order,
    phylum_summary_combined = phylum_summary_combined,
    phylum_order_summary_full = phylum_order_summary_full,
    phylum_order_plot_table = phylum_order_plot_table
  ),
  "microbial_overview_outputs/phylum_class_order_exploration_complete.rds"
)

cat("\n========================================\n")
cat("EXPLORATION COMPLETE — NO GRAPH CREATED\n")
cat("========================================\n")

cat("\nSaved:\n")
cat("microbial_overview_outputs/phylum_summary_for_order_exploration.csv\n")
cat("microbial_overview_outputs/phylum_class_order_full_exploration.csv\n")
cat("microbial_overview_outputs/phylum_class_order_retained_for_possible_plot.csv\n")
cat("microbial_overview_outputs/Supplementary_Table_Sx_phylum_class_order_exploration_MASTER.xlsx\n")



# ============================================================
# OVERALL MICROBIAL COMPOSITION BY SITE × TREATMENT
#
# Figure:
#   Rows  = Community: Prokaryotic / Fungal
#   Facets = Site
#   x-axis = Treatment
#   fill = selected Phylum-Class groups
#
# No text inside bars.
# Legend only.
# Uses final manuscript objects:
#   ps_prok
#   ps_fun
#
# Requires:
#   final_phylum_class_selection_audit.csv
# from previous step, or object class_selection_clean in memory.
# ============================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)

dir.create(
  "microbial_overview_outputs",
  showWarnings = FALSE
)

# ============================================================
# 1. Load class-selection audit if needed
# ============================================================

if (!exists("class_selection_clean")) {
  class_selection_clean <- read_csv(
    "microbial_overview_outputs/final_phylum_class_selection_audit.csv",
    show_col_types = FALSE
  )
}

# ============================================================
# 2. Colours
# ============================================================

class_cols <- c(
  # Prokaryotic — Actinobacteriota
  "Prokaryotic__Actinobacteriota__Acidimicrobiia" = "#03045e",
  "Prokaryotic__Actinobacteriota__Actinobacteria" = "#023e8a",
  "Prokaryotic__Actinobacteriota__Thermoleophilia" = "#80ffdb",
  "Prokaryotic__Actinobacteriota__Other classes" = "#00b4d8",
  
  # Prokaryotic — Proteobacteria
  "Prokaryotic__Proteobacteria__Alphaproteobacteria" = "#ff9e00",
  "Prokaryotic__Proteobacteria__Gammaproteobacteria" = "#ff8500",
  "Prokaryotic__Proteobacteria__Other classes" = "#ff5400",
  
  # Prokaryotic — Halobacterota
  "Prokaryotic__Halobacterota__Halobacteria" = "#ff0000",
  "Prokaryotic__Halobacterota__Other classes" = "#e5383b",
  
  # Prokaryotic — Chloroflexi
  "Prokaryotic__Chloroflexi__Anaerolineae" = "#240046",
  "Prokaryotic__Chloroflexi__Chloroflexia" = "#3c096c",
  "Prokaryotic__Chloroflexi__KD4-96" = "#990b52",
  "Prokaryotic__Chloroflexi__Dehalococcoidia" = "#312244",
  "Prokaryotic__Chloroflexi__Other classes" = "#4d194d",
  
  # Prokaryotic — Other phyla
  "Prokaryotic__Other phyla__Other phyla" = "#f20089",
  
  # Prokaryotic — Bacteroidota
  "Prokaryotic__Bacteroidota__Bacteroidia" = "#006466",
  "Prokaryotic__Bacteroidota__Rhodothermia" = "#144552",
  "Prokaryotic__Bacteroidota__Other classes" = "#212f45",
  
  # Prokaryotic — Gemmatimonadota
  "Prokaryotic__Gemmatimonadota__PAUC43f marine benthic group" = "#b4fadc",
  "Prokaryotic__Gemmatimonadota__BD2-11 terrestrial group" = "#4d86a5",
  "Prokaryotic__Gemmatimonadota__Other classes" = "#cae0e4",
  
  # Prokaryotic — Desulfobacterota
  "Prokaryotic__Desulfobacterota__Desulfuromonadia" = "#9cf945",
  "Prokaryotic__Desulfobacterota__Desulfobulbia" = "#8edf34",
  "Prokaryotic__Desulfobacterota__Desulfobacteria" = "#509724",
  "Prokaryotic__Desulfobacterota__Other classes" = "#1f6924",
  
  # Prokaryotic — Firmicutes
  "Prokaryotic__Firmicutes__Bacilli" = "#ffea00",
  "Prokaryotic__Firmicutes__Clostridia" = "#ffc300",
  "Prokaryotic__Firmicutes__Other classes" = "#ffa200",
  
  # Prokaryotic — Acidobacteriota
  "Prokaryotic__Acidobacteriota__Thermoanaerobaculia" = "#131313",
  "Prokaryotic__Acidobacteriota__Other classes" = "#8fb7b0",
  
  # Prokaryotic — Halanaerobiaeota
  "Prokaryotic__Halanaerobiaeota__Halanaerobiia" = "#597380",
  
  # Fungal — Ascomycota
  "Fungal__Ascomycota__Dothideomycetes" = "#03045e",
  "Fungal__Ascomycota__Sordariomycetes" = "#023e8a",
  "Fungal__Ascomycota__Unclassified class" = "#0f80aa",
  "Fungal__Ascomycota__Eurotiomycetes" = "#00b4d8",
  "Fungal__Ascomycota__Leotiomycetes" = "#3fc1c0",
  "Fungal__Ascomycota__Other classes" = "#a1fcdf",
  
  # Fungal — Basidiomycota
  "Fungal__Basidiomycota__Unclassified class" = "#28005d",
  "Fungal__Basidiomycota__Agaricomycetes" = "#4e0057",
  "Fungal__Basidiomycota__Cystobasidiomycetes" = "#941b85",
  "Fungal__Basidiomycota__Tremellomycetes" = "#901447",
  "Fungal__Basidiomycota__Agaricostilbomycetes" = "#cd0c2b",
  "Fungal__Basidiomycota__Other classes" = "#ff1654",
  
  # Fungal — Unclassified / Other
  "Fungal__Unclassified fungi__Unclassified class" = "#8fb7b0",
  "Fungal__Other phyla__Other phyla" = "#3e5168"
)

# ============================================================
# 3. Helper functions
# ============================================================

clean_tax_label <- function(x, rank_name) {
  
  x <- as.character(x)
  
  x <- case_when(
    is.na(x) ~ paste0("Unclassified_", rank_name),
    x == "" ~ paste0("Unclassified_", rank_name),
    x %in% c(
      "NA", "na",
      "Unclassified", "unclassified",
      "unidentified", "Unidentified"
    ) ~ paste0("Unclassified_", rank_name),
    TRUE ~ x
  )
  
  x %>%
    str_replace("^d__", "") %>%
    str_replace("^k__", "") %>%
    str_replace("^p__", "") %>%
    str_replace("^c__", "") %>%
    str_replace("^o__", "") %>%
    str_replace("^f__", "") %>%
    str_replace("^g__", "") %>%
    str_replace("^s__", "")
}

standardise_treatment <- function(x) {
  
  x <- toupper(trimws(as.character(x)))
  
  case_when(
    x %in% c("NAT", "NATURAL") ~ "Natural",
    x %in% c("BESE") ~ "BESE",
    x %in% c("BARE") ~ "BARE",
    x %in% c("DEG", "DEGRADED", "DEGREST") ~ "Degraded",
    TRUE ~ as.character(x)
  )
}

make_grouped_tax_table <- function(ps_object, community_label) {
  
  cat("\n========================================\n")
  cat("COMMUNITY:", community_label, "\n")
  cat("========================================\n")
  
  cat("Samples:", nsamples(ps_object), "\n")
  cat("ASVs:", ntaxa(ps_object), "\n")
  cat("Reads:", sum(sample_sums(ps_object)), "\n")
  
  # -----------------------------
  # Metadata
  # -----------------------------
  
  meta <- data.frame(sample_data(ps_object))
  meta$SampleID <- rownames(meta)
  
  if (!"Site" %in% names(meta)) {
    stop("Column 'Site' not found in sample_data.")
  }
  
  treatment_col <- case_when(
    "Treatment" %in% names(meta) ~ "Treatment",
    "EcologicalStatus" %in% names(meta) ~ "EcologicalStatus",
    "Ecological_condition" %in% names(meta) ~ "Ecological_condition",
    TRUE ~ NA_character_
  )
  
  if (is.na(treatment_col)) {
    stop(
      "No treatment/ecological-condition column found. Available columns:\n",
      paste(names(meta), collapse = ", ")
    )
  }
  
  meta_clean <- meta %>%
    transmute(
      SampleID = SampleID,
      Site = as.character(Site),
      Treatment = standardise_treatment(.data[[treatment_col]])
    )
  
  # -----------------------------
  # Taxonomy
  # -----------------------------
  
  tax <- data.frame(tax_table(ps_object))
  tax$ASV <- rownames(tax)
  
  if (!all(c("Phylum", "Class") %in% names(tax))) {
    stop(
      "Taxonomy table must contain Phylum and Class. Available columns:\n",
      paste(names(tax), collapse = ", ")
    )
  }
  
  tax_clean <- tax %>%
    transmute(
      ASV = ASV,
      Phylum = clean_tax_label(Phylum, "Phylum"),
      Class = clean_tax_label(Class, "Class")
    ) %>%
    mutate(
      Phylum = case_when(
        community_label == "Fungal" &
          Phylum == "Unclassified_Phylum" ~ "Unclassified fungi",
        TRUE ~ Phylum
      ),
      Class = case_when(
        Class == "Unclassified_Class" ~ "Unclassified class",
        TRUE ~ Class
      )
    )
  
  # -----------------------------
  # Class selection mapper
  # -----------------------------
  
  mapper <- class_selection_clean %>%
    filter(Community == community_label) %>%
    select(
      Community,
      Phylum,
      Class,
      Phylum_selected,
      Retain_class_final
    ) %>%
    distinct()
  
  # -----------------------------
  # Relative abundance per sample
  # -----------------------------
  
  ps_rel <- transform_sample_counts(
    ps_object,
    function(x) x / sum(x)
  )
  
  otu <- as(otu_table(ps_rel), "matrix")
  
  if (taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  
  otu_long <- as.data.frame(otu) %>%
    rownames_to_column("SampleID") %>%
    pivot_longer(
      cols = -SampleID,
      names_to = "ASV",
      values_to = "Relative_abundance"
    )
  
  # -----------------------------
  # Collapse to selected Phylum-Class groups
  # -----------------------------
  
  sample_tax_grouped <- otu_long %>%
    left_join(tax_clean, by = "ASV") %>%
    mutate(Community = community_label) %>%
    left_join(
      mapper,
      by = c("Community", "Phylum", "Class")
    ) %>%
    mutate(
      Phylum_selected = if_else(
        is.na(Phylum_selected),
        FALSE,
        Phylum_selected
      ),
      Retain_class_final = if_else(
        is.na(Retain_class_final),
        FALSE,
        Retain_class_final
      ),
      Phylum_plot = case_when(
        Phylum_selected ~ Phylum,
        TRUE ~ "Other phyla"
      ),
      Class_plot = case_when(
        !Phylum_selected ~ "Other phyla",
        Retain_class_final ~ Class,
        TRUE ~ "Other classes"
      ),
      Segment_ID = paste(
        Community,
        Phylum_plot,
        Class_plot,
        sep = "__"
      )
    ) %>%
    left_join(meta_clean, by = "SampleID") %>%
    group_by(
      Community,
      SampleID,
      Site,
      Treatment,
      Phylum_plot,
      Class_plot,
      Segment_ID
    ) %>%
    summarise(
      Relative_abundance = sum(Relative_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  sample_tax_grouped
}

# ============================================================
# 4. Build grouped composition table
# ============================================================

prok_grouped <- make_grouped_tax_table(
  ps_object = ps_prok,
  community_label = "Prokaryotic"
)

fungal_grouped <- make_grouped_tax_table(
  ps_object = ps_fun,
  community_label = "Fungal"
)

sample_tax_grouped <- bind_rows(
  prok_grouped,
  fungal_grouped
)

# ============================================================
# 5. Summarise by Community × Site × Treatment
# ============================================================

plot_data <- sample_tax_grouped %>%
  group_by(
    Community,
    Site,
    Treatment,
    Phylum_plot,
    Class_plot,
    Segment_ID
  ) %>%
  summarise(
    Mean_RA_percent = mean(Relative_abundance, na.rm = TRUE) * 100,
    Median_RA_percent = median(Relative_abundance, na.rm = TRUE) * 100,
    Prevalence_percent = mean(Relative_abundance > 0, na.rm = TRUE) * 100,
    n_samples = n(),
    .groups = "drop"
  )

# ============================================================
# 6. ORDER SEGMENTS BY OVERALL ABUNDANCE
# Bottom to top:
#   most abundant phylum -> least abundant phylum
#   within each phylum: most abundant class -> least abundant class
# ============================================================

treatment_order <- c("Natural", "BESE", "BARE", "Degraded")
site_order <- c("ACP", "BASS", "TRUG")

# Use the final overall table as the reference for ordering
# This avoids alphabetical stacking.
order_reference <- final_phylum_class_table %>%
  mutate(
    Segment_ID = paste(
      Community,
      Phylum_plot,
      Class_plot,
      sep = "__"
    )
  ) %>%
  arrange(
    Community,
    desc(Phylum_total_percent),
    desc(Mean_RA_percent)
  ) %>%
  group_by(
    Community
  ) %>%
  mutate(
    Global_segment_rank = row_number()
  ) %>%
  ungroup() %>%
  select(
    Community,
    Phylum_plot,
    Class_plot,
    Segment_ID,
    Phylum_total_percent,
    Global_segment_rank
  )

# Join ordering onto the site × treatment plot data
plot_data_ordered <- plot_data %>%
  left_join(
    order_reference,
    by = c(
      "Community",
      "Phylum_plot",
      "Class_plot",
      "Segment_ID"
    )
  ) %>%
  mutate(
    Community = factor(
      Community,
      levels = c("Prokaryotic", "Fungal")
    ),
    Site = factor(
      Site,
      levels = site_order
    ),
    Treatment = factor(
      Treatment,
      levels = treatment_order
    ),
    x_pos = as.numeric(Treatment),
    Segment_ID = factor(
      Segment_ID,
      levels = order_reference$Segment_ID
    )
  ) %>%
  arrange(
    Community,
    Site,
    Treatment,
    Global_segment_rank
  ) %>%
  group_by(
    Community,
    Site,
    Treatment
  ) %>%
  mutate(
    ymin = lag(cumsum(Mean_RA_percent), default = 0),
    ymax = cumsum(Mean_RA_percent)
  ) %>%
  ungroup()

# ============================================================
# 7. Check ordered bar totals
# ============================================================

bar_total_check <- plot_data_ordered %>%
  group_by(
    Community,
    Site,
    Treatment
  ) %>%
  summarise(
    total_percent = sum(Mean_RA_percent),
    n_segments = n(),
    .groups = "drop"
  ) %>%
  arrange(
    Community,
    Site,
    Treatment
  )

cat("\n========================================\n")
cat("ORDERED BAR TOTAL CHECK\n")
cat("========================================\n")

print(
  bar_total_check,
  n = Inf
)

# ============================================================
# 8. Legend labels in same order as stacked bars
# ============================================================

legend_labels <- order_reference %>%
  mutate(
    Legend_label = case_when(
      Phylum_plot == Class_plot ~ Phylum_plot,
      TRUE ~ paste0(Phylum_plot, " — ", Class_plot)
    )
  )

legend_breaks <- legend_labels$Segment_ID

legend_names <- legend_labels$Legend_label
names(legend_names) <- legend_breaks

missing_cols <- setdiff(
  unique(as.character(plot_data_ordered$Segment_ID)),
  names(class_cols)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing colours for:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ============================================================
# 9. Plot with manual stacking
# This guarantees correct order.
# No text inside bars.
# ============================================================

bar_width <- 0.72

overall_site_treatment_plot <- ggplot() +
  
  geom_rect(
    data = plot_data_ordered,
    aes(
      xmin = x_pos - bar_width / 2,
      xmax = x_pos + bar_width / 2,
      ymin = ymin,
      ymax = ymax,
      fill = Segment_ID
    ),
    colour = "black",
    linewidth = 0.10
  ) +
  
  facet_grid(
    Community ~ Site,
    scales = "free_x",
    space = "free_x"
  ) +
  
  scale_x_continuous(
    breaks = seq_along(treatment_order),
    labels = treatment_order,
    expand = c(0.04, 0.04)
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%"),
    expand = c(0, 0)
  ) +
  
  scale_fill_manual(
    values = class_cols,
    breaks = legend_breaks,
    labels = legend_names,
    name = "Phylum — class"
  ) +
  
  labs(
    x = NULL,
    y = "Mean relative abundance",
    title = "Overall microbial composition by site and ecological condition",
    subtitle = "Class-level groups nested within dominant phyla; bulk-soil and rhizosphere samples combined"
  ) +
  
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = FALSE,
      keyheight = unit(0.36, "cm"),
      keywidth = unit(0.36, "cm")
    )
  ) +
  
  coord_cartesian(
    clip = "off"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0,
      size = 13
    ),
    plot.subtitle = element_text(
      hjust = 0,
      size = 10
    ),
    strip.background = element_rect(
      fill = "grey92",
      colour = "grey40",
      linewidth = 0.3
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 9
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.title.y = element_text(
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(
      size = 7.0
    ),
    legend.box.spacing = unit(0.15, "cm"),
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )

overall_site_treatment_plot

# ============================================================
# 10. Save ordered figure and data
# ============================================================

ggsave(
  filename = "microbial_overview_outputs/Figure_X_overall_composition_site_treatment_ORDERED.png",
  plot = overall_site_treatment_plot,
  width = 13,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "microbial_overview_outputs/Figure_X_overall_composition_site_treatment_ORDERED.pdf",
  plot = overall_site_treatment_plot,
  width = 13,
  height = 7.5,
  units = "in"
)

write_csv(
  plot_data_ordered,
  "microbial_overview_outputs/Figure_X_overall_composition_site_treatment_ORDERED_plot_data.csv"
)

write_csv(
  bar_total_check,
  "microbial_overview_outputs/Figure_X_overall_composition_site_treatment_ORDERED_bar_total_check.csv"
)

write_csv(
  legend_labels,
  "microbial_overview_outputs/Figure_X_overall_composition_site_treatment_ORDERED_legend_labels.csv"
)

cat("\nSaved ordered figure:\n")
cat("microbial_overview_outputs/Figure_X_overall_composition_site_treatment_ORDERED.png\n")
cat("microbial_overview_outputs/Figure_X_overall_composition_site_treatment_ORDERED.pdf\n")
