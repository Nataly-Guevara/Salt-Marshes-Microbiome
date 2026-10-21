# ============================================================
# R3.4 / M6 PLANT-CONDITION COMMUNITY TEST
#
# Question:
#   Does visible plant condition explain whole rhizosphere
#   microbial community composition in BESE and BARE plots?
#
# Test:
#   Bray-Curtis rhizosphere community composition ~ Site + Treatment + Green foliage
#
# Communities:
#   1. Prokaryotic rhizosphere
#   2. Fungal rhizosphere
#
# Required:
#   ps_prok
#   ps_fun
#   plant_condition object OR plant_condition.csv file
#
# Outputs:
#   plant_condition_integration_outputs/
# ============================================================

library(phyloseq)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(vegan)
library(permute)

set.seed(123)

dir.create(
  "plant_condition_integration_outputs",
  showWarnings = FALSE
)

# ============================================================
# 1. Helper functions
# ============================================================

find_col <- function(df, candidates) {
  exact <- candidates[candidates %in% names(df)]
  if (length(exact) > 0) return(exact[1])
  
  lower_names <- tolower(names(df))
  lower_candidates <- tolower(candidates)
  hit <- match(lower_candidates, lower_names)
  hit <- hit[!is.na(hit)]
  
  if (length(hit) > 0) return(names(df)[hit[1]])
  
  stop(
    "Could not find any of these columns:\n",
    paste(candidates, collapse = ", "),
    "\n\nAvailable columns:\n",
    paste(names(df), collapse = ", ")
  )
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

standardise_habitat <- function(x) {
  x <- toupper(trimws(as.character(x)))
  
  case_when(
    x %in% c("BS", "BULK", "BULK SOIL", "BULK_SOIL") ~ "BS",
    x %in% c("RZ", "RHIZOSPHERE") ~ "RZ",
    TRUE ~ as.character(x)
  )
}

clean_block <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace("^Block", "") %>%
    str_replace("^BLOCK", "") %>%
    str_replace("^block", "") %>%
    str_replace("^B", "") %>%
    str_replace("^b", "") %>%
    str_trim()
}

anova_to_tibble <- function(x, analysis_name, test_type) {
  out <- as.data.frame(x)
  out$Term <- rownames(out)
  rownames(out) <- NULL
  
  out %>%
    as_tibble() %>%
    relocate(Term) %>%
    mutate(
      Analysis = analysis_name,
      Test_type = test_type,
      .before = 1
    )
}

get_sample_metadata <- function(ps_object) {
  meta <- data.frame(sample_data(ps_object))
  meta$SampleID <- rownames(meta)
  
  site_col <- find_col(
    meta,
    c("Site", "SITE", "site")
  )
  
  treatment_col <- find_col(
    meta,
    c(
      "Treatment",
      "TREATMENT",
      "EcologicalStatus",
      "Ecological_condition",
      "EcologicalCondition"
    )
  )
  
  habitat_col <- find_col(
    meta,
    c(
      "SoilPortion",
      "Soil_habitat",
      "SoilHabitat",
      "Compartment",
      "Habitat"
    )
  )
  
  block_col <- find_col(
    meta,
    c(
      "Block",
      "BLOCK",
      "block",
      "Plot",
      "PLOT",
      "plot",
      "PlotID",
      "Plot_ID"
    )
  )
  
  meta %>%
    mutate(
      SampleID = SampleID,
      Site_clean = as.character(.data[[site_col]]),
      Treatment_clean = standardise_treatment(.data[[treatment_col]]),
      Soil_habitat_clean = standardise_habitat(.data[[habitat_col]]),
      Block_clean = clean_block(.data[[block_col]])
    )
}

get_relative_abundance_matrix <- function(ps_object, sample_order) {
  ps_rel <- transform_sample_counts(
    ps_object,
    function(x) x / sum(x)
  )
  
  otu <- as(otu_table(ps_rel), "matrix")
  
  if (taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  
  otu <- otu[sample_order, , drop = FALSE]
  otu <- otu[, colSums(otu, na.rm = TRUE) > 0, drop = FALSE]
  
  otu
}

# ============================================================
# 2. Load plant-condition data
# ============================================================

load_plant_condition <- function() {
  
  if (exists("plant_condition")) {
    message("Using plant_condition object already in memory.")
    return(plant_condition)
  }
  
  candidate_files <- c(
    "plant_condition.csv",
    "plant_condition_for_microbiome.csv",
    "vegetation_outputs/plant_condition.csv",
    "vegetation_outputs/plant_condition_for_microbiome.csv",
    "plant_condition_integration_outputs/Sx_M6_plant_condition_cleaned.csv"
  )
  
  existing <- candidate_files[file.exists(candidate_files)]
  
  if (length(existing) == 0) {
    stop(
      "No plant_condition object or CSV file found.\n",
      "Expected one of:\n",
      paste(candidate_files, collapse = "\n")
    )
  }
  
  message("Loading plant-condition file: ", existing[1])
  read_csv(existing[1], show_col_types = FALSE)
}

prepare_plant_condition <- function(plant_df) {
  
  site_col <- find_col(
    plant_df,
    c("SITE", "Site", "site")
  )
  
  treatment_col <- find_col(
    plant_df,
    c("TREATMENT", "Treatment", "treatment")
  )
  
  block_col <- find_col(
    plant_df,
    c(
      "Block",
      "BLOCK",
      "block",
      "Plot",
      "PLOT",
      "plot",
      "PlotID",
      "Plot_ID"
    )
  )
  
  green_col <- find_col(
    plant_df,
    c(
      "Green_foliage",
      "green_foliage",
      "GREEN_foliage",
      "Green_percent",
      "green_percent",
      "GREEN_percent",
      "Green",
      "GREEN",
      "green"
    )
  )
  
  plant_clean <- plant_df %>%
    transmute(
      Site = as.character(.data[[site_col]]),
      Treatment = standardise_treatment(.data[[treatment_col]]),
      Block_clean = clean_block(.data[[block_col]]),
      Green_foliage = as.numeric(.data[[green_col]])
    ) %>%
    filter(
      Treatment %in% c("BESE", "BARE")
    ) %>%
    mutate(
      Green_foliage = if_else(
        Green_foliage <= 1,
        Green_foliage * 100,
        Green_foliage
      )
    )
  
  cat("\n========================================\n")
  cat("PLANT-CONDITION TABLE CHECK\n")
  cat("========================================\n")
  
  print(
    plant_clean %>%
      count(Site, Treatment, name = "n_plots"),
    n = Inf
  )
  
  plant_clean
}

# ============================================================
# 3. Prepare microbiome + plant-condition data
# ============================================================

prepare_M6_microbiome <- function(ps_object,
                                  community_name,
                                  plant_clean) {
  
  meta <- get_sample_metadata(ps_object)
  
  meta_rz <- meta %>%
    transmute(
      SampleID = SampleID,
      Site = Site_clean,
      Treatment = Treatment_clean,
      Soil_habitat = Soil_habitat_clean,
      Block_clean = Block_clean
    ) %>%
    filter(
      Soil_habitat == "RZ",
      Treatment %in% c("BESE", "BARE")
    ) %>%
    left_join(
      plant_clean,
      by = c("Site", "Treatment", "Block_clean")
    ) %>%
    drop_na(Green_foliage) %>%
    mutate(
      Community = community_name,
      Green_z = as.numeric(scale(Green_foliage))
    )
  
  ps_sub <- prune_samples(meta_rz$SampleID, ps_object)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  comm <- get_relative_abundance_matrix(
    ps_object = ps_sub,
    sample_order = meta_rz$SampleID
  )
  
  list(
    meta = meta_rz,
    comm = comm
  )
}

# ============================================================
# 4. Run PERMANOVA plant-condition test
# ============================================================

run_M6_plant_condition <- function(ps_object,
                                   community_name,
                                   plant_clean) {
  
  prepared <- prepare_M6_microbiome(
    ps_object = ps_object,
    community_name = community_name,
    plant_clean = plant_clean
  )
  
  meta_m6 <- prepared$meta
  comm <- prepared$comm
  
  cat("\n========================================\n")
  cat("M6 PLANT-CONDITION COMMUNITY TEST:", community_name, "\n")
  cat("Samples:", nrow(meta_m6), "\n")
  cat("ASVs:", ncol(comm), "\n")
  cat("========================================\n")
  
  group_counts <- meta_m6 %>%
    count(
      Community,
      Site,
      Treatment,
      name = "n_samples"
    )
  
  print(group_counts, n = Inf)
  
  bray <- vegan::vegdist(
    comm,
    method = "bray"
  )
  
  # Sequential terms.
  # Site is included first to account for spatial structure.
  permanova_sequential <- vegan::adonis2(
    bray ~ Site + Treatment + Green_z,
    data = meta_m6,
    permutations = 9999,
    strata = meta_m6$Site
  )
  
  # Marginal terms.
  # Green_z here is the key reviewer-relevant test:
  # does green foliage explain composition after Site and Treatment?
  permanova_marginal <- vegan::adonis2(
    bray ~ Site + Treatment + Green_z,
    data = meta_m6,
    permutations = 9999,
    strata = meta_m6$Site,
    by = "margin"
  )
  
  seq_table <- anova_to_tibble(
    permanova_sequential,
    analysis_name = paste0("M6 ", community_name),
    test_type = "sequential_terms"
  )
  
  marginal_table <- anova_to_tibble(
    permanova_marginal,
    analysis_name = paste0("M6 ", community_name),
    test_type = "marginal_terms"
  ) %>%
    mutate(
      p_adjust_BH = if_else(
        Term == "Residual",
        NA_real_,
        p.adjust(`Pr(>F)`, method = "BH")
      ),
      Significant_BH_0.05 = if_else(
        is.na(p_adjust_BH),
        NA,
        p_adjust_BH < 0.05
      )
    )
  
  # Optional treatment dispersion check.
  # Useful if Treatment is significant.
  dispersion <- vegan::betadisper(
    bray,
    group = meta_m6$Treatment
  )
  
  dispersion_test <- vegan::permutest(
    dispersion,
    permutations = 9999,
    pairwise = TRUE
  )
  
  dispersion_table <- as.data.frame(dispersion_test$tab)
  dispersion_table$Term <- rownames(dispersion_table)
  rownames(dispersion_table) <- NULL
  
  dispersion_table <- dispersion_table %>%
    as_tibble() %>%
    mutate(
      Analysis = paste0("M6 ", community_name),
      Test_type = "PERMDISP_Treatment",
      .before = 1
    )
  
  list(
    metadata_used = meta_m6,
    group_counts = group_counts,
    permanova_sequential = seq_table,
    permanova_marginal = marginal_table,
    dispersion_treatment = dispersion_table
  )
}

# ============================================================
# 5. Run analyses
# ============================================================

plant_raw <- load_plant_condition()

plant_clean <- prepare_plant_condition(
  plant_raw
)

write_csv(
  plant_clean,
  "plant_condition_integration_outputs/Sx_M6_plant_condition_cleaned.csv"
)

M6_results <- list(
  Prokaryotic = run_M6_plant_condition(
    ps_object = ps_prok,
    community_name = "Prokaryotic",
    plant_clean = plant_clean
  ),
  
  Fungal = run_M6_plant_condition(
    ps_object = ps_fun,
    community_name = "Fungal",
    plant_clean = plant_clean
  )
)

# ============================================================
# 6. Combine outputs
# ============================================================

M6_metadata_used <- bind_rows(
  lapply(M6_results, `[[`, "metadata_used")
)

M6_group_counts <- bind_rows(
  lapply(M6_results, `[[`, "group_counts")
)

M6_permanova_sequential <- bind_rows(
  lapply(M6_results, `[[`, "permanova_sequential")
)

M6_permanova_marginal <- bind_rows(
  lapply(M6_results, `[[`, "permanova_marginal")
)

M6_dispersion_treatment <- bind_rows(
  lapply(M6_results, `[[`, "dispersion_treatment")
)

# ============================================================
# 7. Save supplementary-ready CSV tables
# ============================================================

write_csv(
  M6_metadata_used,
  "plant_condition_integration_outputs/Sx_M6_metadata_used.csv"
)

write_csv(
  M6_group_counts,
  "plant_condition_integration_outputs/Sx_M6_group_counts.csv"
)

write_csv(
  M6_permanova_sequential,
  "plant_condition_integration_outputs/Sx_M6_rhizosphere_composition_PERMANOVA_sequential.csv"
)

write_csv(
  M6_permanova_marginal,
  "plant_condition_integration_outputs/Sx_M6_rhizosphere_composition_PERMANOVA_marginal_BH.csv"
)

write_csv(
  M6_dispersion_treatment,
  "plant_condition_integration_outputs/Sx_M6_rhizosphere_composition_PERMDISP_treatment.csv"
)

# ============================================================
# 8. Save Excel workbook with separate tabs
# ============================================================

if (requireNamespace("openxlsx", quietly = TRUE)) {
  
  wb <- openxlsx::createWorkbook()
  
  table_list <- list(
    README = tibble(
      Note = c(
        "R3.4 / M6 plant-condition community test.",
        "Analysis restricted to BESE and BARE rhizosphere samples.",
        "Green foliage was used as the primary visible plant-condition variable.",
        "PERMANOVA tested Bray-Curtis rhizosphere community composition as a function of Site, Treatment and Green foliage.",
        "Permutations were restricted within Site.",
        "The key reviewer-relevant test is the marginal Green_z term after accounting for Site and Treatment.",
        "PERMDISP by Treatment is included as a dispersion diagnostic."
      )
    ),
    Plant_condition_cleaned = plant_clean,
    Metadata_used = M6_metadata_used,
    Group_counts = M6_group_counts,
    PERMANOVA_sequential = M6_permanova_sequential,
    PERMANOVA_marginal_BH = M6_permanova_marginal,
    PERMDISP_treatment = M6_dispersion_treatment
  )
  
  for (nm in names(table_list)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, table_list[[nm]])
    openxlsx::setColWidths(
      wb,
      nm,
      cols = 1:ncol(table_list[[nm]]),
      widths = "auto"
    )
  }
  
  openxlsx::saveWorkbook(
    wb,
    "plant_condition_integration_outputs/Sx_M6_plant_condition_community_tables.xlsx",
    overwrite = TRUE
  )
}

# ============================================================
# 9. Print key outputs
# ============================================================

cat("\n========================================\n")
cat("M6 GROUP COUNTS\n")
cat("========================================\n")
print(M6_group_counts, n = Inf, width = Inf)

cat("\n========================================\n")
cat("M6 PERMANOVA SEQUENTIAL\n")
cat("========================================\n")
print(M6_permanova_sequential, n = Inf, width = Inf)

cat("\n========================================\n")
cat("M6 PERMANOVA MARGINAL + BH\n")
cat("========================================\n")
print(M6_permanova_marginal, n = Inf, width = Inf)

cat("\n========================================\n")
cat("M6 PERMDISP TREATMENT\n")
cat("========================================\n")
print(M6_dispersion_treatment, n = Inf, width = Inf)

cat("\nDONE. M6 plant-condition community tables saved in:\n")
cat("plant_condition_integration_outputs/\n")