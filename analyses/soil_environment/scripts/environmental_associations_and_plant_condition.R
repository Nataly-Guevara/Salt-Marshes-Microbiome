# ============================================================
# ENVIRONMENTAL ASSOCIATIONS + M6 PLANT-CONDITION TEST
# Chapter 1 revision
#
# Required objects:
#   ps_prok
#   ps_fun
#
# Optional object for M6:
#   plant_condition
#
# Outputs: 
#   environmental_association_outputs/
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

dir.create("environmental_association_outputs", showWarnings = FALSE)
dir.create("plant_condition_integration_outputs", showWarnings = FALSE)

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
    "\nAvailable columns:\n",
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
    str_replace("^B", "") %>%
    str_trim()
}

get_sample_metadata <- function(ps_object) {
  meta <- data.frame(sample_data(ps_object))
  meta$SampleID <- rownames(meta)
  
  site_col <- find_col(meta, c("Site", "SITE", "site"))
  treatment_col <- find_col(
    meta,
    c("Treatment", "TREATMENT", "EcologicalStatus", "Ecological_condition")
  )
  habitat_col <- find_col(
    meta,
    c("SoilPortion", "Soil_habitat", "SoilHabitat", "Compartment", "Habitat")
  )
  
  meta %>%
    mutate(
      SampleID = SampleID,
      Site_clean = as.character(.data[[site_col]]),
      Treatment_clean = standardise_treatment(.data[[treatment_col]]),
      Soil_habitat_clean = standardise_habitat(.data[[habitat_col]])
    )
}

get_relative_abundance_matrix <- function(ps_object, sample_order) {
  ps_rel <- transform_sample_counts(ps_object, function(x) x / sum(x))
  
  otu <- as(otu_table(ps_rel), "matrix")
  
  if (taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  
  otu <- otu[sample_order, , drop = FALSE]
  otu <- otu[, colSums(otu, na.rm = TRUE) > 0, drop = FALSE]
  
  otu
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

make_env_correlation_table <- function(meta_env, analysis_name) {
  env_vars <- c(
    "Moisture_z",
    "Temperature_z",
    "pH_z",
    "NH4_log_z",
    "PO4_log_z",
    "Redox_z"
  )
  
  var_labels <- c(
    Moisture_z = "Moisture",
    Temperature_z = "Temperature",
    pH_z = "pH",
    NH4_log_z = "log NH4",
    PO4_log_z = "log PO4",
    Redox_z = "Redox"
  )
  
  pairs <- combn(env_vars, 2, simplify = FALSE)
  
  bind_rows(lapply(pairs, function(pair) {
    test <- suppressWarnings(
      cor.test(
        meta_env[[pair[1]]],
        meta_env[[pair[2]]],
        method = "spearman",
        exact = FALSE
      )
    )
    
    tibble(
      Analysis = analysis_name,
      Variable_1 = var_labels[[pair[1]]],
      Variable_2 = var_labels[[pair[2]]],
      Spearman_rho = unname(test$estimate),
      p_value = test$p.value
    )
  })) %>%
    mutate(
      p_adjust_BH = p.adjust(p_value, method = "BH"),
      High_correlation_flag = abs(Spearman_rho) >= 0.70
    )
}

prepare_environmental_data <- function(ps_object,
                                       community_name,
                                       habitat_keep,
                                       analysis_name) {
  
  meta <- get_sample_metadata(ps_object)
  
  env_cols <- c(
    Moisture = find_col(meta, c("Moisture_avg", "Moisture", "Soil_moisture")),
    Temperature = find_col(meta, c("Temperature_avg", "Temperature", "Temp")),
    pH = find_col(meta, c("pH_avg", "pH", "PH")),
    NH4 = find_col(meta, c("NH4", "Ammonium", "NH4_avg")),
    PO4 = find_col(meta, c("PO4", "Phosphate", "PO4_avg")),
    Redox = find_col(meta, c("Redox_avg", "Redox", "Eh"))
  )
  
  meta_env <- meta %>%
    filter(Soil_habitat_clean == habitat_keep) %>%
    transmute(
      SampleID = SampleID,
      Community = community_name,
      Analysis = analysis_name,
      Site = Site_clean,
      Treatment = Treatment_clean,
      Soil_habitat = Soil_habitat_clean,
      Moisture = as.numeric(.data[[env_cols["Moisture"]]]),
      Temperature = as.numeric(.data[[env_cols["Temperature"]]]),
      pH = as.numeric(.data[[env_cols["pH"]]]),
      NH4 = as.numeric(.data[[env_cols["NH4"]]]),
      PO4 = as.numeric(.data[[env_cols["PO4"]]]),
      Redox = as.numeric(.data[[env_cols["Redox"]]])
    )
  
  missingness <- meta_env %>%
    summarise(
      n_samples_initial = n(),
      across(
        c(Moisture, Temperature, pH, NH4, PO4, Redox),
        ~ sum(is.na(.x)),
        .names = "missing_{.col}"
      )
    ) %>%
    pivot_longer(
      cols = starts_with("missing_"),
      names_to = "Variable",
      values_to = "n_missing"
    ) %>%
    mutate(
      Variable = str_replace(Variable, "^missing_", ""),
      Analysis = analysis_name,
      Community = community_name,
      Soil_habitat = habitat_keep,
      .before = 1
    )
  
  meta_env_complete <- meta_env %>%
    drop_na(Moisture, Temperature, pH, NH4, PO4, Redox)
  
  if (nrow(meta_env_complete) < 6) {
    stop("Too few complete samples for ", analysis_name)
  }
  
  nh4_transform <- if (any(meta_env_complete$NH4 <= 0, na.rm = TRUE)) {
    "log1p"
  } else {
    "log"
  }
  
  po4_transform <- if (any(meta_env_complete$PO4 <= 0, na.rm = TRUE)) {
    "log1p"
  } else {
    "log"
  }
  
  meta_env_complete <- meta_env_complete %>%
    mutate(
      NH4_log = if (nh4_transform == "log1p") log1p(NH4) else log(NH4),
      PO4_log = if (po4_transform == "log1p") log1p(PO4) else log(PO4),
      Moisture_z = as.numeric(scale(Moisture)),
      Temperature_z = as.numeric(scale(Temperature)),
      pH_z = as.numeric(scale(pH)),
      NH4_log_z = as.numeric(scale(NH4_log)),
      PO4_log_z = as.numeric(scale(PO4_log)),
      Redox_z = as.numeric(scale(Redox))
    )
  
  range_table <- meta_env_complete %>%
    summarise(
      n_samples = n(),
      across(
        c(Moisture, Temperature, pH, NH4, PO4, Redox),
        list(
          min = ~ min(.x, na.rm = TRUE),
          median = ~ median(.x, na.rm = TRUE),
          mean = ~ mean(.x, na.rm = TRUE),
          max = ~ max(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      )
    ) %>%
    pivot_longer(
      cols = -n_samples,
      names_to = "Variable_stat",
      values_to = "Value"
    ) %>%
    separate(
      Variable_stat,
      into = c("Variable", "Statistic"),
      sep = "_(?=[^_]+$)"
    ) %>%
    mutate(
      Analysis = analysis_name,
      Community = community_name,
      Soil_habitat = habitat_keep,
      NH4_transform = nh4_transform,
      PO4_transform = po4_transform,
      .before = 1
    )
  
  ps_sub <- prune_samples(meta_env_complete$SampleID, ps_object)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  comm <- get_relative_abundance_matrix(
    ps_object = ps_sub,
    sample_order = meta_env_complete$SampleID
  )
  
  list(
    meta = meta_env_complete,
    comm = comm,
    missingness = missingness,
    ranges = range_table
  )
}

# ============================================================
# Corrected environmental db-RDA function
# Fixes recode() masking error
# ============================================================

run_environmental_dbRDA <- function(ps_object,
                                    community_name,
                                    habitat_keep,
                                    analysis_name) {
  
  prepared <- prepare_environmental_data(
    ps_object = ps_object,
    community_name = community_name,
    habitat_keep = habitat_keep,
    analysis_name = analysis_name
  )
  
  meta_env <- prepared$meta
  comm <- prepared$comm
  
  cat("\n========================================\n")
  cat("ENVIRONMENTAL ANALYSIS:", analysis_name, "\n")
  cat("Samples:", nrow(meta_env), "\n")
  cat("ASVs:", ncol(comm), "\n")
  cat("========================================\n")
  
  # ----------------------------------------------------------
  # Partial db-RDA
  # Environmental variables fitted after conditioning on Site
  # ----------------------------------------------------------
  
  env_formula <- comm ~
    Moisture_z +
    Temperature_z +
    pH_z +
    NH4_log_z +
    PO4_log_z +
    Redox_z +
    Condition(Site)
  
  db <- vegan::dbrda(
    formula = env_formula,
    data = meta_env,
    distance = "bray",
    add = "lingoes"
  )
  
  overall_test <- vegan::anova.cca(
    db,
    permutations = 9999,
    strata = meta_env$Site
  )
  
  marginal_test <- vegan::anova.cca(
    db,
    by = "margin",
    permutations = 9999,
    strata = meta_env$Site
  )
  
  term_test <- vegan::anova.cca(
    db,
    by = "terms",
    permutations = 9999,
    strata = meta_env$Site
  )
  
  r2 <- vegan::RsquareAdj(db)
  
  r2_table <- tibble(
    Analysis = analysis_name,
    Community = community_name,
    Soil_habitat = habitat_keep,
    n_samples = nrow(meta_env),
    n_ASVs = ncol(comm),
    R2 = unname(r2$r.squared),
    Adjusted_R2 = unname(r2$adj.r.squared)
  )
  
  vif_table <- tibble(
    Analysis = analysis_name,
    Community = community_name,
    Soil_habitat = habitat_keep,
    Variable = names(vegan::vif.cca(db)),
    VIF = as.numeric(vegan::vif.cca(db))
  ) %>%
    mutate(
      High_VIF_flag = VIF >= 10
    )
  
  correlation_table <- make_env_correlation_table(
    meta_env = meta_env,
    analysis_name = analysis_name
  ) %>%
    mutate(
      Community = community_name,
      Soil_habitat = habitat_keep,
      .after = Analysis
    )
  
  high_correlation_table <- correlation_table %>%
    filter(High_correlation_flag)
  
  # ----------------------------------------------------------
  # envfit on Bray-Curtis PCoA
  # This is descriptive, separate from the partial db-RDA test.
  # ----------------------------------------------------------
  
  bray <- vegan::vegdist(comm, method = "bray")
  
  pcoa <- cmdscale(
    bray,
    eig = TRUE,
    k = 2,
    add = TRUE
  )
  
  pcoa_scores <- as.data.frame(pcoa$points)
  colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
  pcoa_scores$SampleID <- rownames(pcoa_scores)
  
  pcoa_scores <- pcoa_scores %>%
    left_join(
      meta_env %>%
        select(
          SampleID,
          Community,
          Analysis,
          Site,
          Treatment,
          Soil_habitat
        ),
      by = "SampleID"
    ) %>%
    mutate(
      PCoA1_percent = round(
        100 * pcoa$eig[1] / sum(pcoa$eig[pcoa$eig > 0]),
        2
      ),
      PCoA2_percent = round(
        100 * pcoa$eig[2] / sum(pcoa$eig[pcoa$eig > 0]),
        2
      )
    )
  
  envfit_vars <- meta_env %>%
    select(
      Moisture_z,
      Temperature_z,
      pH_z,
      NH4_log_z,
      PO4_log_z,
      Redox_z
    )
  
  ef <- vegan::envfit(
    pcoa_scores %>% select(PCoA1, PCoA2),
    envfit_vars,
    permutations = 9999,
    strata = meta_env$Site
  )
  
  ef_vectors <- as.data.frame(
    vegan::scores(ef, display = "vectors")
  )
  
  ef_vectors$Variable_raw <- rownames(ef_vectors)
  rownames(ef_vectors) <- NULL
  
  variable_labels <- c(
    Moisture_z = "Moisture",
    Temperature_z = "Temperature",
    pH_z = "pH",
    NH4_log_z = "log NH4",
    PO4_log_z = "log PO4",
    Redox_z = "Redox"
  )
  
  envfit_table <- ef_vectors %>%
    as_tibble() %>%
    rename(
      Axis1 = PCoA1,
      Axis2 = PCoA2
    ) %>%
    mutate(
      Analysis = analysis_name,
      Community = community_name,
      Soil_habitat = habitat_keep,
      Variable = unname(variable_labels[Variable_raw]),
      r2 = ef$vectors$r,
      p_value = ef$vectors$pvals,
      p_adjust_BH = p.adjust(p_value, method = "BH")
    ) %>%
    select(
      Analysis,
      Community,
      Soil_habitat,
      Variable,
      Variable_raw,
      Axis1,
      Axis2,
      r2,
      p_value,
      p_adjust_BH
    )
  
  list(
    metadata_used = meta_env,
    missingness = prepared$missingness,
    ranges = prepared$ranges,
    correlations = correlation_table,
    high_correlations = high_correlation_table,
    dbrda_r2 = r2_table,
    dbrda_overall = anova_to_tibble(
      overall_test,
      analysis_name,
      "overall"
    ),
    dbrda_terms = anova_to_tibble(
      term_test,
      analysis_name,
      "sequential_terms"
    ),
    dbrda_marginal = anova_to_tibble(
      marginal_test,
      analysis_name,
      "marginal_terms"
    ),
    vif = vif_table,
    envfit = envfit_table,
    pcoa_scores = pcoa_scores
  )
}

env_results <- list(
  Bulk_prokaryotic = run_environmental_dbRDA(
    ps_object = ps_prok,
    community_name = "Prokaryotic",
    habitat_keep = "BS",
    analysis_name = "Bulk soil prokaryotic"
  ),
  
  Bulk_fungal = run_environmental_dbRDA(
    ps_object = ps_fun,
    community_name = "Fungal",
    habitat_keep = "BS",
    analysis_name = "Bulk soil fungal"
  ),
  
  Rhizosphere_prokaryotic = run_environmental_dbRDA(
    ps_object = ps_prok,
    community_name = "Prokaryotic",
    habitat_keep = "RZ",
    analysis_name = "Rhizosphere prokaryotic"
  ),
  
  Rhizosphere_fungal = run_environmental_dbRDA(
    ps_object = ps_fun,
    community_name = "Fungal",
    habitat_keep = "RZ",
    analysis_name = "Rhizosphere fungal"
  )
)


# ============================================================
# Add BH-adjusted p-values to marginal environmental tests
# ============================================================

env_dbrda_marginal_BH <- env_dbrda_marginal %>%
  group_by(Analysis) %>%
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
  ) %>%
  ungroup()

write_csv(
  env_dbrda_marginal_BH,
  "environmental_association_outputs/Sx_environmental_partial_dbRDA_marginal_tests_BH.csv"
)

print(env_dbrda_marginal_BH, n = Inf, width = Inf)

# ============================================================
# SUPPLEMENTARY HEATMAP
# Environmental predictors associated with microbial composition
#
# Input required:
#   env_dbrda_marginal
#
# Output:
#   environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap.png
#   environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap.pdf
#   environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap_data.csv
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)

dir.create(
  "environmental_association_outputs",
  showWarnings = FALSE
)

# ============================================================
# 1. Labels and order
# ============================================================

term_labels <- c(
  Moisture_z = "Moisture",
  Temperature_z = "Temperature",
  pH_z = "pH",
  NH4_log_z = "log NH4",
  PO4_log_z = "log PO4",
  Redox_z = "Redox"
)

analysis_order <- c(
  "Bulk soil prokaryotic",
  "Bulk soil fungal",
  "Rhizosphere prokaryotic",
  "Rhizosphere fungal"
)

variable_order <- c(
  "Moisture",
  "Temperature",
  "pH",
  "log NH4",
  "log PO4",
  "Redox"
)

# ============================================================
# 2. Prepare heatmap data
# ============================================================

heatmap_data <- env_dbrda_marginal %>%
  filter(Term != "Residual") %>%
  mutate(
    Variable = unname(term_labels[Term]),
    
    Community = case_when(
      str_detect(Analysis, "prokaryotic") ~ "Prokaryotic",
      str_detect(Analysis, "fungal") ~ "Fungal",
      TRUE ~ NA_character_
    ),
    
    Soil_habitat = case_when(
      str_detect(Analysis, "Bulk soil") ~ "BS",
      str_detect(Analysis, "Rhizosphere") ~ "RZ",
      TRUE ~ NA_character_
    ),
    
    Analysis = factor(
      Analysis,
      levels = rev(analysis_order)
    ),
    
    Variable = factor(
      Variable,
      levels = variable_order
    )
  ) %>%
  group_by(Analysis) %>%
  mutate(
    p_adjust_BH = p.adjust(`Pr(>F)`, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    Significance = case_when(
      p_adjust_BH < 0.001 ~ "***",
      p_adjust_BH < 0.01  ~ "**",
      p_adjust_BH < 0.05  ~ "*",
      TRUE ~ ""
    ),
    
    F_label = sprintf("%.2f", F),
    
    Cell_label = paste0(
      F_label,
      Significance
    ),
    
    Significant_BH_0.05 = p_adjust_BH < 0.05
  ) %>%
  select(
    Analysis,
    Community,
    Soil_habitat,
    Variable,
    Term,
    Df,
    SumOfSqs,
    F,
    p_value = `Pr(>F)`,
    p_adjust_BH,
    Significant_BH_0.05,
    Significance,
    Cell_label
  )

# Save heatmap data
write_csv(
  heatmap_data,
  "environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap_data.csv"
)

print(
  heatmap_data,
  n = Inf,
  width = Inf
)

# ============================================================
# 3. Heatmap
# Fill = marginal F statistic
# Text = F statistic + BH-adjusted significance stars
# ============================================================

environmental_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = Variable,
    y = Analysis,
    fill = F
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = Cell_label),
    size = 3.5,
    colour = "grey15"
  ) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#3E7281",
    name = "Marginal\nF statistic"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Environmental predictors associated with microbial community composition",
    subtitle = "Partial db-RDA marginal tests after accounting for Site"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0
    ),
    plot.subtitle = element_text(
      size = 10,
      hjust = 0
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "grey20",
      size = 10
    ),
    axis.text.y = element_text(
      colour = "grey20",
      size = 10
    ),
    panel.grid = element_blank(),
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )

environmental_heatmap

# ============================================================
# 4. Save figure
# ============================================================

ggsave(
  filename = "environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap.png",
  plot = environmental_heatmap,
  width = 8.5,
  height = 4.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap.pdf",
  plot = environmental_heatmap,
  width = 8.5,
  height = 4.8,
  units = "in",
  bg = "white"
)

cat("\nSaved heatmap:\n")
cat("environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap.png\n")
cat("environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap.pdf\n")
cat("environmental_association_outputs/Sx_environmental_marginal_dbRDA_heatmap_data.csv\n")

# ============================================================
# R3.20 — ENVIRONMENT × ECOLOGICAL CONDITION
# STEP 1 — DESIGN / ESTIMABILITY CHECK
#
# Question:
# Can ecological condition and the six measured environmental
# variables be included simultaneously after accounting for Site?
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

stopifnot(
  exists("env_results")
)

analysis_names <- names(env_results)

for (nm in analysis_names) {
  
  cat("\n\n========================================\n")
  cat("ANALYSIS:", nm, "\n")
  cat("========================================\n")
  
  dat <- env_results[[nm]]$metadata_used %>%
    as.data.frame() %>%
    mutate(
      Site = factor(Site),
      Treatment = factor(Treatment)
    )
  
  # ----------------------------------------------------------
  # 1. SAMPLE SIZE
  # ----------------------------------------------------------
  
  cat(
    "Samples:",
    nrow(dat),
    "\n\n"
  )
  
  # ----------------------------------------------------------
  # 2. SITE × ECOLOGICAL CONDITION COUNTS
  # ----------------------------------------------------------
  
  cat("Site × ecological condition:\n")
  
  print(
    with(
      dat,
      table(
        Site,
        Treatment
      )
    )
  )
  
  # ----------------------------------------------------------
  # 3. CONDITION TOTALS
  # ----------------------------------------------------------
  
  cat("\nEcological-condition totals:\n")
  
  print(
    dat %>%
      count(
        Treatment,
        name = "n"
      )
  )
  
  # ----------------------------------------------------------
  # 4. FULL DESIGN MATRIX
  #
  # Site + ecological condition + all six environmental vars
  # ----------------------------------------------------------
  
  X <- model.matrix(
    ~ Site +
      Treatment +
      Moisture_z +
      Temperature_z +
      pH_z +
      NH4_log_z +
      PO4_log_z +
      Redox_z,
    data = dat
  )
  
  matrix_rank <- qr(X)$rank
  matrix_columns <- ncol(X)
  
  cat("\nDesign matrix columns:",
      matrix_columns,
      "\n")
  
  cat(
    "Design matrix rank:",
    matrix_rank,
    "\n"
  )
  
  cat(
    "Full rank:",
    matrix_rank == matrix_columns,
    "\n"
  )
  
  # ----------------------------------------------------------
  # 5. CONDITION NUMBER
  #
  # Larger values indicate increasingly severe collinearity.
  # This is diagnostic only.
  # ----------------------------------------------------------
  
  kappa_value <- kappa(
    X,
    exact = TRUE
  )
  
  cat(
    "Condition number:",
    round(kappa_value, 2),
    "\n"
  )
}

# ============================================================
# R3.20 — DESIGN / COLLINEARITY PREFLIGHT
# Does the dataset support modelling ecological condition
# and environmental variables simultaneously?
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
})

# ------------------------------------------------------------
# 1. Helper to find metadata columns
# ------------------------------------------------------------

find_col_preflight <- function(df, candidates) {
  
  exact <- candidates[
    candidates %in% names(df)
  ]
  
  if (length(exact) > 0) {
    return(exact[1])
  }
  
  lower_names <- tolower(names(df))
  lower_candidates <- tolower(candidates)
  
  hit <- match(
    lower_candidates,
    lower_names
  )
  
  hit <- hit[!is.na(hit)]
  
  if (length(hit) > 0) {
    return(
      names(df)[hit[1]]
    )
  }
  
  stop(
    "Could not find one of: ",
    paste(candidates, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 2. Prepare one habitat × kingdom dataset
# ------------------------------------------------------------

prepare_preflight <- function(
    ps,
    habitat,
    analysis_name
) {
  
  dat <- data.frame(
    sample_data(ps),
    check.names = FALSE
  )
  
  site_col <- find_col_preflight(
    dat,
    c("Site", "SITE", "site")
  )
  
  treatment_col <- find_col_preflight(
    dat,
    c(
      "Treatment",
      "TREATMENT",
      "EcologicalStatus",
      "Ecological_condition"
    )
  )
  
  habitat_col <- find_col_preflight(
    dat,
    c(
      "SoilPortion",
      "Soil_habitat",
      "SoilHabitat",
      "Compartment",
      "Habitat"
    )
  )
  
  moisture_col <- find_col_preflight(
    dat,
    c(
      "Moisture_avg",
      "Moisture",
      "Soil_moisture"
    )
  )
  
  temp_col <- find_col_preflight(
    dat,
    c(
      "Temperature_avg",
      "Temperature",
      "Temp"
    )
  )
  
  ph_col <- find_col_preflight(
    dat,
    c(
      "pH_avg",
      "pH",
      "PH"
    )
  )
  
  nh4_col <- find_col_preflight(
    dat,
    c(
      "NH4",
      "Ammonium",
      "NH4_avg"
    )
  )
  
  po4_col <- find_col_preflight(
    dat,
    c(
      "PO4",
      "Phosphate",
      "PO4_avg"
    )
  )
  
  redox_col <- find_col_preflight(
    dat,
    c(
      "Redox_avg",
      "Redox",
      "Eh"
    )
  )
  
  # Standardise habitat
  soil_habitat <- toupper(
    trimws(
      as.character(
        dat[[habitat_col]]
      )
    )
  )
  
  soil_habitat[
    soil_habitat %in%
      c(
        "BULK",
        "BULK SOIL",
        "BULK_SOIL"
      )
  ] <- "BS"
  
  soil_habitat[
    soil_habitat == "RHIZOSPHERE"
  ] <- "RZ"
  
  # Standardise condition
  condition <- toupper(
    trimws(
      as.character(
        dat[[treatment_col]]
      )
    )
  )
  
  condition[
    condition %in% c("NAT", "NATURAL")
  ] <- "Natural"
  
  condition[
    condition == "BESE"
  ] <- "BESE"
  
  condition[
    condition == "BARE"
  ] <- "BARE"
  
  condition[
    condition %in%
      c("DEG", "DEGRADED", "DEGREST")
  ] <- "Degraded"
  
  out <- data.frame(
    Site = as.character(
      dat[[site_col]]
    ),
    Treatment = condition,
    Soil_habitat = soil_habitat,
    Moisture = as.numeric(
      dat[[moisture_col]]
    ),
    Temperature = as.numeric(
      dat[[temp_col]]
    ),
    pH = as.numeric(
      dat[[ph_col]]
    ),
    NH4 = as.numeric(
      dat[[nh4_col]]
    ),
    PO4 = as.numeric(
      dat[[po4_col]]
    ),
    Redox = as.numeric(
      dat[[redox_col]]
    )
  )
  
  out <- out %>%
    filter(
      Soil_habitat == habitat
    ) %>%
    tidyr::drop_na(
      Moisture,
      Temperature,
      pH,
      NH4,
      PO4,
      Redox
    )
  
  # Match transformations used previously
  out$NH4_log <- if (
    any(out$NH4 <= 0)
  ) {
    log1p(out$NH4)
  } else {
    log(out$NH4)
  }
  
  out$PO4_log <- if (
    any(out$PO4 <= 0)
  ) {
    log1p(out$PO4)
  } else {
    log(out$PO4)
  }
  
  out <- out %>%
    mutate(
      Site = factor(Site),
      Treatment = factor(Treatment),
      Moisture_z =
        as.numeric(scale(Moisture)),
      Temperature_z =
        as.numeric(scale(Temperature)),
      pH_z =
        as.numeric(scale(pH)),
      NH4_log_z =
        as.numeric(scale(NH4_log)),
      PO4_log_z =
        as.numeric(scale(PO4_log)),
      Redox_z =
        as.numeric(scale(Redox))
    )
  
  cat("\n\n========================================\n")
  cat("ANALYSIS:", analysis_name, "\n")
  cat("========================================\n")
  
  cat(
    "Samples:",
    nrow(out),
    "\n\n"
  )
  
  cat("Site × ecological condition:\n")
  
  print(
    with(
      out,
      table(
        Site,
        Treatment
      )
    )
  )
  
  cat("\nCondition totals:\n")
  
  print(
    out %>%
      count(
        Treatment,
        name = "n"
      )
  )
  
  # Full combined design
  X <- model.matrix(
    ~ Site +
      Treatment +
      Moisture_z +
      Temperature_z +
      pH_z +
      NH4_log_z +
      PO4_log_z +
      Redox_z,
    data = out
  )
  
  cat(
    "\nDesign matrix columns:",
    ncol(X),
    "\n"
  )
  
  cat(
    "Design matrix rank:",
    qr(X)$rank,
    "\n"
  )
  
  cat(
    "Full rank:",
    qr(X)$rank == ncol(X),
    "\n"
  )
  
  cat(
    "Condition number:",
    round(
      kappa(
        X,
        exact = TRUE
      ),
      2
    ),
    "\n"
  )
  
  invisible(out)
}

# ------------------------------------------------------------
# 3. Run all four analyses
# ------------------------------------------------------------

preflight_bulk_prok <- prepare_preflight(
  ps_prok,
  "BS",
  "Bulk soil prokaryotic"
)

preflight_bulk_fun <- prepare_preflight(
  ps_fun,
  "BS",
  "Bulk soil fungal"
)

preflight_rz_prok <- prepare_preflight(
  ps_prok,
  "RZ",
  "Rhizosphere prokaryotic"
)

preflight_rz_fun <- prepare_preflight(
  ps_fun,
  "RZ",
  "Rhizosphere fungal"
)
# ============================================================
# R3.20 — RECIPROCAL PARTIAL db-RDA
# STEP 2A — BULK-SOIL PROKARYOTES ONLY
#
# H1:
# Environment | Site + Ecological condition
#
# H2:
# Ecological condition | Site + Environment
#
# These estimate statistically unique associations.
# They do NOT establish causality.
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(vegan)
})

set.seed(123)

stopifnot(
  exists("ps_prok")
)

# ============================================================
# 1. METADATA
# ============================================================

meta_joint <- data.frame(
  sample_data(ps_prok),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_joint$SampleID <- rownames(meta_joint)

# Find relevant columns
site_col <- find_col_preflight(
  meta_joint,
  c("Site", "SITE", "site")
)

treatment_col <- find_col_preflight(
  meta_joint,
  c(
    "Treatment",
    "TREATMENT",
    "EcologicalStatus",
    "Ecological_condition"
  )
)

habitat_col <- find_col_preflight(
  meta_joint,
  c(
    "SoilPortion",
    "Soil_habitat",
    "SoilHabitat",
    "Compartment",
    "Habitat"
  )
)

moisture_col <- find_col_preflight(
  meta_joint,
  c(
    "Moisture_avg",
    "Moisture",
    "Soil_moisture"
  )
)

temp_col <- find_col_preflight(
  meta_joint,
  c(
    "Temperature_avg",
    "Temperature",
    "Temp"
  )
)

ph_col <- find_col_preflight(
  meta_joint,
  c(
    "pH_avg",
    "pH",
    "PH"
  )
)

nh4_col <- find_col_preflight(
  meta_joint,
  c(
    "NH4",
    "Ammonium",
    "NH4_avg"
  )
)

po4_col <- find_col_preflight(
  meta_joint,
  c(
    "PO4",
    "Phosphate",
    "PO4_avg"
  )
)

redox_col <- find_col_preflight(
  meta_joint,
  c(
    "Redox_avg",
    "Redox",
    "Eh"
  )
)

# ============================================================
# 2. CLEAN BULK-SOIL METADATA
# ============================================================

dat_joint <- meta_joint %>%
  transmute(
    SampleID = SampleID,
    
    Site = as.character(
      .data[[site_col]]
    ),
    
    Treatment = toupper(
      trimws(
        as.character(
          .data[[treatment_col]]
        )
      )
    ),
    
    Soil_habitat = toupper(
      trimws(
        as.character(
          .data[[habitat_col]]
        )
      )
    ),
    
    Moisture = as.numeric(
      .data[[moisture_col]]
    ),
    
    Temperature = as.numeric(
      .data[[temp_col]]
    ),
    
    pH = as.numeric(
      .data[[ph_col]]
    ),
    
    NH4 = as.numeric(
      .data[[nh4_col]]
    ),
    
    PO4 = as.numeric(
      .data[[po4_col]]
    ),
    
    Redox = as.numeric(
      .data[[redox_col]]
    )
  )

# Standardise labels
dat_joint$Treatment[
  dat_joint$Treatment %in%
    c("NAT", "NATURAL")
] <- "Natural"

dat_joint$Treatment[
  dat_joint$Treatment == "BESE"
] <- "BESE"

dat_joint$Treatment[
  dat_joint$Treatment == "BARE"
] <- "BARE"

dat_joint$Treatment[
  dat_joint$Treatment %in%
    c(
      "DEG",
      "DEGRADED",
      "DEGREST"
    )
] <- "Degraded"

dat_joint$Soil_habitat[
  dat_joint$Soil_habitat %in%
    c(
      "BULK",
      "BULK SOIL",
      "BULK_SOIL"
    )
] <- "BS"

dat_joint <- dat_joint %>%
  filter(
    Soil_habitat == "BS"
  ) %>%
  drop_na(
    Moisture,
    Temperature,
    pH,
    NH4,
    PO4,
    Redox
  )

# ============================================================
# 3. ENVIRONMENTAL TRANSFORMATIONS
# Identical logic to previous db-RDA
# ============================================================

if (any(dat_joint$NH4 <= 0)) {
  dat_joint$NH4_log <- log1p(
    dat_joint$NH4
  )
} else {
  dat_joint$NH4_log <- log(
    dat_joint$NH4
  )
}

if (any(dat_joint$PO4 <= 0)) {
  dat_joint$PO4_log <- log1p(
    dat_joint$PO4
  )
} else {
  dat_joint$PO4_log <- log(
    dat_joint$PO4
  )
}

dat_joint <- dat_joint %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Natural",
        "BESE",
        "BARE",
        "Degraded"
      )
    ),
    
    Moisture_z =
      as.numeric(scale(Moisture)),
    
    Temperature_z =
      as.numeric(scale(Temperature)),
    
    pH_z =
      as.numeric(scale(pH)),
    
    NH4_log_z =
      as.numeric(scale(NH4_log)),
    
    PO4_log_z =
      as.numeric(scale(PO4_log)),
    
    Redox_z =
      as.numeric(scale(Redox))
  )

# ============================================================
# 4. COMMUNITY MATRIX
# ============================================================

ps_joint <- prune_samples(
  dat_joint$SampleID,
  ps_prok
)

ps_joint <- prune_taxa(
  taxa_sums(ps_joint) > 0,
  ps_joint
)

ps_joint_rel <- transform_sample_counts(
  ps_joint,
  function(x) x / sum(x)
)

comm_joint <- as(
  otu_table(ps_joint_rel),
  "matrix"
)

if (
  taxa_are_rows(ps_joint_rel)
) {
  comm_joint <- t(
    comm_joint
  )
}

# Exact sample order
comm_joint <- comm_joint[
  dat_joint$SampleID,
  ,
  drop = FALSE
]

comm_joint <- comm_joint[
  ,
  colSums(comm_joint) > 0,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(comm_joint),
    dat_joint$SampleID
  )
)

cat("\n========================================\n")
cat("BULK PROKARYOTIC — JOINT db-RDA\n")
cat("========================================\n")

cat(
  "Samples:",
  nrow(dat_joint),
  "\n"
)

cat(
  "ASVs:",
  ncol(comm_joint),
  "\n\n"
)

print(
  table(
    dat_joint$Site,
    dat_joint$Treatment
  )
)

# ============================================================
# 5. H1
# UNIQUE ENVIRONMENTAL ASSOCIATION
# after Site + ecological condition
# ============================================================

db_env_unique <- vegan::dbrda(
  comm_joint ~
    Moisture_z +
    Temperature_z +
    pH_z +
    NH4_log_z +
    PO4_log_z +
    Redox_z +
    Condition(
      Site +
        Treatment
    ),
  data = dat_joint,
  distance = "bray",
  add = "lingoes"
)

env_unique_overall <- vegan::anova.cca(
  db_env_unique,
  permutations = 9999,
  strata = dat_joint$Site
)

env_unique_margin <- vegan::anova.cca(
  db_env_unique,
  by = "margin",
  permutations = 9999,
  strata = dat_joint$Site
)

env_unique_r2 <- vegan::RsquareAdj(
  db_env_unique
)

cat("\n========================================\n")
cat("A. ENVIRONMENT | SITE + CONDITION\n")
cat("========================================\n")

print(
  env_unique_overall
)

cat("\nMarginal environmental terms:\n")

print(
  env_unique_margin
)

cat("\nR2:\n")

print(
  env_unique_r2
)

cat("\nVIF:\n")

print(
  vegan::vif.cca(
    db_env_unique
  )
)

# ============================================================
# 6. H2
# UNIQUE ECOLOGICAL-CONDITION ASSOCIATION
# after Site + measured environment
# ============================================================

db_condition_unique <- vegan::dbrda(
  comm_joint ~
    Treatment +
    Condition(
      Site +
        Moisture_z +
        Temperature_z +
        pH_z +
        NH4_log_z +
        PO4_log_z +
        Redox_z
    ),
  data = dat_joint,
  distance = "bray",
  add = "lingoes"
)

condition_unique_test <- vegan::anova.cca(
  db_condition_unique,
  permutations = 9999,
  strata = dat_joint$Site
)

condition_unique_r2 <- vegan::RsquareAdj(
  db_condition_unique
)

cat("\n========================================\n")
cat("B. CONDITION | SITE + ENVIRONMENT\n")
cat("========================================\n")

print(
  condition_unique_test
)

cat("\nR2:\n")

print(
  condition_unique_r2
)

# ============================================================
# R3.20 — RECIPROCAL PARTIAL db-RDA
# STEP 2B — BULK-SOIL FUNGI
#
# H1:
# Environment | Site + Ecological condition
#
# H2:
# Ecological condition | Site + Environment
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(vegan)
})

set.seed(123)

stopifnot(
  exists("ps_fun"),
  exists("find_col_preflight")
)

# ============================================================
# 1. METADATA
# ============================================================

meta_joint_fun <- data.frame(
  sample_data(ps_fun),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_joint_fun$SampleID <- rownames(meta_joint_fun)

site_col_fun <- find_col_preflight(
  meta_joint_fun,
  c("Site", "SITE", "site")
)

treatment_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "Treatment",
    "TREATMENT",
    "EcologicalStatus",
    "Ecological_condition"
  )
)

habitat_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "SoilPortion",
    "Soil_habitat",
    "SoilHabitat",
    "Compartment",
    "Habitat"
  )
)

moisture_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "Moisture_avg",
    "Moisture",
    "Soil_moisture"
  )
)

temp_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "Temperature_avg",
    "Temperature",
    "Temp"
  )
)

ph_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "pH_avg",
    "pH",
    "PH"
  )
)

nh4_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "NH4",
    "Ammonium",
    "NH4_avg"
  )
)

po4_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "PO4",
    "Phosphate",
    "PO4_avg"
  )
)

redox_col_fun <- find_col_preflight(
  meta_joint_fun,
  c(
    "Redox_avg",
    "Redox",
    "Eh"
  )
)

# ============================================================
# 2. CLEAN BULK-SOIL METADATA
# ============================================================

dat_joint_fun <- meta_joint_fun %>%
  transmute(
    SampleID = SampleID,
    
    Site = as.character(
      .data[[site_col_fun]]
    ),
    
    Treatment = toupper(
      trimws(
        as.character(
          .data[[treatment_col_fun]]
        )
      )
    ),
    
    Soil_habitat = toupper(
      trimws(
        as.character(
          .data[[habitat_col_fun]]
        )
      )
    ),
    
    Moisture = as.numeric(
      .data[[moisture_col_fun]]
    ),
    
    Temperature = as.numeric(
      .data[[temp_col_fun]]
    ),
    
    pH = as.numeric(
      .data[[ph_col_fun]]
    ),
    
    NH4 = as.numeric(
      .data[[nh4_col_fun]]
    ),
    
    PO4 = as.numeric(
      .data[[po4_col_fun]]
    ),
    
    Redox = as.numeric(
      .data[[redox_col_fun]]
    )
  )

# Standardise treatment labels
dat_joint_fun$Treatment[
  dat_joint_fun$Treatment %in%
    c("NAT", "NATURAL")
] <- "Natural"

dat_joint_fun$Treatment[
  dat_joint_fun$Treatment == "BESE"
] <- "BESE"

dat_joint_fun$Treatment[
  dat_joint_fun$Treatment == "BARE"
] <- "BARE"

dat_joint_fun$Treatment[
  dat_joint_fun$Treatment %in%
    c("DEG", "DEGRADED", "DEGREST")
] <- "Degraded"

# Standardise habitat
dat_joint_fun$Soil_habitat[
  dat_joint_fun$Soil_habitat %in%
    c(
      "BULK",
      "BULK SOIL",
      "BULK_SOIL"
    )
] <- "BS"

dat_joint_fun <- dat_joint_fun %>%
  filter(
    Soil_habitat == "BS"
  ) %>%
  drop_na(
    Moisture,
    Temperature,
    pH,
    NH4,
    PO4,
    Redox
  )

# ============================================================
# 3. ENVIRONMENTAL TRANSFORMATIONS
# Same logic as previous environmental analysis
# ============================================================

if (any(dat_joint_fun$NH4 <= 0)) {
  dat_joint_fun$NH4_log <- log1p(
    dat_joint_fun$NH4
  )
} else {
  dat_joint_fun$NH4_log <- log(
    dat_joint_fun$NH4
  )
}

if (any(dat_joint_fun$PO4 <= 0)) {
  dat_joint_fun$PO4_log <- log1p(
    dat_joint_fun$PO4
  )
} else {
  dat_joint_fun$PO4_log <- log(
    dat_joint_fun$PO4
  )
}

dat_joint_fun <- dat_joint_fun %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Natural",
        "BESE",
        "BARE",
        "Degraded"
      )
    ),
    
    Moisture_z =
      as.numeric(scale(Moisture)),
    
    Temperature_z =
      as.numeric(scale(Temperature)),
    
    pH_z =
      as.numeric(scale(pH)),
    
    NH4_log_z =
      as.numeric(scale(NH4_log)),
    
    PO4_log_z =
      as.numeric(scale(PO4_log)),
    
    Redox_z =
      as.numeric(scale(Redox))
  )

# ============================================================
# 4. COMMUNITY MATRIX
# ============================================================

ps_joint_fun <- prune_samples(
  dat_joint_fun$SampleID,
  ps_fun
)

ps_joint_fun <- prune_taxa(
  taxa_sums(ps_joint_fun) > 0,
  ps_joint_fun
)

ps_joint_rel_fun <- transform_sample_counts(
  ps_joint_fun,
  function(x) x / sum(x)
)

comm_joint_fun <- as(
  otu_table(ps_joint_rel_fun),
  "matrix"
)

if (taxa_are_rows(ps_joint_rel_fun)) {
  comm_joint_fun <- t(
    comm_joint_fun
  )
}

comm_joint_fun <- comm_joint_fun[
  dat_joint_fun$SampleID,
  ,
  drop = FALSE
]

comm_joint_fun <- comm_joint_fun[
  ,
  colSums(comm_joint_fun) > 0,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(comm_joint_fun),
    dat_joint_fun$SampleID
  )
)

cat("\n========================================\n")
cat("BULK FUNGAL — JOINT db-RDA\n")
cat("========================================\n")

cat(
  "Samples:",
  nrow(dat_joint_fun),
  "\n"
)

cat(
  "ASVs:",
  ncol(comm_joint_fun),
  "\n\n"
)

print(
  table(
    dat_joint_fun$Site,
    dat_joint_fun$Treatment
  )
)

# ============================================================
# 5. H1
# ENVIRONMENT | SITE + ECOLOGICAL CONDITION
# ============================================================

db_env_unique_fun <- vegan::dbrda(
  comm_joint_fun ~
    Moisture_z +
    Temperature_z +
    pH_z +
    NH4_log_z +
    PO4_log_z +
    Redox_z +
    Condition(
      Site +
        Treatment
    ),
  data = dat_joint_fun,
  distance = "bray",
  add = "lingoes"
)

env_unique_overall_fun <- vegan::anova.cca(
  db_env_unique_fun,
  permutations = 9999,
  strata = dat_joint_fun$Site
)

env_unique_margin_fun <- vegan::anova.cca(
  db_env_unique_fun,
  by = "margin",
  permutations = 9999,
  strata = dat_joint_fun$Site
)

env_unique_r2_fun <- vegan::RsquareAdj(
  db_env_unique_fun
)

cat("\n========================================\n")
cat("A. ENVIRONMENT | SITE + CONDITION\n")
cat("========================================\n")

print(
  env_unique_overall_fun
)

cat("\nMarginal environmental terms:\n")

print(
  env_unique_margin_fun
)

cat("\nR2:\n")

print(
  env_unique_r2_fun
)

cat("\nVIF:\n")

print(
  vegan::vif.cca(
    db_env_unique_fun
  )
)

# ============================================================
# 6. BH CORRECTION FOR SIX MARGINAL ENVIRONMENTAL TESTS
# ============================================================

env_margin_fun_table <- as.data.frame(
  env_unique_margin_fun
)

env_margin_fun_table$Variable <-
  rownames(env_margin_fun_table)

rownames(env_margin_fun_table) <- NULL

env_margin_fun_table <- env_margin_fun_table %>%
  filter(
    Variable != "Residual"
  ) %>%
  mutate(
    p_BH = p.adjust(
      `Pr(>F)`,
      method = "BH"
    ),
    Significant_BH =
      p_BH < 0.05
  )

cat("\n========================================\n")
cat("BH-ADJUSTED MARGINAL ENVIRONMENTAL TESTS\n")
cat("========================================\n")

print(
  env_margin_fun_table,
  row.names = FALSE
)

# ============================================================
# 7. H2
# ECOLOGICAL CONDITION | SITE + ENVIRONMENT
# ============================================================

db_condition_unique_fun <- vegan::dbrda(
  comm_joint_fun ~
    Treatment +
    Condition(
      Site +
        Moisture_z +
        Temperature_z +
        pH_z +
        NH4_log_z +
        PO4_log_z +
        Redox_z
    ),
  data = dat_joint_fun,
  distance = "bray",
  add = "lingoes"
)

condition_unique_test_fun <- vegan::anova.cca(
  db_condition_unique_fun,
  permutations = 9999,
  strata = dat_joint_fun$Site
)

condition_unique_r2_fun <- vegan::RsquareAdj(
  db_condition_unique_fun
)

cat("\n========================================\n")
cat("B. CONDITION | SITE + ENVIRONMENT\n")
cat("========================================\n")

print(
  condition_unique_test_fun
)

cat("\nR2:\n")

print(
  condition_unique_r2_fun
)


# ============================================================
# R3.20 — RECIPROCAL PARTIAL db-RDA
# STEP 2C — RHIZOSPHERE PROKARYOTES
#
# H1:
# Environment | Site + Ecological condition
#
# H2:
# Ecological condition | Site + Environment
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(vegan)
})

set.seed(123)

stopifnot(
  exists("ps_prok"),
  exists("find_col_preflight")
)

# ============================================================
# 1. METADATA
# ============================================================

meta_joint_rz_prok <- data.frame(
  sample_data(ps_prok),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_joint_rz_prok$SampleID <-
  rownames(meta_joint_rz_prok)

site_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c("Site", "SITE", "site")
)

treatment_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "Treatment",
    "TREATMENT",
    "EcologicalStatus",
    "Ecological_condition"
  )
)

habitat_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "SoilPortion",
    "Soil_habitat",
    "SoilHabitat",
    "Compartment",
    "Habitat"
  )
)

moisture_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "Moisture_avg",
    "Moisture",
    "Soil_moisture"
  )
)

temp_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "Temperature_avg",
    "Temperature",
    "Temp"
  )
)

ph_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "pH_avg",
    "pH",
    "PH"
  )
)

nh4_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "NH4",
    "Ammonium",
    "NH4_avg"
  )
)

po4_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "PO4",
    "Phosphate",
    "PO4_avg"
  )
)

redox_col_rz_prok <- find_col_preflight(
  meta_joint_rz_prok,
  c(
    "Redox_avg",
    "Redox",
    "Eh"
  )
)

# ============================================================
# 2. CLEAN RHIZOSPHERE METADATA
# ============================================================

dat_joint_rz_prok <- meta_joint_rz_prok %>%
  transmute(
    SampleID = SampleID,
    
    Site = as.character(
      .data[[site_col_rz_prok]]
    ),
    
    Treatment = toupper(
      trimws(
        as.character(
          .data[[treatment_col_rz_prok]]
        )
      )
    ),
    
    Soil_habitat = toupper(
      trimws(
        as.character(
          .data[[habitat_col_rz_prok]]
        )
      )
    ),
    
    Moisture = as.numeric(
      .data[[moisture_col_rz_prok]]
    ),
    
    Temperature = as.numeric(
      .data[[temp_col_rz_prok]]
    ),
    
    pH = as.numeric(
      .data[[ph_col_rz_prok]]
    ),
    
    NH4 = as.numeric(
      .data[[nh4_col_rz_prok]]
    ),
    
    PO4 = as.numeric(
      .data[[po4_col_rz_prok]]
    ),
    
    Redox = as.numeric(
      .data[[redox_col_rz_prok]]
    )
  )

# Standardise treatment
dat_joint_rz_prok$Treatment[
  dat_joint_rz_prok$Treatment %in%
    c("NAT", "NATURAL")
] <- "Natural"

dat_joint_rz_prok$Treatment[
  dat_joint_rz_prok$Treatment == "BESE"
] <- "BESE"

dat_joint_rz_prok$Treatment[
  dat_joint_rz_prok$Treatment == "BARE"
] <- "BARE"

# Standardise habitat
dat_joint_rz_prok$Soil_habitat[
  dat_joint_rz_prok$Soil_habitat == "RHIZOSPHERE"
] <- "RZ"

dat_joint_rz_prok <- dat_joint_rz_prok %>%
  filter(
    Soil_habitat == "RZ",
    Treatment %in% c(
      "Natural",
      "BESE",
      "BARE"
    )
  ) %>%
  drop_na(
    Moisture,
    Temperature,
    pH,
    NH4,
    PO4,
    Redox
  )

# ============================================================
# 3. ENVIRONMENTAL TRANSFORMATIONS
# ============================================================

if (any(dat_joint_rz_prok$NH4 <= 0)) {
  dat_joint_rz_prok$NH4_log <-
    log1p(dat_joint_rz_prok$NH4)
} else {
  dat_joint_rz_prok$NH4_log <-
    log(dat_joint_rz_prok$NH4)
}

if (any(dat_joint_rz_prok$PO4 <= 0)) {
  dat_joint_rz_prok$PO4_log <-
    log1p(dat_joint_rz_prok$PO4)
} else {
  dat_joint_rz_prok$PO4_log <-
    log(dat_joint_rz_prok$PO4)
}

dat_joint_rz_prok <- dat_joint_rz_prok %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    ),
    
    Moisture_z =
      as.numeric(scale(Moisture)),
    
    Temperature_z =
      as.numeric(scale(Temperature)),
    
    pH_z =
      as.numeric(scale(pH)),
    
    NH4_log_z =
      as.numeric(scale(NH4_log)),
    
    PO4_log_z =
      as.numeric(scale(PO4_log)),
    
    Redox_z =
      as.numeric(scale(Redox))
  )

# ============================================================
# 4. COMMUNITY MATRIX
# ============================================================

ps_joint_rz_prok <- prune_samples(
  dat_joint_rz_prok$SampleID,
  ps_prok
)

ps_joint_rz_prok <- prune_taxa(
  taxa_sums(ps_joint_rz_prok) > 0,
  ps_joint_rz_prok
)

ps_joint_rz_prok_rel <-
  transform_sample_counts(
    ps_joint_rz_prok,
    function(x) x / sum(x)
  )

comm_joint_rz_prok <- as(
  otu_table(ps_joint_rz_prok_rel),
  "matrix"
)

if (taxa_are_rows(ps_joint_rz_prok_rel)) {
  comm_joint_rz_prok <-
    t(comm_joint_rz_prok)
}

comm_joint_rz_prok <-
  comm_joint_rz_prok[
    dat_joint_rz_prok$SampleID,
    ,
    drop = FALSE
  ]

comm_joint_rz_prok <-
  comm_joint_rz_prok[
    ,
    colSums(comm_joint_rz_prok) > 0,
    drop = FALSE
  ]

stopifnot(
  identical(
    rownames(comm_joint_rz_prok),
    dat_joint_rz_prok$SampleID
  )
)

cat("\n========================================\n")
cat("RHIZOSPHERE PROKARYOTIC — JOINT db-RDA\n")
cat("========================================\n")

cat(
  "Samples:",
  nrow(dat_joint_rz_prok),
  "\n"
)

cat(
  "ASVs:",
  ncol(comm_joint_rz_prok),
  "\n\n"
)

print(
  table(
    dat_joint_rz_prok$Site,
    dat_joint_rz_prok$Treatment
  )
)

# ============================================================
# 5. H1
# ENVIRONMENT | SITE + ECOLOGICAL CONDITION
# ============================================================

db_env_unique_rz_prok <- vegan::dbrda(
  comm_joint_rz_prok ~
    Moisture_z +
    Temperature_z +
    pH_z +
    NH4_log_z +
    PO4_log_z +
    Redox_z +
    Condition(
      Site +
        Treatment
    ),
  data = dat_joint_rz_prok,
  distance = "bray",
  add = "lingoes"
)

env_unique_overall_rz_prok <-
  vegan::anova.cca(
    db_env_unique_rz_prok,
    permutations = 9999,
    strata = dat_joint_rz_prok$Site
  )

env_unique_margin_rz_prok <-
  vegan::anova.cca(
    db_env_unique_rz_prok,
    by = "margin",
    permutations = 9999,
    strata = dat_joint_rz_prok$Site
  )

env_unique_r2_rz_prok <-
  vegan::RsquareAdj(
    db_env_unique_rz_prok
  )

cat("\n========================================\n")
cat("A. ENVIRONMENT | SITE + CONDITION\n")
cat("========================================\n")

print(
  env_unique_overall_rz_prok
)

cat("\nMarginal environmental terms:\n")

print(
  env_unique_margin_rz_prok
)

cat("\nR2:\n")

print(
  env_unique_r2_rz_prok
)

cat("\nVIF:\n")

print(
  vegan::vif.cca(
    db_env_unique_rz_prok
  )
)

# ============================================================
# 6. BH CORRECTION
# ============================================================

env_margin_rz_prok_table <-
  as.data.frame(
    env_unique_margin_rz_prok
  )

env_margin_rz_prok_table$Variable <-
  rownames(
    env_margin_rz_prok_table
  )

rownames(
  env_margin_rz_prok_table
) <- NULL

env_margin_rz_prok_table <-
  env_margin_rz_prok_table %>%
  filter(
    Variable != "Residual"
  ) %>%
  mutate(
    p_BH = p.adjust(
      `Pr(>F)`,
      method = "BH"
    ),
    
    Significant_BH =
      p_BH < 0.05
  )

cat("\n========================================\n")
cat("BH-ADJUSTED MARGINAL ENVIRONMENTAL TESTS\n")
cat("========================================\n")

print(
  env_margin_rz_prok_table,
  row.names = FALSE
)

# ============================================================
# 7. H2
# ECOLOGICAL CONDITION | SITE + ENVIRONMENT
# ============================================================

db_condition_unique_rz_prok <- vegan::dbrda(
  comm_joint_rz_prok ~
    Treatment +
    Condition(
      Site +
        Moisture_z +
        Temperature_z +
        pH_z +
        NH4_log_z +
        PO4_log_z +
        Redox_z
    ),
  data = dat_joint_rz_prok,
  distance = "bray",
  add = "lingoes"
)

condition_unique_test_rz_prok <-
  vegan::anova.cca(
    db_condition_unique_rz_prok,
    permutations = 9999,
    strata = dat_joint_rz_prok$Site
  )

condition_unique_r2_rz_prok <-
  vegan::RsquareAdj(
    db_condition_unique_rz_prok
  )

cat("\n========================================\n")
cat("B. CONDITION | SITE + ENVIRONMENT\n")
cat("========================================\n")

print(
  condition_unique_test_rz_prok
)

cat("\nR2:\n")

print(
  condition_unique_r2_rz_prok
)


# ============================================================
# R3.20 — RECIPROCAL PARTIAL db-RDA
# STEP 2D — RHIZOSPHERE FUNGI
#
# H1:
# Environment | Site + Ecological condition
#
# H2:
# Ecological condition | Site + Environment
# ============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(vegan)
})

set.seed(123)

stopifnot(
  exists("ps_fun"),
  exists("find_col_preflight")
)

# ============================================================
# 1. METADATA
# ============================================================

meta_joint_rz_fun <- data.frame(
  sample_data(ps_fun),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

meta_joint_rz_fun$SampleID <-
  rownames(meta_joint_rz_fun)

site_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c("Site", "SITE", "site")
)

treatment_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "Treatment",
    "TREATMENT",
    "EcologicalStatus",
    "Ecological_condition"
  )
)

habitat_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "SoilPortion",
    "Soil_habitat",
    "SoilHabitat",
    "Compartment",
    "Habitat"
  )
)

moisture_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "Moisture_avg",
    "Moisture",
    "Soil_moisture"
  )
)

temp_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "Temperature_avg",
    "Temperature",
    "Temp"
  )
)

ph_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "pH_avg",
    "pH",
    "PH"
  )
)

nh4_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "NH4",
    "Ammonium",
    "NH4_avg"
  )
)

po4_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "PO4",
    "Phosphate",
    "PO4_avg"
  )
)

redox_col_rz_fun <- find_col_preflight(
  meta_joint_rz_fun,
  c(
    "Redox_avg",
    "Redox",
    "Eh"
  )
)

# ============================================================
# 2. CLEAN RHIZOSPHERE METADATA
# ============================================================

dat_joint_rz_fun <- meta_joint_rz_fun %>%
  transmute(
    SampleID = SampleID,
    
    Site = as.character(
      .data[[site_col_rz_fun]]
    ),
    
    Treatment = toupper(
      trimws(
        as.character(
          .data[[treatment_col_rz_fun]]
        )
      )
    ),
    
    Soil_habitat = toupper(
      trimws(
        as.character(
          .data[[habitat_col_rz_fun]]
        )
      )
    ),
    
    Moisture = as.numeric(
      .data[[moisture_col_rz_fun]]
    ),
    
    Temperature = as.numeric(
      .data[[temp_col_rz_fun]]
    ),
    
    pH = as.numeric(
      .data[[ph_col_rz_fun]]
    ),
    
    NH4 = as.numeric(
      .data[[nh4_col_rz_fun]]
    ),
    
    PO4 = as.numeric(
      .data[[po4_col_rz_fun]]
    ),
    
    Redox = as.numeric(
      .data[[redox_col_rz_fun]]
    )
  )

# Standardise treatment
dat_joint_rz_fun$Treatment[
  dat_joint_rz_fun$Treatment %in%
    c("NAT", "NATURAL")
] <- "Natural"

dat_joint_rz_fun$Treatment[
  dat_joint_rz_fun$Treatment == "BESE"
] <- "BESE"

dat_joint_rz_fun$Treatment[
  dat_joint_rz_fun$Treatment == "BARE"
] <- "BARE"

# Standardise habitat
dat_joint_rz_fun$Soil_habitat[
  dat_joint_rz_fun$Soil_habitat == "RHIZOSPHERE"
] <- "RZ"

dat_joint_rz_fun <- dat_joint_rz_fun %>%
  filter(
    Soil_habitat == "RZ",
    Treatment %in% c(
      "Natural",
      "BESE",
      "BARE"
    )
  ) %>%
  drop_na(
    Moisture,
    Temperature,
    pH,
    NH4,
    PO4,
    Redox
  )

# ============================================================
# 3. ENVIRONMENTAL TRANSFORMATIONS
# ============================================================

if (any(dat_joint_rz_fun$NH4 <= 0)) {
  dat_joint_rz_fun$NH4_log <-
    log1p(dat_joint_rz_fun$NH4)
} else {
  dat_joint_rz_fun$NH4_log <-
    log(dat_joint_rz_fun$NH4)
}

if (any(dat_joint_rz_fun$PO4 <= 0)) {
  dat_joint_rz_fun$PO4_log <-
    log1p(dat_joint_rz_fun$PO4)
} else {
  dat_joint_rz_fun$PO4_log <-
    log(dat_joint_rz_fun$PO4)
}

dat_joint_rz_fun <- dat_joint_rz_fun %>%
  mutate(
    Site = factor(
      Site,
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    ),
    
    Moisture_z =
      as.numeric(scale(Moisture)),
    
    Temperature_z =
      as.numeric(scale(Temperature)),
    
    pH_z =
      as.numeric(scale(pH)),
    
    NH4_log_z =
      as.numeric(scale(NH4_log)),
    
    PO4_log_z =
      as.numeric(scale(PO4_log)),
    
    Redox_z =
      as.numeric(scale(Redox))
  )

# ============================================================
# 4. COMMUNITY MATRIX
# ============================================================

ps_joint_rz_fun <- prune_samples(
  dat_joint_rz_fun$SampleID,
  ps_fun
)

ps_joint_rz_fun <- prune_taxa(
  taxa_sums(ps_joint_rz_fun) > 0,
  ps_joint_rz_fun
)

ps_joint_rz_fun_rel <-
  transform_sample_counts(
    ps_joint_rz_fun,
    function(x) x / sum(x)
  )

comm_joint_rz_fun <- as(
  otu_table(ps_joint_rz_fun_rel),
  "matrix"
)

if (taxa_are_rows(ps_joint_rz_fun_rel)) {
  comm_joint_rz_fun <-
    t(comm_joint_rz_fun)
}

comm_joint_rz_fun <-
  comm_joint_rz_fun[
    dat_joint_rz_fun$SampleID,
    ,
    drop = FALSE
  ]

comm_joint_rz_fun <-
  comm_joint_rz_fun[
    ,
    colSums(comm_joint_rz_fun) > 0,
    drop = FALSE
  ]

stopifnot(
  identical(
    rownames(comm_joint_rz_fun),
    dat_joint_rz_fun$SampleID
  )
)

cat("\n========================================\n")
cat("RHIZOSPHERE FUNGAL — JOINT db-RDA\n")
cat("========================================\n")

cat(
  "Samples:",
  nrow(dat_joint_rz_fun),
  "\n"
)

cat(
  "ASVs:",
  ncol(comm_joint_rz_fun),
  "\n\n"
)

print(
  table(
    dat_joint_rz_fun$Site,
    dat_joint_rz_fun$Treatment
  )
)

# ============================================================
# 5. H1
# ENVIRONMENT | SITE + ECOLOGICAL CONDITION
# ============================================================

db_env_unique_rz_fun <- vegan::dbrda(
  comm_joint_rz_fun ~
    Moisture_z +
    Temperature_z +
    pH_z +
    NH4_log_z +
    PO4_log_z +
    Redox_z +
    Condition(
      Site +
        Treatment
    ),
  data = dat_joint_rz_fun,
  distance = "bray",
  add = "lingoes"
)

env_unique_overall_rz_fun <-
  vegan::anova.cca(
    db_env_unique_rz_fun,
    permutations = 9999,
    strata = dat_joint_rz_fun$Site
  )

env_unique_margin_rz_fun <-
  vegan::anova.cca(
    db_env_unique_rz_fun,
    by = "margin",
    permutations = 9999,
    strata = dat_joint_rz_fun$Site
  )

env_unique_r2_rz_fun <-
  vegan::RsquareAdj(
    db_env_unique_rz_fun
  )

cat("\n========================================\n")
cat("A. ENVIRONMENT | SITE + CONDITION\n")
cat("========================================\n")

print(
  env_unique_overall_rz_fun
)

cat("\nMarginal environmental terms:\n")

print(
  env_unique_margin_rz_fun
)

cat("\nR2:\n")

print(
  env_unique_r2_rz_fun
)

cat("\nVIF:\n")

print(
  vegan::vif.cca(
    db_env_unique_rz_fun
  )
)

# ============================================================
# 6. BH CORRECTION
# ============================================================

env_margin_rz_fun_table <-
  as.data.frame(
    env_unique_margin_rz_fun
  )

env_margin_rz_fun_table$Variable <-
  rownames(
    env_margin_rz_fun_table
  )

rownames(
  env_margin_rz_fun_table
) <- NULL

env_margin_rz_fun_table <-
  env_margin_rz_fun_table %>%
  filter(
    Variable != "Residual"
  ) %>%
  mutate(
    p_BH = p.adjust(
      `Pr(>F)`,
      method = "BH"
    ),
    
    Significant_BH =
      p_BH < 0.05
  )

cat("\n========================================\n")
cat("BH-ADJUSTED MARGINAL ENVIRONMENTAL TESTS\n")
cat("========================================\n")

print(
  env_margin_rz_fun_table,
  row.names = FALSE
)

# ============================================================
# 7. H2
# ECOLOGICAL CONDITION | SITE + ENVIRONMENT
# ============================================================

db_condition_unique_rz_fun <- vegan::dbrda(
  comm_joint_rz_fun ~
    Treatment +
    Condition(
      Site +
        Moisture_z +
        Temperature_z +
        pH_z +
        NH4_log_z +
        PO4_log_z +
        Redox_z
    ),
  data = dat_joint_rz_fun,
  distance = "bray",
  add = "lingoes"
)

condition_unique_test_rz_fun <-
  vegan::anova.cca(
    db_condition_unique_rz_fun,
    permutations = 9999,
    strata = dat_joint_rz_fun$Site
  )

condition_unique_r2_rz_fun <-
  vegan::RsquareAdj(
    db_condition_unique_rz_fun
  )

cat("\n========================================\n")
cat("B. CONDITION | SITE + ENVIRONMENT\n")
cat("========================================\n")

print(
  condition_unique_test_rz_fun
)

cat("\nR2:\n")

print(
  condition_unique_r2_rz_fun
)

# ============================================================
# SAVE R3.20 — RECIPROCAL PARTIAL db-RDA RESULTS
#
# Saves:
# 1. Overall reciprocal partial db-RDA summary
# 2. Marginal environmental tests + BH correction
# 3. VIF diagnostics
# 4. Complete model objects for reproducibility
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(vegan)
})

# ------------------------------------------------------------
# 1. Check required objects exist
# ------------------------------------------------------------

required_objects <- c(
  
  # Bulk prokaryotic
  "db_env_unique",
  "env_unique_overall",
  "env_unique_margin",
  "env_unique_r2",
  "db_condition_unique",
  "condition_unique_test",
  "condition_unique_r2",
  
  # Bulk fungal
  "db_env_unique_fun",
  "env_unique_overall_fun",
  "env_unique_margin_fun",
  "env_unique_r2_fun",
  "db_condition_unique_fun",
  "condition_unique_test_fun",
  "condition_unique_r2_fun",
  
  # Rhizosphere prokaryotic
  "db_env_unique_rz_prok",
  "env_unique_overall_rz_prok",
  "env_unique_margin_rz_prok",
  "env_unique_r2_rz_prok",
  "db_condition_unique_rz_prok",
  "condition_unique_test_rz_prok",
  "condition_unique_r2_rz_prok",
  
  # Rhizosphere fungal
  "db_env_unique_rz_fun",
  "env_unique_overall_rz_fun",
  "env_unique_margin_rz_fun",
  "env_unique_r2_rz_fun",
  "db_condition_unique_rz_fun",
  "condition_unique_test_rz_fun",
  "condition_unique_r2_rz_fun"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1)
  )
]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects:\n",
    paste(
      missing_objects,
      collapse = "\n"
    )
  )
}

# ------------------------------------------------------------
# 2. Output folders
# ------------------------------------------------------------

dir.create(
  "environmental_association_outputs/R3_20_reciprocal_dbRDA",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "environmental_association_outputs/R3_20_reciprocal_dbRDA/models",
  recursive = TRUE,
  showWarnings = FALSE
)

out_dir <-
  "environmental_association_outputs/R3_20_reciprocal_dbRDA"

# ============================================================
# 3. Helper functions
# ============================================================

extract_overall <- function(
    test_object,
    r2_object,
    analysis,
    community,
    habitat,
    model_type
) {
  
  tab <- as.data.frame(
    test_object
  )
  
  tibble(
    Analysis = analysis,
    Community = community,
    Soil_habitat = habitat,
    Model = model_type,
    
    Df = tab$Df[1],
    
    SumOfSqs =
      tab$SumOfSqs[1],
    
    F =
      tab$F[1],
    
    p_value =
      tab$`Pr(>F)`[1],
    
    R2 =
      unname(
        r2_object$r.squared
      ),
    
    Adjusted_R2 =
      unname(
        r2_object$adj.r.squared
      )
  )
}


extract_marginal_env <- function(
    test_object,
    analysis,
    community,
    habitat
) {
  
  tab <- as.data.frame(
    test_object
  )
  
  tab$Variable <-
    rownames(tab)
  
  rownames(tab) <- NULL
  
  tab %>%
    as_tibble() %>%
    filter(
      Variable != "Residual"
    ) %>%
    mutate(
      Analysis = analysis,
      Community = community,
      Soil_habitat = habitat,
      
      p_adjust_BH =
        p.adjust(
          `Pr(>F)`,
          method = "BH"
        ),
      
      Significant_BH_0.05 =
        p_adjust_BH < 0.05,
      
      .before = 1
    ) %>%
    rename(
      p_value = `Pr(>F)`
    )
}


extract_vif <- function(
    db_object,
    analysis,
    community,
    habitat
) {
  
  vif_values <-
    vegan::vif.cca(
      db_object
    )
  
  tibble(
    Analysis = analysis,
    Community = community,
    Soil_habitat = habitat,
    Term = names(vif_values),
    VIF = as.numeric(
      vif_values
    ),
    High_VIF_10 =
      VIF >= 10
  )
}

# ============================================================
# 4. OVERALL RECIPROCAL PARTIAL db-RDA TABLE
# ============================================================

overall_R3_20 <- bind_rows(
  
  # ----------------------------------------------------------
  # Bulk prokaryotic
  # ----------------------------------------------------------
  
  extract_overall(
    env_unique_overall,
    env_unique_r2,
    "Bulk soil prokaryotic",
    "Prokaryotic",
    "BS",
    "Environment | Site + ecological condition"
  ),
  
  extract_overall(
    condition_unique_test,
    condition_unique_r2,
    "Bulk soil prokaryotic",
    "Prokaryotic",
    "BS",
    "Ecological condition | Site + environment"
  ),
  
  # ----------------------------------------------------------
  # Bulk fungal
  # ----------------------------------------------------------
  
  extract_overall(
    env_unique_overall_fun,
    env_unique_r2_fun,
    "Bulk soil fungal",
    "Fungal",
    "BS",
    "Environment | Site + ecological condition"
  ),
  
  extract_overall(
    condition_unique_test_fun,
    condition_unique_r2_fun,
    "Bulk soil fungal",
    "Fungal",
    "BS",
    "Ecological condition | Site + environment"
  ),
  
  # ----------------------------------------------------------
  # Rhizosphere prokaryotic
  # ----------------------------------------------------------
  
  extract_overall(
    env_unique_overall_rz_prok,
    env_unique_r2_rz_prok,
    "Rhizosphere prokaryotic",
    "Prokaryotic",
    "RZ",
    "Environment | Site + ecological condition"
  ),
  
  extract_overall(
    condition_unique_test_rz_prok,
    condition_unique_r2_rz_prok,
    "Rhizosphere prokaryotic",
    "Prokaryotic",
    "RZ",
    "Ecological condition | Site + environment"
  ),
  
  # ----------------------------------------------------------
  # Rhizosphere fungal
  # ----------------------------------------------------------
  
  extract_overall(
    env_unique_overall_rz_fun,
    env_unique_r2_rz_fun,
    "Rhizosphere fungal",
    "Fungal",
    "RZ",
    "Environment | Site + ecological condition"
  ),
  
  extract_overall(
    condition_unique_test_rz_fun,
    condition_unique_r2_rz_fun,
    "Rhizosphere fungal",
    "Fungal",
    "RZ",
    "Ecological condition | Site + environment"
  )
)

write_csv(
  overall_R3_20,
  file.path(
    out_dir,
    "Sx_R3_20_reciprocal_partial_dbRDA_overall.csv"
  )
)

# ============================================================
# 5. MARGINAL ENVIRONMENTAL TESTS
# BH correction performed separately within each analysis
# ============================================================

marginal_environment_R3_20 <- bind_rows(
  
  extract_marginal_env(
    env_unique_margin,
    "Bulk soil prokaryotic",
    "Prokaryotic",
    "BS"
  ),
  
  extract_marginal_env(
    env_unique_margin_fun,
    "Bulk soil fungal",
    "Fungal",
    "BS"
  ),
  
  extract_marginal_env(
    env_unique_margin_rz_prok,
    "Rhizosphere prokaryotic",
    "Prokaryotic",
    "RZ"
  ),
  
  extract_marginal_env(
    env_unique_margin_rz_fun,
    "Rhizosphere fungal",
    "Fungal",
    "RZ"
  )
)

write_csv(
  marginal_environment_R3_20,
  file.path(
    out_dir,
    "Sx_R3_20_environment_unique_marginal_tests_BH.csv"
  )
)

# ============================================================
# 6. VIF DIAGNOSTICS
# Based on Environment | Site + ecological condition model
# ============================================================

vif_R3_20 <- bind_rows(
  
  extract_vif(
    db_env_unique,
    "Bulk soil prokaryotic",
    "Prokaryotic",
    "BS"
  ),
  
  extract_vif(
    db_env_unique_fun,
    "Bulk soil fungal",
    "Fungal",
    "BS"
  ),
  
  extract_vif(
    db_env_unique_rz_prok,
    "Rhizosphere prokaryotic",
    "Prokaryotic",
    "RZ"
  ),
  
  extract_vif(
    db_env_unique_rz_fun,
    "Rhizosphere fungal",
    "Fungal",
    "RZ"
  )
)

write_csv(
  vif_R3_20,
  file.path(
    out_dir,
    "Sx_R3_20_reciprocal_partial_dbRDA_VIF.csv"
  )
)

# ============================================================
# 7. SHORT READER-FACING SUMMARY
# Convenient candidate for a Supplementary table
# ============================================================

summary_R3_20 <- overall_R3_20 %>%
  mutate(
    R2_percent =
      100 * R2,
    
    Adjusted_R2_percent =
      100 * Adjusted_R2
  ) %>%
  select(
    Analysis,
    Community,
    Soil_habitat,
    Model,
    Df,
    F,
    p_value,
    R2_percent,
    Adjusted_R2_percent
  )

write_csv(
  summary_R3_20,
  file.path(
    out_dir,
    "Sx_R3_20_reciprocal_partial_dbRDA_summary_table.csv"
  )
)

# ============================================================
# 8. SIGNIFICANT ENVIRONMENTAL VARIABLES ONLY
# Useful for manuscript checking
# ============================================================

significant_environment_R3_20 <-
  marginal_environment_R3_20 %>%
  filter(
    Significant_BH_0.05
  )

write_csv(
  significant_environment_R3_20,
  file.path(
    out_dir,
    "Sx_R3_20_significant_environmental_variables_BH.csv"
  )
)

# ============================================================
# 9. SAVE COMPLETE MODEL OBJECTS
# ============================================================

R3_20_model_objects <- list(
  
  Bulk_prokaryotic = list(
    environment_model =
      db_env_unique,
    
    environment_overall =
      env_unique_overall,
    
    environment_marginal =
      env_unique_margin,
    
    environment_R2 =
      env_unique_r2,
    
    condition_model =
      db_condition_unique,
    
    condition_overall =
      condition_unique_test,
    
    condition_R2 =
      condition_unique_r2
  ),
  
  Bulk_fungal = list(
    environment_model =
      db_env_unique_fun,
    
    environment_overall =
      env_unique_overall_fun,
    
    environment_marginal =
      env_unique_margin_fun,
    
    environment_R2 =
      env_unique_r2_fun,
    
    condition_model =
      db_condition_unique_fun,
    
    condition_overall =
      condition_unique_test_fun,
    
    condition_R2 =
      condition_unique_r2_fun
  ),
  
  Rhizosphere_prokaryotic = list(
    environment_model =
      db_env_unique_rz_prok,
    
    environment_overall =
      env_unique_overall_rz_prok,
    
    environment_marginal =
      env_unique_margin_rz_prok,
    
    environment_R2 =
      env_unique_r2_rz_prok,
    
    condition_model =
      db_condition_unique_rz_prok,
    
    condition_overall =
      condition_unique_test_rz_prok,
    
    condition_R2 =
      condition_unique_r2_rz_prok
  ),
  
  Rhizosphere_fungal = list(
    environment_model =
      db_env_unique_rz_fun,
    
    environment_overall =
      env_unique_overall_rz_fun,
    
    environment_marginal =
      env_unique_margin_rz_fun,
    
    environment_R2 =
      env_unique_r2_rz_fun,
    
    condition_model =
      db_condition_unique_rz_fun,
    
    condition_overall =
      condition_unique_test_rz_fun,
    
    condition_R2 =
      condition_unique_r2_rz_fun
  )
)

saveRDS(
  R3_20_model_objects,
  file.path(
    out_dir,
    "models",
    "R3_20_reciprocal_partial_dbRDA_all_models.rds"
  )
)

# ============================================================
# 10. PRINT CHECKS
# ============================================================

cat("\n========================================\n")
cat("R3.20 RECIPROCAL db-RDA RESULTS SAVED\n")
cat("========================================\n\n")

cat("Overall summary:\n")
print(
  summary_R3_20,
  n = Inf,
  width = Inf
)

cat("\nBH-significant environmental variables:\n")
print(
  significant_environment_R3_20,
  n = Inf,
  width = Inf
)

cat("\nSaved files:\n")

print(
  list.files(
    out_dir,
    recursive = TRUE
  )
)

