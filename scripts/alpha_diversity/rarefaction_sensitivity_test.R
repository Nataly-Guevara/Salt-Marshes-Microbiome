#!/usr/bin/env Rscript

##############################################################################
# RAREFACTION SENSITIVITY WORKFLOW
#
# Purpose:
# This script evaluates whether alpha-diversity patterns are robust to
# differences in sequencing depth. Rarefaction is used here only as a
# sensitivity analysis and not as the main inferential approach.
#
# Raw reads are available from NCBI SRA BioProject PRJNA1446205.
#
# Run from the repository root:
# Rscript scripts/alpha_diversity/rarefaction_sensitivity_test.R
##############################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(readr)
  library(purrr)
  library(vegan)
})

# ---------------------------
# Paths and user inputs
# ---------------------------
base_path <- "."

outdir <- file.path(base_path, "tables/alpha_diversity/rarefaction_depth_sensitivity")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ps_16S_file <- file.path(base_path, "metadata/dada2/phyloseq_16S.rds")
ps_ITS_file <- file.path(base_path, "metadata/dada2/phyloseq_ITS.rds")

# Candidate depths to test
candidate_depths <- c(10000, 15000, 20000, 25000, 30000, 40000, 50000)

# Set seed for reproducibility
set.seed(123)

# ---------------------------
# Helper: depth summaries
# ---------------------------
make_depth_tables <- function(ps, dataset_name, candidate_depths) {
  depths <- sample_sums(ps)
  
  summary_tbl <- tibble(
    Dataset = dataset_name,
    Min = min(depths),
    Q1 = as.numeric(quantile(depths, 0.25)),
    Median = median(depths),
    Mean = mean(depths),
    Q3 = as.numeric(quantile(depths, 0.75)),
    Max = max(depths),
    N = length(depths)
  )
  
  retention_tbl <- tibble(
    Dataset = dataset_name,
    RarefactionDepth = candidate_depths,
    SamplesRetained = map_int(candidate_depths, ~ sum(depths >= .x)),
    SamplesDropped = length(depths) - SamplesRetained,
    ProportionRetained = SamplesRetained / length(depths)
  )
  
  list(summary = summary_tbl, retention = retention_tbl, depths = depths)
}

# ---------------------------
# Helper: histogram
# ---------------------------
plot_depth_histogram <- function(depths, dataset_name) {
  ggplot(tibble(Depth = depths), aes(x = Depth)) +
    geom_histogram(bins = 30, fill = "grey70", color = "black") +
    labs(
      title = paste(dataset_name, "sequencing-depth distribution"),
      x = "Reads per sample",
      y = "Number of samples"
    ) +
    theme_bw(base_size = 12)
}

# ---------------------------
# Helper: retention plot
# ---------------------------
plot_retention_curve <- function(retention_tbl, dataset_name) {
  ggplot(retention_tbl, aes(x = RarefactionDepth, y = SamplesRetained)) +
    geom_line() +
    geom_point(size = 2) +
    labs(
      title = paste(dataset_name, "sample retention across candidate rarefaction depths"),
      x = "Rarefaction depth",
      y = "Samples retained"
    ) +
    theme_bw(base_size = 12)
}

# ---------------------------
# Helper: rarefaction curves
# vegan::rarecurve works on samples as rows and taxa as columns.
# ---------------------------
make_rarefaction_curve_plot <- function(ps, dataset_name, max_points = 40) {
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  
  max_depth <- max(rowSums(otu))
  step_size <- max(1, floor(max_depth / max_points))
  
  png(
    filename = file.path(outdir, paste0(dataset_name, "_rarefaction_curves.png")),
    width = 1400, height = 1000, res = 150
  )
  rarecurve(
    otu,
    step = step_size,
    sample = min(max_depth, max(rowSums(otu))),
    xlab = "Sequencing depth",
    ylab = "Observed ASVs",
    label = FALSE,
    main = paste(dataset_name, "rarefaction curves")
  )
  dev.off()
}

# ---------------------------
# Helper: run one alpha model
# ---------------------------
fit_alpha_models <- function(alpha_df, metric, status_var, compartment_var, site_var) {
  form <- as.formula(
    paste0(metric, " ~ ", status_var, " * ", compartment_var, " + (1|", site_var, ")")
  )
  
  model <- lmer(form, data = alpha_df)
  
  anova_tbl <- as.data.frame(anova(model)) %>%
    rownames_to_column("Term")
  
  emm <- emmeans(model, specs = as.formula(
    paste0("pairwise ~ ", status_var, " | ", compartment_var)
  ))
  
  emmeans_tbl <- as.data.frame(emm$emmeans)
  contrasts_tbl <- as.data.frame(emm$contrasts)
  
  list(
    model = model,
    anova = anova_tbl,
    emmeans = emmeans_tbl,
    contrasts = contrasts_tbl
  )
}

# ---------------------------
# Helper: analyse one dataset at one rarefaction depth
# ---------------------------
analyse_one_depth <- function(ps,
                              dataset_name,
                              rare_depth,
                              status_var = "EcologicalStatus",
                              compartment_var = "SoilPortion",
                              site_var = "Site",
                              rngseed = 123) {
  
  depths <- sample_sums(ps)
  n_before <- length(depths)
  n_after <- sum(depths >= rare_depth)
  
  ps_rare <- rarefy_even_depth(
    ps,
    sample.size = rare_depth,
    rngseed = rngseed,
    replace = FALSE,
    trimOTUs = TRUE,
    verbose = FALSE
  )
  
  alpha_raw <- estimate_richness(ps, measures = c("Chao1", "Shannon", "Simpson")) %>%
    cbind(as.data.frame(sample_data(ps))) %>%
    as_tibble() %>%
    mutate(Data = "Raw",
           Dataset = dataset_name)
  
  alpha_rare <- estimate_richness(ps_rare, measures = c("Chao1", "Shannon", "Simpson")) %>%
    cbind(as.data.frame(sample_data(ps_rare))) %>%
    as_tibble() %>%
    mutate(Data = paste0("Rarefied_", rare_depth),
           Dataset = dataset_name)
  
  summary_tbl <- bind_rows(alpha_raw, alpha_rare) %>%
    group_by(Dataset, Data, .data[[status_var]], .data[[compartment_var]]) %>%
    summarise(
      n = n(),
      mean_Chao1 = mean(Chao1, na.rm = TRUE),
      sd_Chao1 = sd(Chao1, na.rm = TRUE),
      mean_Shannon = mean(Shannon, na.rm = TRUE),
      sd_Shannon = sd(Shannon, na.rm = TRUE),
      mean_Simpson = mean(Simpson, na.rm = TRUE),
      sd_Simpson = sd(Simpson, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(RarefactionDepth = rare_depth)
  
  # Models are repeated for Chao1 and Shannon because these were the main
  # rarefaction-sensitivity metrics used to evaluate robustness of alpha diversity.
  chao_raw  <- fit_alpha_models(alpha_raw,  "Chao1",   status_var, compartment_var, site_var)
  chao_rare <- fit_alpha_models(alpha_rare, "Chao1",   status_var, compartment_var, site_var)
  shan_raw  <- fit_alpha_models(alpha_raw,  "Shannon", status_var, compartment_var, site_var)
  shan_rare <- fit_alpha_models(alpha_rare, "Shannon", status_var, compartment_var, site_var)
  
  anova_tbl <- bind_rows(
    chao_raw$anova  %>% mutate(Metric = "Chao1",   Data = "Raw"),
    chao_rare$anova %>% mutate(Metric = "Chao1",   Data = paste0("Rarefied_", rare_depth)),
    shan_raw$anova  %>% mutate(Metric = "Shannon", Data = "Raw"),
    shan_rare$anova %>% mutate(Metric = "Shannon", Data = paste0("Rarefied_", rare_depth))
  ) %>%
    mutate(Dataset = dataset_name,
           RarefactionDepth = rare_depth,
           SamplesBefore = n_before,
           SamplesAfter = n_after)
  
  emmeans_tbl <- bind_rows(
    chao_raw$emmeans  %>% mutate(Metric = "Chao1",   Data = "Raw"),
    chao_rare$emmeans %>% mutate(Metric = "Chao1",   Data = paste0("Rarefied_", rare_depth)),
    shan_raw$emmeans  %>% mutate(Metric = "Shannon", Data = "Raw"),
    shan_rare$emmeans %>% mutate(Metric = "Shannon", Data = paste0("Rarefied_", rare_depth))
  ) %>%
    mutate(Dataset = dataset_name,
           RarefactionDepth = rare_depth)
  
  contrasts_tbl <- bind_rows(
    chao_raw$contrasts  %>% mutate(Metric = "Chao1",   Data = "Raw"),
    chao_rare$contrasts %>% mutate(Metric = "Chao1",   Data = paste0("Rarefied_", rare_depth)),
    shan_raw$contrasts  %>% mutate(Metric = "Shannon", Data = "Raw"),
    shan_rare$contrasts %>% mutate(Metric = "Shannon", Data = paste0("Rarefied_", rare_depth))
  ) %>%
    mutate(Dataset = dataset_name,
           RarefactionDepth = rare_depth)
  
  list(
    summary = summary_tbl,
    anova = anova_tbl,
    emmeans = emmeans_tbl,
    contrasts = contrasts_tbl
  )
}

# ---------------------------
# Helper: run one dataset across depths
# ---------------------------
run_dataset_workflow <- function(ps,
                                 dataset_name,
                                 candidate_depths,
                                 status_var = "EcologicalStatus",
                                 compartment_var = "SoilPortion",
                                 site_var = "Site",
                                 min_samples_required = 20) {
  
  depth_info <- make_depth_tables(ps, dataset_name, candidate_depths)
  write_csv(depth_info$summary,
            file.path(outdir, paste0(dataset_name, "_depth_summary.csv")))
  write_csv(depth_info$retention,
            file.path(outdir, paste0(dataset_name, "_depth_retention.csv")))
  
  p_hist <- plot_depth_histogram(depth_info$depths, dataset_name)
  ggsave(
    filename = file.path(outdir, paste0(dataset_name, "_depth_histogram.png")),
    plot = p_hist, width = 7, height = 5, dpi = 300
  )
  
  p_ret <- plot_retention_curve(depth_info$retention, dataset_name)
  ggsave(
    filename = file.path(outdir, paste0(dataset_name, "_retention_curve.png")),
    plot = p_ret, width = 7, height = 5, dpi = 300
  )
  
  make_rarefaction_curve_plot(ps, dataset_name)
  
  valid_depths <- depth_info$retention %>%
    filter(SamplesRetained >= min_samples_required) %>%
    pull(RarefactionDepth)
  
  results_list <- map(valid_depths, function(d) {
    message("Running ", dataset_name, " at depth ", d)
    analyse_one_depth(
      ps = ps,
      dataset_name = dataset_name,
      rare_depth = d,
      status_var = status_var,
      compartment_var = compartment_var,
      site_var = site_var,
      rngseed = 123
    )
  })
  
  combined_summary <- bind_rows(map(results_list, "summary"))
  combined_anova <- bind_rows(map(results_list, "anova"))
  combined_emmeans <- bind_rows(map(results_list, "emmeans"))
  combined_contrasts <- bind_rows(map(results_list, "contrasts"))
  
  write_csv(combined_summary,
            file.path(outdir, paste0(dataset_name, "_alpha_group_summaries_all_depths.csv")))
  write_csv(combined_anova,
            file.path(outdir, paste0(dataset_name, "_anova_all_depths.csv")))
  write_csv(combined_emmeans,
            file.path(outdir, paste0(dataset_name, "_emmeans_all_depths.csv")))
  write_csv(combined_contrasts,
            file.path(outdir, paste0(dataset_name, "_contrasts_all_depths.csv")))
  
  interaction_tbl <- combined_anova %>%
    filter(Term == paste0(status_var, ":", compartment_var),
           Metric %in% c("Chao1", "Shannon")) %>%
    select(Dataset, RarefactionDepth, Metric, Data, Term, `F value`, `Pr(>F)`,
           SamplesBefore, SamplesAfter)
  
  write_csv(interaction_tbl,
            file.path(outdir, paste0(dataset_name, "_interaction_stability_table.csv")))
  
  list(
    depth_summary = depth_info$summary,
    depth_retention = depth_info$retention,
    group_summaries = combined_summary,
    anova = combined_anova,
    emmeans = combined_emmeans,
    contrasts = combined_contrasts,
    interaction_stability = interaction_tbl
  )
}

##############################################################################
# RUN BOTH DATASETS
##############################################################################

ps_16S <- readRDS(ps_16S_file)
ps_ITS <- readRDS(ps_ITS_file)

res_16S <- run_dataset_workflow(
  ps = ps_16S,
  dataset_name = "16S",
  candidate_depths = candidate_depths,
  status_var = "EcologicalStatus",
  compartment_var = "SoilPortion",
  site_var = "Site"
)

res_ITS <- run_dataset_workflow(
  ps = ps_ITS,
  dataset_name = "ITS",
  candidate_depths = candidate_depths,
  status_var = "EcologicalStatus",
  compartment_var = "SoilPortion",
  site_var = "Site"
)

##############################################################################
# COMBINED TABLES
##############################################################################

write_csv(
  bind_rows(res_16S$depth_retention, res_ITS$depth_retention),
  file.path(outdir, "combined_depth_retention.csv")
)

write_csv(
  bind_rows(res_16S$interaction_stability, res_ITS$interaction_stability),
  file.path(outdir, "combined_interaction_stability.csv")
)

write_csv(
  bind_rows(res_16S$anova, res_ITS$anova),
  file.path(outdir, "combined_anova_all_depths.csv")
)

write_csv(
  bind_rows(res_16S$contrasts, res_ITS$contrasts),
  file.path(outdir, "combined_contrasts_all_depths.csv")
)

##############################################################################
# Candidate-depth justification table
##############################################################################

justification_tbl <- bind_rows(res_16S$depth_retention, res_ITS$depth_retention) %>%
  pivot_wider(
    names_from = Dataset,
    values_from = c(SamplesRetained, SamplesDropped, ProportionRetained)
  )

write_csv(
  justification_tbl,
  file.path(outdir, "candidate_depth_justification_table.csv")
)

cat("Rarefaction sensitivity workflow completed.\n")
cat("Outputs saved to:", outdir, "\n")
