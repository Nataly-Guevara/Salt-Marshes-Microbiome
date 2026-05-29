#!/usr/bin/env Rscript

###############################################################################
# Vegetation coverage, colour, and height analysis
#
# Purpose:
# This script analyses vegetation cover, colour categories, and plant height
# across saltmarsh restoration treatments and sites.
#
# Run from the repository root:
# Rscript scripts/vegetation/vegetation_coverage_color_height_analysis.R
###############################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(FSA)
  library(ggplot2)
})

# -----------------------------
# Paths
# -----------------------------

base_path <- "."

input_file <- file.path(
  base_path,
  "metadata/vegetation/Vegetation_coverage_color_height.xlsx"
)

out_dir <- file.path(base_path, "tables/vegetation")
fig_dir <- file.path(base_path, "figures/vegetation")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Load data
# -----------------------------

veg <- read_excel(input_file)

veg <- veg %>%
  mutate(
    SITE = as.factor(SITE),
    TREATMENT = as.factor(TREATMENT),
    PLOT = as.factor(PLOT)
  )

# -----------------------------
# Calculate plot-level coverage
# -----------------------------

coverage_cols <- grep("_COV$", names(veg), value = TRUE)

veg <- veg %>%
  rowwise() %>%
  mutate(
    COVERAGE_mean = mean(c_across(all_of(coverage_cols)), na.rm = TRUE),
    COVERAGE_sd = sd(c_across(all_of(coverage_cols)), na.rm = TRUE)
  ) %>%
  ungroup()

# -----------------------------
# Summary statistics by site and treatment
# -----------------------------

summary_site_treatment <- veg %>%
  group_by(SITE, TREATMENT) %>%
  summarise(
    n = n(),
    Coverage_mean = mean(COVERAGE_mean, na.rm = TRUE),
    Coverage_sd = sd(COVERAGE_mean, na.rm = TRUE),
    Height_mean = mean(HEIGHT, na.rm = TRUE),
    Height_sd = sd(HEIGHT, na.rm = TRUE),
    Green_mean = mean(GREEN, na.rm = TRUE),
    Green_sd = sd(GREEN, na.rm = TRUE),
    Red_mean = mean(RED, na.rm = TRUE),
    Red_sd = sd(RED, na.rm = TRUE),
    Brown_mean = mean(BROWN, na.rm = TRUE),
    Brown_sd = sd(BROWN, na.rm = TRUE),
    .groups = "drop"
  )

write.xlsx(
  summary_site_treatment,
  file.path(out_dir, "vegetation_summary_by_site_treatment.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Helper functions
# -----------------------------

run_kruskal <- function(data, response, group_var) {
  formula <- as.formula(paste(response, "~", group_var))
  test <- kruskal.test(formula, data = data)
  
  data.frame(
    Variable = response,
    Grouping = group_var,
    Statistic = as.numeric(test$statistic),
    df = as.numeric(test$parameter),
    P_value = test$p.value
  )
}

run_dunn <- function(data, response, group_var) {
  formula <- as.formula(paste(response, "~", group_var))
  dunnTest(formula, data = data, method = "bonferroni")$res %>%
    mutate(
      Variable = response,
      Grouping = group_var
    )
}

response_vars <- c(
  "COVERAGE_mean",
  "HEIGHT",
  "GREEN",
  "RED",
  "BROWN"
)

# -----------------------------
# Kruskal-Wallis tests by treatment
# -----------------------------

kruskal_treatment <- bind_rows(
  lapply(response_vars, function(x) run_kruskal(veg, x, "TREATMENT"))
)

write.xlsx(
  kruskal_treatment,
  file.path(out_dir, "vegetation_kruskal_by_treatment.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Dunn tests by treatment
# -----------------------------

dunn_treatment <- bind_rows(
  lapply(response_vars, function(x) run_dunn(veg, x, "TREATMENT"))
)

write.xlsx(
  dunn_treatment,
  file.path(out_dir, "vegetation_dunn_by_treatment.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Kruskal-Wallis tests by site
# -----------------------------

kruskal_site <- bind_rows(
  lapply(response_vars, function(x) run_kruskal(veg, x, "SITE"))
)

write.xlsx(
  kruskal_site,
  file.path(out_dir, "vegetation_kruskal_by_site.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Dunn tests by site
# -----------------------------

dunn_site <- bind_rows(
  lapply(response_vars, function(x) run_dunn(veg, x, "SITE"))
)

write.xlsx(
  dunn_site,
  file.path(out_dir, "vegetation_dunn_by_site.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Significant post hoc results only
# -----------------------------

significant_dunn_treatment <- dunn_treatment %>%
  filter(P.adj < 0.05)

significant_dunn_site <- dunn_site %>%
  filter(P.adj < 0.05)

write.xlsx(
  significant_dunn_treatment,
  file.path(out_dir, "vegetation_significant_dunn_by_treatment.xlsx"),
  overwrite = TRUE
)

write.xlsx(
  significant_dunn_site,
  file.path(out_dir, "vegetation_significant_dunn_by_site.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Long-format table for plotting
# -----------------------------

veg_long <- veg %>%
  select(SITE, TREATMENT, PLOT, COVERAGE_mean, HEIGHT, GREEN, RED, BROWN) %>%
  pivot_longer(
    cols = c(COVERAGE_mean, HEIGHT, GREEN, RED, BROWN),
    names_to = "Variable",
    values_to = "Value"
  )

write.xlsx(
  veg_long,
  file.path(out_dir, "vegetation_long_format_for_plots.xlsx"),
  overwrite = TRUE
)

# -----------------------------
# Boxplots by treatment
# -----------------------------

p_treatment <- ggplot(
  veg_long,
  aes(x = TREATMENT, y = Value, fill = TREATMENT)
) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Vegetation variables by treatment",
    x = "Treatment",
    y = "Value"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "vegetation_variables_by_treatment.png"),
  plot = p_treatment,
  width = 9,
  height = 6,
  dpi = 300
)

# -----------------------------
# Boxplots by site
# -----------------------------

p_site <- ggplot(
  veg_long,
  aes(x = SITE, y = Value, fill = SITE)
) +
  geom_boxplot() +
  facet_wrap(~ Variable, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Vegetation variables by site",
    x = "Site",
    y = "Value"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "vegetation_variables_by_site.png"),
  plot = p_site,
  width = 9,
  height = 6,
  dpi = 300
)

# -----------------------------
# Heatmap of mean vegetation variables by site and treatment
# -----------------------------

heatmap_data <- veg_long %>%
  group_by(SITE, TREATMENT, Variable) %>%
  summarise(
    Mean_value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

write.xlsx(
  heatmap_data,
  file.path(out_dir, "vegetation_heatmap_summary.xlsx"),
  overwrite = TRUE
)

p_heatmap <- ggplot(
  heatmap_data,
  aes(x = TREATMENT, y = SITE, fill = Mean_value)
) +
  geom_tile() +
  facet_wrap(~ Variable, scales = "free") +
  theme_minimal() +
  labs(
    title = "Mean vegetation variables by site and treatment",
    x = "Treatment",
    y = "Site",
    fill = "Mean value"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(fig_dir, "vegetation_heatmap_by_site_treatment.png"),
  plot = p_heatmap,
  width = 9,
  height = 6,
  dpi = 300
)

cat("Vegetation coverage, colour, and height analysis completed.\n")
cat("Tables saved to:", out_dir, "\n")
cat("Figures saved to:", fig_dir, "\n")
