# ============================================================
# PCA of soil parameters and nutrients
# Treatment color is kept consistent:
# Natural = green, BESE = orange, Bare = dark blue, Degraded = violet
# Site is represented by point shape
# Dotted 95% confidence ellipses, no fill
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)

soil <- read_excel("metadata_env_parameters.xlsx") %>%
  filter(!is.na(SITE), !is.na(TREATMENT)) %>%
  mutate(
    SITE = factor(
      trimws(toupper(SITE)),
      levels = c("ACP", "BASS", "TRUG")
    ),
    TREATMENT = factor(
      trimws(toupper(TREATMENT)),
      levels = c("NATURAL", "BESE", "BARE", "DEGRADED")
    )
  )

pca_vars <- c(
  "Moisture_avg",
  "Temperature_avg",
  "pH_avg",
  "Redox_avg",
  "Salinity",
  "NO3",
  "NO2",
  "PO4",
  "NH4"
)

soil_pca <- soil %>%
  select(SITE, TREATMENT, all_of(pca_vars)) %>%
  na.omit()

pca_result <- prcomp(
  soil_pca[, pca_vars],
  center = TRUE,
  scale. = TRUE
)

scores <- as.data.frame(pca_result$x[, 1:2]) %>%
  bind_cols(
    soil_pca %>%
      select(SITE, TREATMENT)
  )

loadings <- as.data.frame(pca_result$rotation[, 1:2]) %>%
  mutate(
    Variable = rownames(pca_result$rotation)
  )

explained <- summary(pca_result)$importance[2, 1:2] * 100

treatment_colors <- c(
  "NATURAL" = "#2E8B57",
  "BESE" = "#E58A52",
  "BARE" = "#1F4E79",
  "DEGRADED" = "#8E63A9"
)

arrow_scale <- min(
  max(abs(scores$PC1)),
  max(abs(scores$PC2))
) * 0.95

loading_radius <- max(
  sqrt(loadings$PC1^2 + loadings$PC2^2)
)

loadings <- loadings %>%
  mutate(
    PC1_plot = PC1 / loading_radius * arrow_scale,
    PC2_plot = PC2 / loading_radius * arrow_scale
  )

p <- ggplot(scores, aes(x = PC1, y = PC2)) +
  stat_ellipse(
    aes(colour = TREATMENT),
    level = 0.95,
    type = "norm",
    linewidth = 0.8,
    linetype = "dotted",
    fill = NA,
    show.legend = TRUE
  ) +
  geom_point(
    aes(
      colour = TREATMENT,
      shape = SITE
    ),
    size = 3,
    stroke = 0.7
  ) +
  geom_segment(
    data = loadings,
    aes(
      x = 0,
      y = 0,
      xend = PC1_plot,
      yend = PC2_plot
    ),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.5,
    colour = "black"
  ) +
  geom_text(
    data = loadings,
    aes(
      x = PC1_plot * 1.06,
      y = PC2_plot * 1.06,
      label = Variable
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_colour_manual(
    values = treatment_colors,
    labels = c(
      "NATURAL" = "Natural",
      "BESE" = "BESE",
      "BARE" = "Bare",
      "DEGRADED" = "Degraded"
    )
  ) +
  scale_shape_manual(
    values = c(
      "ACP" = 16,
      "BASS" = 17,
      "TRUG" = 15
    )
  ) +
  labs(
    title = "PCA of Soil Parameters and Nutrients",
    x = paste0("PC1 (", round(explained[1], 1), "%)"),
    y = paste0("PC2 (", round(explained[2], 1), "%)"),
    colour = "Treatment",
    shape = "Site"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "right"
  )

p

ggsave(
  "soil_pca_treatment_palette_R.png",
  p,
  width = 10,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  "soil_pca_treatment_palette_R.pdf",
  p,
  width = 10,
  height = 7
)

write.csv(scores, "soil_pca_scores_R.csv", row.names = FALSE)
write.csv(loadings, "soil_pca_loadings_R.csv", row.names = FALSE)
