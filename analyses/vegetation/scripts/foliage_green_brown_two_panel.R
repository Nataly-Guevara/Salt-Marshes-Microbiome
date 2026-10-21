# ============================================================
# Green and brown foliage in BESE and BARE plots
# Individual plots + treatment means + stratified bootstrap CIs
# ============================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(patchwork)

vegetation <- read_excel(
  "Vegetation_coverage_color_height.xlsx"
) %>%
  rename(
    TREATMENT = `TREATMENT `
  ) %>%
  mutate(
    SITE = trimws(toupper(SITE)),
    TREATMENT = trimws(toupper(TREATMENT))
  ) %>%
  filter(
    TREATMENT %in% c("BESE", "BARE")
  ) %>%
  mutate(
    SITE = factor(
      SITE,
      levels = c("ACP", "BASS", "TRUG")
    ),
    TREATMENT = factor(
      TREATMENT,
      levels = c("BESE", "BARE")
    )
  )

# Site-stratified bootstrap:
# plots are resampled within each Site x Treatment combination.
bootstrap_mean_ci <- function(
    data,
    response,
    n_boot = 20000,
    seed = 123
) {

  set.seed(seed)

  treatments <- levels(data$TREATMENT)

  map_dfr(
    treatments,
    function(treatment_name) {

      treatment_data <- data %>%
        filter(
          TREATMENT == treatment_name
        )

      observed_mean <- mean(
        treatment_data[[response]],
        na.rm = TRUE
      )

      bootstrap_means <- replicate(
        n_boot,
        {
          treatment_data %>%
            group_by(SITE) %>%
            slice_sample(
              n = n(),
              replace = TRUE
            ) %>%
            ungroup() %>%
            summarise(
              value = mean(
                .data[[response]],
                na.rm = TRUE
              )
            ) %>%
            pull(value)
        }
      )

      tibble(
        TREATMENT = factor(
          treatment_name,
          levels = treatments
        ),
        estimate = observed_mean,
        lower_ci = quantile(
          bootstrap_means,
          0.025,
          na.rm = TRUE
        ),
        upper_ci = quantile(
          bootstrap_means,
          0.975,
          na.rm = TRUE
        )
      )
    }
  )
}

make_foliage_panel <- function(
    data,
    response,
    title,
    y_label,
    seed
) {

  plot_stats <- bootstrap_mean_ci(
    data = data,
    response = response,
    seed = seed
  )

  ggplot(
    data,
    aes(
      x = TREATMENT,
      y = .data[[response]]
    )
  ) +
    geom_point(
      aes(
        shape = SITE,
        group = SITE
      ),
      position = position_jitterdodge(
        jitter.width = 0.05,
        jitter.height = 0,
        dodge.width = 0.28,
        seed = seed
      ),
      size = 2.7,
      alpha = 0.75
    ) +
    geom_errorbar(
      data = plot_stats,
      aes(
        x = TREATMENT,
        ymin = lower_ci,
        ymax = upper_ci
      ),
      inherit.aes = FALSE,
      width = 0.10,
      linewidth = 0.8
    ) +
    geom_point(
      data = plot_stats,
      aes(
        x = TREATMENT,
        y = estimate
      ),
      inherit.aes = FALSE,
      shape = 23,
      size = 4,
      stroke = 1
    ) +
    scale_shape_manual(
      values = c(
        "ACP" = 16,
        "BASS" = 17,
        "TRUG" = 15
      )
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      expand = expansion(
        mult = c(0.02, 0.04)
      )
    ) +
    labs(
      title = title,
      x = NULL,
      y = y_label,
      shape = "Site"
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      legend.position = "bottom"
    )
}

green_plot <- make_foliage_panel(
  vegetation,
  response = "GREEN",
  title = "Green foliage",
  y_label = "Proportion of green foliage",
  seed = 123
)

brown_plot <- make_foliage_panel(
  vegetation,
  response = "BROWN",
  title = "Brown foliage",
  y_label = "Proportion of brown foliage",
  seed = 456
)

foliage_figure <- (
  green_plot + brown_plot +
    plot_layout(
      guides = "collect"
    ) +
    plot_annotation(
      tag_levels = "A"
    )
) &
  theme(
    legend.position = "bottom"
  )

foliage_figure

ggsave(
  "foliage_green_brown_two_panel_R.png",
  foliage_figure,
  width = 8,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  "foliage_green_brown_two_panel_R.pdf",
  foliage_figure,
  width = 8,
  height = 4.5
)
