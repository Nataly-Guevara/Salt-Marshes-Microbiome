# ============================================================
# M6 — VEGETATION DESCRIPTIVE FIGURE + PLANT-CONDITION TABLE
# ============================================================

# -----------------------------
# 1. Packages
# -----------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(emmeans)
library(purrr)
library(patchwork)


# -----------------------------
# 2. Import vegetation data
# -----------------------------

vegetation_file <- "Vegetation data for R(1) 2 copy.xlsx"

veg_raw <- readxl::read_excel(
  vegetation_file,
  sheet = 1
)

# Clean column names without needing janitor
names(veg_raw) <- names(veg_raw) %>%
  toupper() %>%
  gsub("[^A-Z0-9]+", "_", .) %>%
  gsub("^_|_$", "", .)


# -----------------------------
# 3. Helper functions
# -----------------------------

parse_number_safe <- function(x) {
  as.numeric(
    gsub(
      "[^0-9\\.\\-]",
      "",
      as.character(x)
    )
  )
}

find_col <- function(names_vector, possible_names) {
  
  hits <- names_vector[
    names_vector %in% possible_names
  ]
  
  if (length(hits) == 0) {
    stop(
      "None of these columns were found: ",
      paste(possible_names, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(names_vector, collapse = ", ")
    )
  }
  
  hits[1]
}


# -----------------------------
# 4. Identify key columns
# -----------------------------

site_col <- find_col(
  names(veg_raw),
  c("SITE")
)

treatment_col <- find_col(
  names(veg_raw),
  c("TREATMENT", "TRATAMIENTO")
)

plot_col <- find_col(
  names(veg_raw),
  c("PLOT", "PLOT_NO", "NO", "BLOCK", "BLOQUE")
)

height_col <- find_col(
  names(veg_raw),
  c("HEIGHT", "MAX_HEIGHT", "MAX_HEIGHT_MM", "MAXIMUM_HEIGHT")
)

green_col <- find_col(
  names(veg_raw),
  c("GREEN", "GREEN_PERCENT", "GREEN_PCT")
)

brown_col <- find_col(
  names(veg_raw),
  c("BROWN", "BROWN_PERCENT", "BROWN_PCT")
)

red_col <- find_col(
  names(veg_raw),
  c("RED", "RED_PERCENT", "RED_PCT")
)

coverage_cols <- grep(
  "COV|COVER",
  names(veg_raw),
  value = TRUE
)

if (length(coverage_cols) == 0) {
  stop(
    "No vegetation coverage columns were found. ",
    "I expected columns such as A1_COV, A2_COV, etc."
  )
}


# -----------------------------
# 5. Create plot-level vegetation table
# -----------------------------

veg_plot <- veg_raw %>%
  
  mutate(
    
    SITE =
      toupper(
        trimws(
          as.character(.data[[site_col]])
        )
      ),
    
    TREATMENT =
      toupper(
        trimws(
          as.character(.data[[treatment_col]])
        )
      ),
    
    PLOT =
      as.integer(
        parse_number_safe(.data[[plot_col]])
      ),
    
    HEIGHT =
      parse_number_safe(.data[[height_col]]),
    
    GREEN =
      parse_number_safe(.data[[green_col]]),
    
    BROWN =
      parse_number_safe(.data[[brown_col]]),
    
    RED =
      parse_number_safe(.data[[red_col]])
  ) %>%
  
 # insert 
  
  mutate(
    
    # Convert colour values to percentages if they were entered as proportions
    GREEN_pct =
      if (max(GREEN, na.rm = TRUE) <= 1) {
        GREEN * 100
      } else {
        GREEN
      },
    
    BROWN_pct =
      if (max(BROWN, na.rm = TRUE) <= 1) {
        BROWN * 100
      } else {
        BROWN
      },
    
    RED_pct =
      if (max(RED, na.rm = TRUE) <= 1) {
        RED * 100
      } else {
        RED
      }
  )
  
  
# Convert all coverage columns to numeric
veg_plot <- veg_plot %>%
  
  mutate(
    across(
      all_of(coverage_cols),
      parse_number_safe
    )
  ) %>%
  
  rowwise() %>%
  
  mutate(
    Mean_cover_pct =
      mean(
        c_across(all_of(coverage_cols)),
        na.rm = TRUE
      )
  ) %>%
  
  ungroup()


# -----------------------------
# 6. Treat BASS-BESE Plot 2 height = 0 as missing
# -----------------------------

veg_plot <- veg_plot %>%
  
  mutate(
    
    HEIGHT_original =
      HEIGHT,
    
    HEIGHT =
      if_else(
        SITE == "BASS" &
          TREATMENT == "BESE" &
          PLOT == 2 &
          HEIGHT == 0,
        NA_real_,
        HEIGHT
      ),
    
    Height_flag =
      if_else(
        SITE == "BASS" &
          TREATMENT == "BESE" &
          PLOT == 2,
        "Height treated as missing; recorded zero inconsistent with cover/foliage",
        NA_character_
      )
  )


# -----------------------------
# 7. Keep restoration treatments only
# -----------------------------

veg_restore <- veg_plot %>%
  
  filter(
    TREATMENT %in% c("BESE", "BARE")
  ) %>%
  
  mutate(
    
    SITE =
      factor(
        SITE,
        levels = c("ACP", "BASS", "TRUG")
      ),
    
    TREATMENT =
      factor(
        TREATMENT,
        levels = c("BESE", "BARE")
      ),
    
    PlotID =
      paste(
        SITE,
        TREATMENT,
        PLOT,
        sep = "_"
      )
  )


# -----------------------------
# 8. Quick audit
# -----------------------------

cat("\n========================================\n")
cat("VEGETATION DATA AUDIT\n")
cat("========================================\n")

print(
  veg_restore %>%
    count(SITE, TREATMENT)
)

cat("\nMissing values:\n")

print(
  veg_restore %>%
    summarise(
      missing_cover = sum(is.na(Mean_cover_pct)),
      missing_height = sum(is.na(HEIGHT)),
      missing_green = sum(is.na(GREEN_pct)),
      missing_brown = sum(is.na(BROWN_pct))
    )
)

cat("\nFlagged records:\n")

print(
  veg_restore %>%
    filter(!is.na(Height_flag)) %>%
    select(SITE, TREATMENT, PLOT, HEIGHT_original, HEIGHT, Height_flag)
)


# -----------------------------
# 9. Descriptive statistics
# -----------------------------

vegetation_summary <- veg_restore %>%
  
  group_by(TREATMENT) %>%
  
  summarise(
    n_cover =
      sum(!is.na(Mean_cover_pct)),
    
    mean_cover_pct =
      mean(Mean_cover_pct, na.rm = TRUE),
    
    sd_cover_pct =
      sd(Mean_cover_pct, na.rm = TRUE),
    
    n_height =
      sum(!is.na(HEIGHT)),
    
    mean_height =
      mean(HEIGHT, na.rm = TRUE),
    
    sd_height =
      sd(HEIGHT, na.rm = TRUE),
    
    n_green =
      sum(!is.na(GREEN_pct)),
    
    mean_green_pct =
      mean(GREEN_pct, na.rm = TRUE),
    
    sd_green_pct =
      sd(GREEN_pct, na.rm = TRUE),
    
    n_brown =
      sum(!is.na(BROWN_pct)),
    
    mean_brown_pct =
      mean(BROWN_pct, na.rm = TRUE),
    
    sd_brown_pct =
      sd(BROWN_pct, na.rm = TRUE),
    
    .groups = "drop"
  )

print(vegetation_summary)

write.csv(
  vegetation_summary,
  "M6_vegetation_summary_BESE_BARE.csv",
  row.names = FALSE
)


# -----------------------------
# 10. Function for vegetation panels
# -----------------------------

make_veg_panel <- function(data,
                           response,
                           y_label,
                           title_text) {
  
  plot_data <- data %>%
    
    filter(
      !is.na(.data[[response]]),
      !is.na(SITE),
      !is.na(TREATMENT)
    )
  
  model <- lm(
    reformulate(
      c("SITE", "TREATMENT"),
      response = response
    ),
    data = plot_data
  )
  
  emm <- emmeans(
    model,
    ~ TREATMENT
  )
  
  emm_ci <- as.data.frame(
    confint(emm)
  ) %>%
    
    transmute(
      TREATMENT,
      estimate = emmean,
      lower_ci = lower.CL,
      upper_ci = upper.CL
    )
  
  pairwise <- as.data.frame(
    pairs(
      emm,
      adjust = "tukey"
    )
  )
  
  p_pair <- pairwise$p.value[1]
  
  letter_table <- emm_ci %>%
    
    arrange(desc(estimate)) %>%
    
    mutate(
      letter =
        if (p_pair < 0.05) {
          c("a", "b")
        } else {
          c("a", "a")
        }
    ) %>%
    
    select(TREATMENT, letter)
  
  plot_stats <- emm_ci %>%
    
    left_join(
      letter_table,
      by = "TREATMENT"
    )
  
  y_range <- diff(
    range(
      c(
        plot_data[[response]],
        plot_stats$lower_ci,
        plot_stats$upper_ci
      ),
      na.rm = TRUE
    )
  )
  
  plot_stats <- plot_stats %>%
    
    mutate(
      letter_y =
        upper_ci + 0.12 * y_range
    )
  
  panel <- ggplot(
    plot_data,
    aes(
      x = TREATMENT,
      y = .data[[response]]
    )
  ) +
    
    geom_point(
      aes(
        fill = TREATMENT,
        shape = SITE
      ),
      position =
        position_jitter(
          width = 0.08,
          height = 0,
          seed = 123
        ),
      size = 3,
      alpha = 0.85,
      colour = "black",
      stroke = 0.45
    ) +
    
    geom_errorbar(
      data = plot_stats,
      aes(
        x = TREATMENT,
        ymin = lower_ci,
        ymax = upper_ci
      ),
      inherit.aes = FALSE,
      width = 0.12,
      linewidth = 0.8,
      colour = "black"
    ) +
    
    geom_point(
      data = plot_stats,
      aes(
        x = TREATMENT,
        y = estimate
      ),
      inherit.aes = FALSE,
      shape = 23,
      size = 4.2,
      fill = "white",
      colour = "black",
      stroke = 1
    ) +
    
    geom_text(
      data = plot_stats,
      aes(
        x = TREATMENT,
        y = letter_y,
        label = letter
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 4.5
    ) +
    
    scale_fill_manual(
      values = c(
        "BESE" = "#FF5A5F",
        "BARE" = "#FFA630"
      ),
      drop = FALSE
    ) +
    
    scale_shape_manual(
      values = c(
        "ACP" = 21,
        "BASS" = 24,
        "TRUG" = 22
      ),
      drop = FALSE
    ) +
    
    scale_x_discrete(
      labels = c(
        "BESE" = "BESE",
        "BARE" = "BARE"
      )
    ) +
    
    labs(
      title = title_text,
      x = NULL,
      y = y_label,
      fill = "Treatment",
      shape = "Site"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 13,
          hjust = 0.5
        ),
      
      axis.text.x =
        element_text(
          angle = 35,
          hjust = 1,
          colour = "black"
        ),
      
      axis.text.y =
        element_text(
          colour = "black"
        ),
      
      axis.title.y =
        element_text(
          colour = "black"
        ),
      
      legend.position =
        "bottom"
    )
  
  list(
    plot = panel,
    model = model,
    emmeans = emm_ci,
    pairwise = pairwise
  )
}


# -----------------------------
# 11. Run vegetation panels
# -----------------------------

veg_cover_panel <- make_veg_panel(
  data = veg_restore,
  response = "Mean_cover_pct",
  y_label = "Mean vegetation cover (%)",
  title_text = "Vegetation cover"
)

veg_height_panel <- make_veg_panel(
  data = veg_restore,
  response = "HEIGHT",
  y_label = "Maximum height",
  title_text = "Height"
)

veg_green_panel <- make_veg_panel(
  data = veg_restore,
  response = "GREEN_pct",
  y_label = "Green foliage (%)",
  title_text = "Green foliage"
)


# -----------------------------
# 12. Combine vegetation figure
# -----------------------------

vegetation_figure <- (
  veg_cover_panel$plot |
    veg_height_panel$plot |
    veg_green_panel$plot
) +
  
  patchwork::plot_layout(
    guides = "collect"
  ) +
  
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    legend.position = "bottom"
  )

vegetation_figure


# -----------------------------
# 13. Save vegetation figure
# -----------------------------

ggsave(
  filename = "M6_vegetation_BESE_BARE_three_panel.png",
  plot = vegetation_figure,
  width = 9,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "M6_vegetation_BESE_BARE_three_panel.pdf",
  plot = vegetation_figure,
  width = 9,
  height = 4,
  units = "in"
)


# -----------------------------
# 14. Save plant-condition table for microbiome integration
# -----------------------------

plant_condition <- veg_restore %>%
  
  transmute(
    SITE,
    TREATMENT,
    PLOT,
    PlotID,
    Mean_cover_pct,
    HEIGHT,
    GREEN_pct,
    BROWN_pct,
    RED_pct,
    Green_fraction = GREEN_pct / 100,
    Brown_fraction = BROWN_pct / 100,
    Height_flag
  )

write.csv(
  plant_condition,
  "M6_plant_condition_BESE_BARE_plot_level.csv",
  row.names = FALSE
)


# ============================================================
# M6 — MICROBIAL RESEMBLANCE TO NATURAL
#
# This section requires your distance-to-Natural tables.
#
# Required columns in each distance table:
# SITE, TREATMENT, PLOT, Distance_to_Natural
#
# Example expected objects:
# distance_prok_rz
# distance_fungi_rz
# ============================================================

standardise_distance_table <- function(distance_data) {
  
  names(distance_data) <- names(distance_data) %>%
    toupper() %>%
    gsub("[^A-Z0-9]+", "_", .) %>%
    gsub("^_|_$", "", .)
  
  site_distance_col <- find_col(
    names(distance_data),
    c("SITE")
  )
  
  treatment_distance_col <- find_col(
    names(distance_data),
    c("TREATMENT", "ECOLOGICAL_CONDITION", "CONDITION")
  )
  
  plot_distance_col <- find_col(
    names(distance_data),
    c("PLOT", "PLOT_NO", "NO", "BLOCK", "BLOQUE")
  )
  
  distance_col <- find_col(
    names(distance_data),
    c("DISTANCE_TO_NATURAL", "DISTANCE_NATURAL", "BRAY_DISTANCE_TO_NATURAL", "DISTANCE")
  )
  
  distance_data %>%
    
    transmute(
      SITE =
        toupper(
          trimws(
            as.character(.data[[site_distance_col]])
          )
        ),
      
      TREATMENT =
        toupper(
          trimws(
            as.character(.data[[treatment_distance_col]])
          )
        ),
      
      PLOT =
        as.integer(
          parse_number_safe(.data[[plot_distance_col]])
        ),
      
      Distance_to_Natural =
        as.numeric(.data[[distance_col]])
    ) %>%
    
    filter(
      TREATMENT %in% c("BESE", "BARE")
    ) %>%
    
    mutate(
      SITE =
        factor(
          SITE,
          levels = c("ACP", "BASS", "TRUG")
        ),
      
      TREATMENT =
        factor(
          TREATMENT,
          levels = c("BESE", "BARE")
        )
    )
}


fit_green_distance_model <- function(distance_data,
                                     kingdom_label) {
  
  distance_clean <- standardise_distance_table(
    distance_data
  )
  
  model_data <- distance_clean %>%
    
    left_join(
      plant_condition,
      by = c(
        "SITE",
        "TREATMENT",
        "PLOT"
      )
    ) %>%
    
    filter(
      !is.na(Distance_to_Natural),
      !is.na(GREEN_pct),
      !is.na(SITE),
      !is.na(TREATMENT)
    ) %>%
    
    mutate(
      Green_z =
        as.numeric(
          scale(GREEN_pct)
        )
    )
  
  model <- lm(
    Distance_to_Natural ~ SITE + TREATMENT + Green_z,
    data = model_data
  )
  
  model_test <- drop1(
    model,
    test = "F"
  )
  
  coef_table <- summary(model)$coefficients
  
  green_effect <- data.frame(
    Kingdom = kingdom_label,
    Estimate_per_1SD_green =
      coef_table["Green_z", "Estimate"],
    SE =
      coef_table["Green_z", "Std. Error"],
    t_value =
      coef_table["Green_z", "t value"],
    p_value =
      coef_table["Green_z", "Pr(>|t|)"],
    row.names = NULL
  )
  
  plot <- ggplot(
    model_data,
    aes(
      x = GREEN_pct,
      y = Distance_to_Natural
    )
  ) +
    
    geom_point(
      aes(
        fill = TREATMENT,
        shape = SITE
      ),
      size = 3,
      alpha = 0.85,
      colour = "black",
      stroke = 0.45
    ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      colour = "black",
      linewidth = 0.8
    ) +
    
    scale_fill_manual(
      values = c(
        "BESE" = "#FF5A5F",
        "BARE" = "#FFA630"
      ),
      drop = FALSE
    ) +
    
    scale_shape_manual(
      values = c(
        "ACP" = 21,
        "BASS" = 24,
        "TRUG" = 22
      ),
      drop = FALSE
    ) +
    
    labs(
      title = paste0(
        kingdom_label,
        " rhizosphere"
      ),
      x = "Green foliage (%)",
      y = "Distance to Natural reference",
      fill = "Treatment",
      shape = "Site"
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title =
        element_text(
          face = "bold",
          hjust = 0.5
        ),
      
      axis.text =
        element_text(
          colour = "black"
        ),
      
      legend.position =
        "bottom"
    )
  
  list(
    model_data = model_data,
    model = model,
    model_test = model_test,
    green_effect = green_effect,
    plot = plot
  )
}


# -----------------------------
# 15. Run distance-to-Natural models
# -----------------------------
# This will run only if these objects already exist in your R session.

if (
  exists("distance_prok_rz") &&
  exists("distance_fungi_rz")
) {
  
  prok_green_distance <- fit_green_distance_model(
    distance_data = distance_prok_rz,
    kingdom_label = "Prokaryotic"
  )
  
  fungi_green_distance <- fit_green_distance_model(
    distance_data = distance_fungi_rz,
    kingdom_label = "Fungal"
  )
  
  cat("\n========================================\n")
  cat("PROKARYOTIC DISTANCE MODEL\n")
  cat("========================================\n")
  
  print(
    prok_green_distance$model_test
  )
  
  print(
    summary(
      prok_green_distance$model
    )
  )
  
  cat("\n========================================\n")
  cat("FUNGAL DISTANCE MODEL\n")
  cat("========================================\n")
  
  print(
    fungi_green_distance$model_test
  )
  
  print(
    summary(
      fungi_green_distance$model
    )
  )
  
  green_distance_effects <- bind_rows(
    prok_green_distance$green_effect,
    fungi_green_distance$green_effect
  )
  
  write.csv(
    green_distance_effects,
    "M6_green_foliage_distance_to_Natural_effects.csv",
    row.names = FALSE
  )
  
  green_distance_figure <- (
    prok_green_distance$plot |
      fungi_green_distance$plot
  ) +
    
    patchwork::plot_layout(
      guides = "collect"
    ) +
    
    patchwork::plot_annotation(
      tag_levels = "A"
    ) &
    
    theme(
      legend.position = "bottom"
    )
  
  green_distance_figure
  
  ggsave(
    filename = "M6_green_foliage_distance_to_Natural.png",
    plot = green_distance_figure,
    width = 7.5,
    height = 4,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  ggsave(
    filename = "M6_green_foliage_distance_to_Natural.pdf",
    plot = green_distance_figure,
    width = 7.5,
    height = 4,
    units = "in"
  )
  
} else {
  
  message(
    "\nDistance-to-Natural objects were not found.\n",
    "The vegetation figure and plant-condition table were created.\n",
    "To run the microbial resemblance models, load or create these objects first:\n",
    "  distance_prok_rz\n",
    "  distance_fungi_rz\n",
    "Each must contain SITE, TREATMENT, PLOT, and Distance_to_Natural.\n"
  )
}

