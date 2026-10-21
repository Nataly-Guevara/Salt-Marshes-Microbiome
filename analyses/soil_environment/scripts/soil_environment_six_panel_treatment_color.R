# ============================================================
# SOIL ENVIRONMENTAL CONDITIONS
# SIX-PANEL FIGURE + SITE-ADJUSTED STATISTICS
#
# Model:
#   Environmental parameter ~ SITE + TREATMENT
#
# SITE      = fixed effect
# TREATMENT = fixed effect
#
# Figure:
#   raw observations
#   Site represented by shape
#   Treatment represented by colour
#   Site-adjusted treatment EMMs
#   95% CI
#   Tukey compact-letter groups
#
# Outputs:
#   figure PNG + PDF
#   EMM table
#   Tukey pairwise table
#   overall treatment-test table
#   complete RDS analysis object
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(readxl)
library(dplyr)
library(tibble)
library(ggplot2)
library(emmeans)
library(multcompView)
library(purrr)
library(patchwork)


# ============================================================
# 2. IMPORT DATA
# ============================================================

soil <- read_excel(
  "metadata_env_parameters.xlsx"
) %>%
  
  filter(
    !is.na(SITE),
    !is.na(TREATMENT)
  ) %>%
  
  mutate(
    
    SITE = factor(
      trimws(
        toupper(SITE)
      ),
      levels = c(
        "ACP",
        "BASS",
        "TRUG"
      )
    ),
    
    TREATMENT = factor(
      trimws(
        toupper(TREATMENT)
      ),
      levels = c(
        "NATURAL",
        "BESE",
        "BARE",
        "DEGRADED"
      )
    )
  )


# ============================================================
# 3. CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "SITE",
  "TREATMENT",
  "Moisture_avg",
  "Temperature_avg",
  "pH_avg",
  "NH4",
  "PO4",
  "Redox_avg"
)

missing_columns <- setdiff(
  required_columns,
  names(soil)
)

if (length(missing_columns) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


# ============================================================
# 4. VARIABLE SPECIFICATIONS
# ============================================================

variable_spec <- tribble(
  
  ~response,
  ~title,
  ~y_label,
  ~transformation,
  ~multiply_by,
  ~log_axis,
  ~zero_line,
  
  "Moisture_avg",
  "Moisture",
  "Soil moisture (%)",
  "none",
  100,
  FALSE,
  FALSE,
  
  "Temperature_avg",
  "Temperature",
  "Soil temperature (°C)",
  "none",
  1,
  FALSE,
  FALSE,
  
  "pH_avg",
  "pH",
  "Soil pH",
  "none",
  1,
  FALSE,
  FALSE,
  
  "NH4",
  "Ammonium",
  expression(NH[4]^"+"~"concentration"),
  "log",
  1,
  TRUE,
  FALSE,
  
  "PO4",
  "Phosphate",
  expression(PO[4]^"3-"~"concentration"),
  "log",
  1,
  TRUE,
  FALSE,
  
  "Redox_avg",
  "Redox",
  "Redox potential (mV)",
  "none",
  1,
  FALSE,
  TRUE
)


# ============================================================
# 5. FUNCTION FOR ONE ENVIRONMENTAL VARIABLE
# ============================================================

make_soil_panel <- function(
    data,
    response,
    title,
    y_label,
    transformation = "none",
    multiply_by = 1,
    log_axis = FALSE,
    zero_line = FALSE
) {
  
  
  # ----------------------------------------------------------
  # A. Colour palette
  # 
  # ----------------------------------------------------------
  
  treatment_colors <- c(
    "NATURAL" = "#9BC1BC",
    "BESE" = "#FF5A5F",
    "BARE" = "#FFA630",
    "DEGRADED" = "#5D576B"
  )
  
  
  # ----------------------------------------------------------
  # B. Site shapes
  # ----------------------------------------------------------
  
  site_shapes <- c(
    "ACP" = 16,
    "BASS" = 17,
    "TRUG" = 15
  )
  
  
  # ----------------------------------------------------------
  # C. Prepare data
  # ----------------------------------------------------------
  
  plot_data <- data %>%
    
    filter(
      !is.na(.data[[response]]),
      !is.na(SITE),
      !is.na(TREATMENT)
    ) %>%
    
    mutate(
      
      .raw_y = as.numeric(
        .data[[response]]
      ),
      
      .display_y =
        .raw_y *
        multiply_by
    )
  
  
  if (nrow(plot_data) == 0) {
    
    stop(
      response,
      ": no usable observations."
    )
  }
  
  
  # ----------------------------------------------------------
  # D. Statistical transformation
  # ----------------------------------------------------------
  
  if (transformation == "log") {
    
    if (
      any(
        plot_data$.raw_y <= 0,
        na.rm = TRUE
      )
    ) {
      
      stop(
        response,
        " contains zero or negative values and cannot be log-transformed."
      )
    }
    
    
    plot_data <- plot_data %>%
      
      mutate(
        .model_y =
          log(.raw_y)
      )
    
    
    back_transform <- function(x) {
      exp(x) *
        multiply_by
    }
    
    
  } else {
    
    
    plot_data <- plot_data %>%
      
      mutate(
        .model_y =
          .raw_y
      )
    
    
    back_transform <- function(x) {
      x *
        multiply_by
    }
  }
  
  
  # ----------------------------------------------------------
  # E. SITE-ADJUSTED FIXED-EFFECT MODEL
  # ----------------------------------------------------------
  
  model <- stats::lm(
    
    .model_y ~
      SITE +
      TREATMENT,
    
    data =
      plot_data
  )
  
  
  # ----------------------------------------------------------
  # F. Marginal overall treatment test
  #
  # drop1() compares:
  # SITE + TREATMENT
  # versus
  # SITE only
  #
  # Therefore this tests Treatment after accounting for Site.
  # ----------------------------------------------------------
  
  marginal_tests <- stats::drop1(
    model,
    test = "F"
  )
  
  marginal_tests <- as.data.frame(
    marginal_tests
  ) %>%
    rownames_to_column(
      "Term"
    )
  
  treatment_test <- marginal_tests %>%
    filter(
      Term == "TREATMENT"
    )
  
  
  # ----------------------------------------------------------
  # G. Site-adjusted treatment EMMs
  # ----------------------------------------------------------
  
  emm <- emmeans::emmeans(
    model,
    ~ TREATMENT
  )
  
  
  # ----------------------------------------------------------
  # H. 95% confidence intervals
  # ----------------------------------------------------------
  
  emm_ci_raw <- as.data.frame(
    
    confint(
      emm,
      level = 0.95
    )
  )
  
  
  emm_ci <- emm_ci_raw %>%
    
    transmute(
      
      TREATMENT,
      
      estimate =
        back_transform(
          emmean
        ),
      
      lower_ci =
        back_transform(
          lower.CL
        ),
      
      upper_ci =
        back_transform(
          upper.CL
        )
    )
  
  
  # ----------------------------------------------------------
  # I. Tukey-adjusted pairwise comparisons
  # ----------------------------------------------------------
  
  pairwise_table <- as.data.frame(
    
    emmeans::contrast(
      emm,
      method = "pairwise",
      adjust = "tukey"
    )
  )
  
  
  # ----------------------------------------------------------
  # J. Tukey compact-letter display
  # ----------------------------------------------------------
  
  p_values <-
    pairwise_table$p.value
  
  
  names(p_values) <- gsub(
    "\\s+",
    "",
    pairwise_table$contrast
  )
  
  
  letter_vector <- multcompView::multcompLetters(
    p_values,
    threshold = 0.05
  )$Letters
  
  
  letter_table <- tibble(
    
    TREATMENT = factor(
      names(
        letter_vector
      ),
      levels =
        levels(
          data$TREATMENT
        )
    ),
    
    tukey_group =
      unname(
        letter_vector
      )
  )
  
  
  plot_stats <- emm_ci %>%
    
    left_join(
      letter_table,
      by = "TREATMENT"
    )
  
  
  # ----------------------------------------------------------
  # K. Highest raw observation per Treatment
  # ----------------------------------------------------------
  
  raw_max <- plot_data %>%
    
    group_by(
      TREATMENT
    ) %>%
    
    summarise(
      
      raw_max = max(
        .display_y,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # L. Plot range
  # ----------------------------------------------------------
  
  display_range <- diff(
    
    range(
      
      c(
        plot_data$.display_y,
        plot_stats$lower_ci,
        plot_stats$upper_ci
      ),
      
      na.rm = TRUE
    )
  )
  
  
  if (
    !is.finite(display_range) ||
    display_range == 0
  ) {
    
    display_range <- 1
  }
  
  
  # ----------------------------------------------------------
  # M. Tukey-letter positions
  # Increase separation from EMM diamonds / confidence intervals
  # ----------------------------------------------------------
  
  plot_stats <- plot_stats %>%
    
    left_join(
      raw_max,
      by = "TREATMENT"
    ) %>%
    
    mutate(
      
      top_value = pmax(
        raw_max,
        upper_ci,
        na.rm = TRUE
      ),
      
      letter_y =
        if (log_axis) {
          
          # More separation on logarithmic axes
          top_value * 1.55
          
        } else {
          
          # More separation on linear axes
          top_value +
            0.10 * display_range
        }
    )
  
  # ----------------------------------------------------------
  # N. Raw point positioning
  # ----------------------------------------------------------
  
  raw_position <- position_jitterdodge(
    
    jitter.width =
      0.045,
    
    jitter.height =
      0,
    
    dodge.width =
      0.50,
    
    seed =
      123
  )
  
  
  # ----------------------------------------------------------
  # O. Base plot
  # ----------------------------------------------------------
  
  panel <- ggplot(
    
    plot_data,
    
    aes(
      x = TREATMENT,
      y = .display_y,
      colour = TREATMENT,
      shape = SITE
    )
  )
  
  
  # ----------------------------------------------------------
  # P. Zero line for Redox
  # ----------------------------------------------------------
  
  if (zero_line) {
    
    panel <- panel +
      
      geom_hline(
        
        yintercept = 0,
        
        linetype = "dashed",
        
        linewidth = 0.4,
        
        colour = "grey70"
      )
  }
  
  
  # ----------------------------------------------------------
  # Q. Raw observations + EMMs
  # ----------------------------------------------------------
  
  panel <- panel +
    
    geom_point(
      
      position =
        raw_position,
      
      size =
        2.6,
      
      alpha =
        0.80,
      
      stroke =
        0.4
    ) +
    
    
    geom_errorbar(
      
      data =
        plot_stats,
      
      aes(
        x = TREATMENT,
        ymin = lower_ci,
        ymax = upper_ci,
        colour = TREATMENT
      ),
      
      inherit.aes =
        FALSE,
      
      width =
        0.12,
      
      linewidth =
        0.8
    ) +
    
    
    geom_point(
      
      data =
        plot_stats,
      
      aes(
        x = TREATMENT,
        y = estimate,
        fill = TREATMENT
      ),
      
      inherit.aes =
        FALSE,
      
      shape =
        23,
      
      size =
        4,
      
      stroke =
        1,
      
      colour =
        "black"
    ) +
    
    
    geom_text(
      
      data =
        plot_stats,
      
      aes(
        x = TREATMENT,
        y = letter_y,
        label = tukey_group
      ),
      
      inherit.aes =
        FALSE,
      
      fontface =
        "bold",
      
      size =
        4
    ) +
    
    
    # --------------------------------------------------------
  # Treatment colours
  # --------------------------------------------------------
  
  scale_colour_manual(
    
    values =
      treatment_colors,
    
    guide =
      "none"
  ) +
    
    
    scale_fill_manual(
      
      values =
        treatment_colors,
      
      guide =
        "none"
    ) +
    
    
    # --------------------------------------------------------
  # Site shapes
  # --------------------------------------------------------
  
  scale_shape_manual(
    
    values =
      site_shapes,
    
    drop =
      FALSE
  ) +
    
    
    # --------------------------------------------------------
  # Treatment labels
  # --------------------------------------------------------
  
  scale_x_discrete(
    
    labels = c(
      
      "NATURAL" =
        "Natural",
      
      "BESE" =
        "BESE",
      
      "BARE" =
        "Bare",
      
      "DEGRADED" =
        "Degraded"
    )
  ) +
    
    
    labs(
      
      title =
        title,
      
      x =
        NULL,
      
      y =
        y_label,
      
      shape =
        "Site"
    ) +
    
    
    coord_cartesian(
      clip = "off"
    ) +
    
    
    theme_classic(
      base_size = 11
    ) +
    
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = 12
        ),
      
      axis.text.x =
        element_text(
          angle = 25,
          hjust = 1
        ),
      
      legend.position =
        "bottom"
    )
  
  
  # ----------------------------------------------------------
  # R. Axis
  # ----------------------------------------------------------
  
  if (log_axis) {
    
    panel <- panel +
      
      scale_y_log10(
        
        expand =
          expansion(
            mult = c(
              0.08,
              0.20
            )
          )
      )
    
    
  } else {
    
    
    panel <- panel +
      
      scale_y_continuous(
        
        expand =
          expansion(
            mult = c(
              0.06,
              0.18
            )
          )
      )
  }
  
  
  # ----------------------------------------------------------
  # S. Return
  # ----------------------------------------------------------
  
  list(
    
    plot =
      panel,
    
    model =
      model,
    
    treatment_test =
      treatment_test,
    
    emmeans_raw =
      emm_ci_raw,
    
    emmeans =
      plot_stats,
    
    pairwise =
      pairwise_table
  )
}


# ============================================================
# 6. RUN ALL SIX VARIABLES
# ============================================================

results <- purrr::pmap(
  
  variable_spec,
  
  function(
    response,
    title,
    y_label,
    transformation,
    multiply_by,
    log_axis,
    zero_line
  ) {
    
    make_soil_panel(
      
      data =
        soil,
      
      response =
        response,
      
      title =
        title,
      
      y_label =
        y_label,
      
      transformation =
        transformation,
      
      multiply_by =
        multiply_by,
      
      log_axis =
        log_axis,
      
      zero_line =
        zero_line
    )
  }
)


names(results) <-
  variable_spec$response


# ============================================================
# 7. EXTRACT PLOTS
# ============================================================
# ============================================================
# 7. EXTRACT PLOTS AND FORCE PANEL ASPECT RATIO
# ============================================================

plot_list <- purrr::map(
  variable_spec$response,
  function(x) {
    
    results[[x]]$plot +
      
      theme(
        # Makes each individual panel taller
        aspect.ratio = 0.72
      )
  }
)


# ============================================================
# 8. FORCE EXACT 2-COLUMN × 3-ROW LAYOUT
#
# A Moisture      | B Temperature
# C pH            | D Ammonium
# E Phosphate     | F Redox
# ============================================================

layout_design <- "
AB
CD
EF
"

final_figure <- patchwork::wrap_plots(
  plot_list,
  design = layout_design,
  guides = "collect"
) +
  
  patchwork::plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    legend.position = "bottom"
  )


# ============================================================
# 9. SAVE— PORTRAIT OUTPUT
# ============================================================

ggsave(
  filename =
    "soil_environment_six_panel_treatment_color_R.png",
  
  plot =
    final_figure,
  
  width =
    8,
  
  height =
    12,
  
  units =
    "in",
  
  dpi =
    600,
  
  bg =
    "white"
)


ggsave(
  filename =
    "soil_environment_six_panel_treatment_color_R.pdf",
  
  plot =
    final_figure,
  
  width =
    8,
  
  height =
    12,
  
  units =
    "in"
)


# ============================================================
# 10. DISPLAY
# ============================================================

final_figure




 

# ============================================================
# 11. COMBINE EMM RESULTS
# ============================================================

emm_results <- purrr::imap_dfr(
  
  results,
  
  function(
    result,
    variable_name
  ) {
    
    result$emmeans %>%
      
      mutate(
        
        Variable =
          variable_name,
        
        .before =
          1
      )
  }
)


# ============================================================
# 12. COMBINE TUKEY PAIRWISE RESULTS
# ============================================================

pairwise_results <- purrr::imap_dfr(
  
  results,
  
  function(
    result,
    variable_name
  ) {
    
    result$pairwise %>%
      
      mutate(
        
        Variable =
          variable_name,
        
        .before =
          1
      )
  }
)


# ============================================================
# 13. COMBINE OVERALL TREATMENT TESTS
# ============================================================

overall_treatment_tests <- purrr::imap_dfr(
  
  results,
  
  function(
    result,
    variable_name
  ) {
    
    result$treatment_test %>%
      
      mutate(
        
        Variable =
          variable_name,
        
        .before =
          1
      )
  }
)


# ============================================================
# 14. SAVE TABLES
# ============================================================

write.csv(
  
  emm_results,
  
  "soil_environment_six_panel_treatment_color_EMMs.csv",
  
  row.names =
    FALSE
)


write.csv(
  
  pairwise_results,
  
  "soil_environment_six_panel_treatment_color_Tukey.csv",
  
  row.names =
    FALSE
)


write.csv(
  
  overall_treatment_tests,
  
  "soil_environment_six_panel_overall_treatment_tests.csv",
  
  row.names =
    FALSE
)


# ============================================================
# 15. PRINT OVERALL TREATMENT TESTS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "SITE-ADJUSTED OVERALL TREATMENT TESTS\n"
)

cat(
  "========================================\n"
)

print(
  overall_treatment_tests,
  row.names = FALSE
)


# ============================================================
# 16. PRINT TUKEY PAIRWISE RESULTS
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "SITE-ADJUSTED TUKEY PAIRWISE RESULTS\n"
)

cat(
  "========================================\n"
)

print(
  pairwise_results,
  row.names = FALSE
)


# ============================================================
# 17. SAVE COMPLETE ANALYSIS
# ============================================================

saveRDS(
  
  results,
  
  file =
    "soil_environment_six_panel_complete_analysis.rds"
)


cat(
  "\n========================================\n"
)

cat(
  "SOIL ENVIRONMENTAL ANALYSIS COMPLETE\n"
)

cat(
  "========================================\n"
)
