# ============================================================
# FINAL HORIZONTAL STACKED BARPLOT
# Site × Treatment microbial composition
#
# Features:
#   - Horizontal stacked bars
#   - Separate Prokaryotic and Fungal panels
#   - Separate legends
#   - Phyla ordered logically
#   - Classes ordered within each phylum by abundance
#   - Taxon names only in legends
#
# Requires:
#   plot_data
#
# plot_data must contain:
#   Community, Site, Treatment, Phylum_plot, Class_plot,
#   Segment_ID, Mean_RA_percent
# ============================================================

library(dplyr)
library(ggplot2)
library(readr)
library(grid)
library(patchwork)

dir.create(
  "microbial_overview_outputs",
  showWarnings = FALSE
)

# ============================================================
# 1. Load plot_data if not already in memory
# ============================================================

if (!exists("plot_data")) {
  plot_data <- read_csv(
    "microbial_overview_outputs/Figure_X_overall_composition_site_treatment_plot_data.csv",
    show_col_types = FALSE
  )
}

# ============================================================
# 2. Palette
# Phylum = colour family
# Class = shade within that family
# ============================================================

class_cols <- c(
  
  # ==========================================================
  # PROKARYOTIC COMMUNITY
  # ==========================================================
  
  # Actinobacteriota — blue slate family
  "Prokaryotic__Actinobacteriota__Acidimicrobiia" = "#2F5F6B",
  "Prokaryotic__Actinobacteriota__Actinobacteria" = "#3E7281",
  "Prokaryotic__Actinobacteriota__Thermoleophilia" = "#619FB0",
  "Prokaryotic__Actinobacteriota__Other classes" = "#A9CFD8",
  
  # Proteobacteria — muted teal family
  "Prokaryotic__Proteobacteria__Alphaproteobacteria" = "#4F8C75",
  "Prokaryotic__Proteobacteria__Gammaproteobacteria" = "#81B399",
  "Prokaryotic__Proteobacteria__Other classes" = "#C7DDD1",
  
  # Halobacterota — deep blue family
  "Prokaryotic__Halobacterota__Halobacteria" = "#264653",
  "Prokaryotic__Halobacterota__Other classes" = "#93A8B0",
  
  # Chloroflexi — pacific blue family
  "Prokaryotic__Chloroflexi__Anaerolineae" = "#356B7A",
  "Prokaryotic__Chloroflexi__Chloroflexia" = "#619FB0",
  "Prokaryotic__Chloroflexi__KD4-96" = "#8AB9C5",
  "Prokaryotic__Chloroflexi__Dehalococcoidia" = "#B9D7DE",
  "Prokaryotic__Chloroflexi__Other classes" = "#DCECEF",
  
  # Bacteroidota — coral family
  "Prokaryotic__Bacteroidota__Bacteroidia" = "#D95863",
  "Prokaryotic__Bacteroidota__Rhodothermia" = "#EC6973",
  "Prokaryotic__Bacteroidota__Other classes" = "#F7B3B8",
  
  # Gemmatimonadota — salmon family
  "Prokaryotic__Gemmatimonadota__PAUC43f marine benthic group" = "#E6786D",
  "Prokaryotic__Gemmatimonadota__BD2-11 terrestrial group" = "#FE887C",
  "Prokaryotic__Gemmatimonadota__Other classes" = "#FFC2BA",
  
  # Desulfobacterota — ash green family
  "Prokaryotic__Desulfobacterota__Desulfuromonadia" = "#6F7F6A",
  "Prokaryotic__Desulfobacterota__Desulfobulbia" = "#ACBDA3",
  "Prokaryotic__Desulfobacterota__Desulfobacteria" = "#C8D4C1",
  "Prokaryotic__Desulfobacterota__Other classes" = "#E5EAE1",
  
  # Firmicutes — graphite family
  "Prokaryotic__Firmicutes__Bacilli" = "#4D5656",
  "Prokaryotic__Firmicutes__Clostridia" = "#7B8584",
  "Prokaryotic__Firmicutes__Other classes" = "#C7CECC",
  
  # Acidobacteriota — muted blue-grey family
  "Prokaryotic__Acidobacteriota__Thermoanaerobaculia" = "#5B7286",
  "Prokaryotic__Acidobacteriota__Other classes" = "#B6C6CF",
  
  # Halanaerobiaeota — dark blue slate
  "Prokaryotic__Halanaerobiaeota__Halanaerobiia" = "#223D47",
  
  # Other prokaryotic phyla — bone/neutral
  "Prokaryotic__Other phyla__Other phyla" = "#E1DACE",
  
  
  # ==========================================================
  # FUNGAL COMMUNITY
  # ==========================================================
  
  # Ascomycota — salmon-pink family
  "Fungal__Ascomycota__Dothideomycetes" = "#C84E59",
  "Fungal__Ascomycota__Sordariomycetes" = "#FF928B",
  "Fungal__Ascomycota__Unclassified class" = "#FBA8A0",
  "Fungal__Ascomycota__Eurotiomycetes" = "#FEC3A6",
  "Fungal__Ascomycota__Leotiomycetes" = "#FFD7C2",
  "Fungal__Ascomycota__Other classes" = "#FBE5DA",
  
  # Basidiomycota — orange-peach family
  "Fungal__Basidiomycota__Unclassified class" = "#D98256",
  "Fungal__Basidiomycota__Agaricomycetes" = "#FFAC81",
  "Fungal__Basidiomycota__Cystobasidiomycetes" = "#FFC19F",
  "Fungal__Basidiomycota__Tremellomycetes" = "#FEC3A6",
  "Fungal__Basidiomycota__Agaricostilbomycetes" = "#FFE0CC",
  "Fungal__Basidiomycota__Other classes" = "#FFF0E6",
  
  # Unclassified fungi — tea green
  "Fungal__Unclassified fungi__Unclassified class" = "#CDEAC0",
  
  # Other fungal phyla — vanilla
  "Fungal__Other phyla__Other phyla" = "#EFE9AE"
)

# ============================================================
# 3. Site, treatment and phylum orders
# Correct order:
#   ACP Natural
#   ACP BESE
#   ACP BARE
#   ACP Degraded
#   BASS Natural
#   BASS BESE
#   BASS BARE
#   BASS Degraded
#   TRUG Natural
#   TRUG BESE
#   TRUG BARE
#   TRUG Degraded
# ============================================================

site_order <- c("ACP", "BASS", "TRUG")

treatment_order <- c(
  "Natural",
  "BESE",
  "BARE",
  "Degraded"
)

group_order <- c(
  "ACP — Natural",
  "ACP — BESE",
  "ACP — BARE",
  "ACP — Degraded",
  "BASS — Natural",
  "BASS — BESE",
  "BASS — BARE",
  "BASS — Degraded",
  "TRUG — Natural",
  "TRUG — BESE",
  "TRUG — BARE",
  "TRUG — Degraded"
)

# ggplot puts the first factor level at the bottom in horizontal plots.
# Reversing makes ACP appear at the top.
group_order_plot <- rev(group_order)

prok_phylum_order <- c(
  "Actinobacteriota",
  "Proteobacteria",
  "Halobacterota",
  "Chloroflexi",
  "Bacteroidota",
  "Gemmatimonadota",
  "Desulfobacterota",
  "Firmicutes",
  "Acidobacteriota",
  "Halanaerobiaeota",
  "Other phyla"
)

fungal_phylum_order <- c(
  "Ascomycota",
  "Basidiomycota",
  "Unclassified fungi",
  "Other phyla"
)


# ============================================================
# 4. Function to make one horizontal stacked plot
# ============================================================

make_horizontal_taxa_plot <- function(dat,
                                      community_name,
                                      plot_title,
                                      phylum_order,
                                      legend_ncol = 1) {
  
  dat_use <- dat %>%
    filter(Community == community_name) %>%
    mutate(
      Site = factor(Site, levels = site_order),
      Treatment = factor(Treatment, levels = treatment_order),
      Group_label = paste(Site, Treatment, sep = " — "),
      Group_label = factor(
        Group_label,
        levels = group_order_plot
      ),
      y_pos = as.numeric(Group_label),
      Segment_ID = as.character(Segment_ID),
      Phylum_rank = match(Phylum_plot, phylum_order)
    )
  
  # ----------------------------------------------------------
  # Segment order:
  #   phylum first
  #   class within phylum second
  # ----------------------------------------------------------
  
  segment_order <- dat_use %>%
    group_by(
      Phylum_plot,
      Class_plot,
      Segment_ID
    ) %>%
    summarise(
      Overall_mean_percent = mean(Mean_RA_percent, na.rm = TRUE),
      Phylum_rank = first(Phylum_rank),
      .groups = "drop"
    ) %>%
    arrange(
      Phylum_rank,
      desc(Overall_mean_percent)
    ) %>%
    group_by(Phylum_plot) %>%
    mutate(
      Class_rank_within_phylum = row_number()
    ) %>%
    ungroup() %>%
    arrange(
      Phylum_rank,
      Class_rank_within_phylum
    ) %>%
    mutate(
      Segment_rank = row_number(),
      Legend_label = case_when(
        Phylum_plot == Class_plot ~ Phylum_plot,
        TRUE ~ paste0(Phylum_plot, " — ", Class_plot)
      )
    )
  
  legend_breaks <- segment_order$Segment_ID
  
  legend_labels <- segment_order$Legend_label
  names(legend_labels) <- legend_breaks
  
  missing_cols <- setdiff(
    legend_breaks,
    names(class_cols)
  )
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing colours for:\n",
        paste(missing_cols, collapse = "\n")
      )
    )
  }
  
  # ----------------------------------------------------------
  # Manual stacking
  # Left to right = phylum order, then class order.
  # ----------------------------------------------------------
  
  dat_ordered <- dat_use %>%
    left_join(
      segment_order %>%
        select(
          Segment_ID,
          Phylum_rank,
          Class_rank_within_phylum,
          Segment_rank
        ),
      by = "Segment_ID"
    ) %>%
    mutate(
      Segment_ID = factor(
        Segment_ID,
        levels = legend_breaks
      )
    ) %>%
    arrange(
      Group_label,
      Segment_rank
    ) %>%
    group_by(Group_label) %>%
    mutate(
      xmin = lag(cumsum(Mean_RA_percent), default = 0),
      xmax = cumsum(Mean_RA_percent),
      ymin = y_pos - 0.38,
      ymax = y_pos + 0.38
    ) %>%
    ungroup()
  
  # ----------------------------------------------------------
  # Check totals
  # ----------------------------------------------------------
  
  bar_check <- dat_ordered %>%
    group_by(
      Community,
      Site,
      Treatment,
      Group_label
    ) %>%
    summarise(
      total_percent = sum(Mean_RA_percent),
      n_segments = n(),
      .groups = "drop"
    ) %>%
    arrange(Group_label)
  
  cat("\n========================================\n")
  cat("BAR TOTAL CHECK:", community_name, "\n")
  cat("========================================\n")
  print(bar_check, n = Inf, width = Inf)
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot() +
    geom_rect(
      data = dat_ordered,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = Segment_ID
      ),
      colour = "white",
      linewidth = 0.18
    ) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 25),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = seq_along(group_order_plot),
      labels = group_order_plot,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_fill_manual(
      values = class_cols[legend_breaks],
      breaks = legend_breaks,
      labels = legend_labels,
      name = paste0(community_name, "\nphylum — class")
    ) +
    labs(
      title = plot_title,
      x = "Mean relative abundance",
      y = NULL
    ) +
    guides(
      fill = guide_legend(
        ncol = legend_ncol,
        byrow = FALSE,
        keyheight = unit(0.32, "cm"),
        keywidth = unit(0.42, "cm")
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = "grey86",
        linewidth = 0.25
      ),
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      axis.text.y = element_text(
        colour = "grey20",
        size = 8.8
      ),
      axis.text.x = element_text(
        colour = "grey25",
        size = 9
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 10
      ),
      legend.position = "right",
      legend.title = element_text(
        face = "bold",
        size = 8.5
      ),
      legend.text = element_text(
        size = 6.7
      ),
      legend.key = element_rect(
        fill = "white",
        colour = NA
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
        t = 6,
        r = 6,
        b = 6,
        l = 6
      )
    )
  
  list(
    plot = p,
    plot_data = dat_ordered,
    legend_table = segment_order,
    bar_check = bar_check
  )
}

# ============================================================
# 5. Make plots
# ============================================================

prok_out <- make_horizontal_taxa_plot(
  dat = plot_data,
  community_name = "Prokaryotic",
  plot_title = "A. Prokaryotic community",
  phylum_order = prok_phylum_order,
  legend_ncol = 2
)

fungal_out <- make_horizontal_taxa_plot(
  dat = plot_data,
  community_name = "Fungal",
  plot_title = "B. Fungal community",
  phylum_order = fungal_phylum_order,
  legend_ncol = 1
)

# ============================================================
# 6. Combine plots, keeping separate legends
# ============================================================

overall_horizontal_plot <-
  prok_out$plot /
  fungal_out$plot +
  patchwork::plot_layout(
    heights = c(1.25, 1),
    guides = "keep"
  ) +
  patchwork::plot_annotation(
    title = "Overall microbial composition by site and ecological condition",
    subtitle = "Class-level groups nested within dominant phyla; bulk-soil and rhizosphere samples combined",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        size = 14
      ),
      plot.subtitle = element_text(
        hjust = 0,
        size = 10
      ),
      plot.background = element_rect(
        fill = "white",
        colour = NA
      )
    )
  )

overall_horizontal_plot

# ============================================================
# 7. Save figure and data
# ============================================================

ggsave(
  filename = "microbial_overview_outputs/Figure_X_overall_composition_HORIZONTAL_PHYLA_ORDERED.png",
  plot = overall_horizontal_plot,
  width = 14,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = "microbial_overview_outputs/Figure_X_overall_composition_HORIZONTAL_PHYLA_ORDERED.pdf",
  plot = overall_horizontal_plot,
  width = 14,
  height = 10,
  units = "in",
  bg = "white"
)

write_csv(
  prok_out$plot_data,
  "microbial_overview_outputs/Figure_X_prokaryotic_HORIZONTAL_PHYLA_ORDERED_plot_data.csv"
)

write_csv(
  fungal_out$plot_data,
  "microbial_overview_outputs/Figure_X_fungal_HORIZONTAL_PHYLA_ORDERED_plot_data.csv"
)

write_csv(
  prok_out$legend_table,
  "microbial_overview_outputs/Figure_X_prokaryotic_HORIZONTAL_PHYLA_ORDERED_legend_table.csv"
)

write_csv(
  fungal_out$legend_table,
  "microbial_overview_outputs/Figure_X_fungal_HORIZONTAL_PHYLA_ORDERED_legend_table.csv"
)

write_csv(
  prok_out$bar_check,
  "microbial_overview_outputs/Figure_X_prokaryotic_HORIZONTAL_PHYLA_ORDERED_bar_check.csv"
)

write_csv(
  fungal_out$bar_check,
  "microbial_overview_outputs/Figure_X_fungal_HORIZONTAL_PHYLA_ORDERED_bar_check.csv"
)

cat("\nSaved final phylum-ordered horizontal figure:\n")
cat("microbial_overview_outputs/Figure_X_overall_composition_HORIZONTAL_PHYLA_ORDERED.png\n")
cat("microbial_overview_outputs/Figure_X_overall_composition_HORIZONTAL_PHYLA_ORDERED.pdf\n")
