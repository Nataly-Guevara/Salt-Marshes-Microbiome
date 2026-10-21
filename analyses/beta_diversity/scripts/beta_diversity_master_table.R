# ============================================================
# FINAL BETA-DIVERSITY MASTER TABLE
# H1 + H2 + H3
#
# H1 = Bulk soil:
#      Natural, BESE, BARE, Degraded
#
# H2 = Rhizosphere:
#      Natural, BESE, BARE
#
# H3 = Paired bulk soil vs rhizosphere:
#      Natural, BESE, BARE
#
# Outputs:
# 1. Complete PERMANOVA master table
# 2. Complete PERMDISP master table
# 3. Combined beta-diversity master table
# 4. Main-text statistics table
#
# IMPORTANT:
# Site R2 is retained descriptively.
# Site p-values are NOT interpreted when permutations
# are restricted within Site.
#
# For H3, Site and ecological condition are between-plot
# variables. With permutations restricted within PlotID,
# inferential p-values are used for Soil habitat and the
# Condition × Habitat interaction.
# ============================================================

suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(dplyr)
  library(tibble)
})

set.seed(123)


# ============================================================
# 0. HELPER FUNCTIONS
# ============================================================

get_first_object <- function(candidates) {
  
  found <- candidates[
    vapply(
      candidates,
      exists,
      logical(1),
      inherits = TRUE
    )
  ]
  
  if (length(found) == 0) {
    stop(
      "Could not find any of these objects:\n",
      paste(candidates, collapse = ", ")
    )
  }
  
  get(
    found[1],
    inherits = TRUE
  )
}


find_column_beta <- function(
    dat,
    candidates,
    required = TRUE
) {
  
  found <- candidates[
    candidates %in% names(dat)
  ]
  
  if (length(found) == 0) {
    
    if (required) {
      stop(
        "Could not find column. Tried: ",
        paste(candidates, collapse = ", ")
      )
    }
    
    return(NA_character_)
  }
  
  found[1]
}


clean_condition <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  x[
    x %in% c(
      "NAT",
      "NATURAL"
    )
  ] <- "Natural"
  
  x[
    x == "BESE"
  ] <- "BESE"
  
  x[
    x == "BARE"
  ] <- "BARE"
  
  x[
    x %in% c(
      "DEG",
      "DEGRADED",
      "DEGREST"
    )
  ] <- "Degraded"
  
  x
}


clean_habitat <- function(x) {
  
  x <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  x[
    x %in% c(
      "BS",
      "BULK",
      "BULK SOIL",
      "BULK_SOIL"
    )
  ] <- "Bulk soil"
  
  x[
    x %in% c(
      "RZ",
      "RHIZOSPHERE",
      "RHIZOSPHERE SOIL",
      "RHIZOSPHERE_SOIL"
    )
  ] <- "Rhizosphere"
  
  x
}


align_metadata <- function(
    bray,
    meta
) {
  
  meta <- as.data.frame(
    meta,
    stringsAsFactors = FALSE
  )
  
  ids <- attr(
    bray,
    "Labels"
  )
  
  if (is.null(ids)) {
    stop(
      "Bray-Curtis object has no sample labels."
    )
  }
  
  if (
    !all(ids %in% rownames(meta))
  ) {
    
    sample_col <- find_column_beta(
      meta,
      c(
        "SampleID",
        "Sample_ID",
        "sample",
        "Sample"
      ),
      required = FALSE
    )
    
    if (!is.na(sample_col)) {
      rownames(meta) <-
        as.character(
          meta[[sample_col]]
        )
    }
  }
  
  if (
    !all(ids %in% rownames(meta))
  ) {
    stop(
      "Metadata cannot be matched to Bray-Curtis sample IDs."
    )
  }
  
  meta <- meta[
    ids,
    ,
    drop = FALSE
  ]
  
  stopifnot(
    identical(
      rownames(meta),
      ids
    )
  )
  
  meta
}


standardise_H12_meta <- function(
    bray,
    meta
) {
  
  meta <- align_metadata(
    bray,
    meta
  )
  
  site_col <- find_column_beta(
    meta,
    c(
      "Site",
      "SITE",
      "site"
    )
  )
  
  condition_col <- find_column_beta(
    meta,
    c(
      "Condition",
      "Ecological_condition",
      "EcologicalCondition",
      "Treatment"
    )
  )
  
  meta$Site_master <-
    factor(
      as.character(
        meta[[site_col]]
      )
    )
  
  meta$Condition_master <-
    factor(
      clean_condition(
        meta[[condition_col]]
      )
    )
  
  meta
}


standardise_H3_meta <- function(
    bray,
    meta
) {
  
  meta <- align_metadata(
    bray,
    meta
  )
  
  site_col <- find_column_beta(
    meta,
    c(
      "Site",
      "SITE",
      "site"
    )
  )
  
  condition_col <- find_column_beta(
    meta,
    c(
      "Ecological_condition",
      "Condition",
      "Treatment"
    )
  )
  
  habitat_col <- find_column_beta(
    meta,
    c(
      "Soil_habitat",
      "SoilHabitat",
      "SoilPortion",
      "Habitat"
    )
  )
  
  plot_col <- find_column_beta(
    meta,
    c(
      "PlotID",
      "Plot_ID",
      "Plot"
    )
  )
  
  meta$Site_master <-
    factor(
      as.character(
        meta[[site_col]]
      )
    )
  
  meta$Condition_master <-
    factor(
      clean_condition(
        meta[[condition_col]]
      ),
      levels = c(
        "Natural",
        "BESE",
        "BARE"
      )
    )
  
  meta$Habitat_master <-
    factor(
      clean_habitat(
        meta[[habitat_col]]
      ),
      levels = c(
        "Bulk soil",
        "Rhizosphere"
      )
    )
  
  meta$Plot_master <-
    factor(
      as.character(
        meta[[plot_col]]
      )
    )
  
  meta
}


make_site_permutations <- function(
    meta,
    nperm = 9999
) {
  
  control <- permute::how(
    blocks = meta$Site_master,
    nperm = nperm
  )
  
  permute::shuffleSet(
    n = nrow(meta),
    nset = nperm,
    control = control
  )
}


make_plot_permutations <- function(
    meta,
    nperm = 9999
) {
  
  control <- permute::how(
    blocks = meta$Plot_master,
    nperm = nperm
  )
  
  permute::shuffleSet(
    n = nrow(meta),
    nset = nperm,
    control = control
  )
}


extract_row <- function(
    fit,
    term
) {
  
  tab <- as.data.frame(
    fit
  )
  
  if (
    !term %in% rownames(tab)
  ) {
    return(
      tibble(
        df = NA_real_,
        pseudo_F = NA_real_,
        R2 = NA_real_,
        p_raw = NA_real_
      )
    )
  }
  
  tibble(
    df =
      tab[
        term,
        "Df"
      ],
    
    pseudo_F =
      tab[
        term,
        "F"
      ],
    
    R2 =
      tab[
        term,
        "R2"
      ],
    
    p_raw =
      tab[
        term,
        "Pr(>F)"
      ]
  )
}


# ============================================================
# 1. LOCATE FINAL H1 OBJECTS
# ============================================================

bray_H1_prok <- get_first_object(
  c(
    "bray_bac_bulk",
    "bray_prok_bulk",
    "bray_bulk_prok"
  )
)

meta_H1_prok <- get_first_object(
  c(
    "meta_bac_bulk_df",
    "meta_bac_bulk",
    "meta_prok_bulk"
  )
)

bray_H1_fun <- get_first_object(
  c(
    "bray_fun_bulk",
    "bray_fungi_bulk",
    "bray_bulk_fun"
  )
)

meta_H1_fun <- get_first_object(
  c(
    "meta_fun_bulk_df",
    "meta_fun_bulk",
    "meta_fungi_bulk"
  )
)


# ============================================================
# 2. LOCATE FINAL H2 OBJECTS
# ============================================================

if (
  exists("bac") &&
  is.list(bac) &&
  all(
    c(
      "bray",
      "meta"
    ) %in% names(bac)
  )
) {
  
  bray_H2_prok <- bac$bray
  meta_H2_prok <- bac$meta
  
} else {
  
  bray_H2_prok <- get_first_object(
    c(
      "bray_bac_rz",
      "bray_prok_rz",
      "bray_rz_prok"
    )
  )
  
  meta_H2_prok <- get_first_object(
    c(
      "meta_bac_rz_df",
      "meta_bac_rz",
      "meta_prok_rz"
    )
  )
}


if (
  exists("fun") &&
  is.list(fun) &&
  all(
    c(
      "bray",
      "meta"
    ) %in% names(fun)
  )
) {
  
  bray_H2_fun <- fun$bray
  meta_H2_fun <- fun$meta
  
} else {
  
  bray_H2_fun <- get_first_object(
    c(
      "bray_fun_rz",
      "bray_fungi_rz",
      "bray_rz_fun"
    )
  )
  
  meta_H2_fun <- get_first_object(
    c(
      "meta_fun_rz_df",
      "meta_fun_rz",
      "meta_fungi_rz"
    )
  )
}


# ============================================================
# 3. LOCATE FINAL H3 OBJECTS
# ============================================================

stopifnot(
  exists("bray_C_prok"),
  exists("meta_C"),
  exists("bray_C_fun"),
  exists("meta_C_fun")
)

bray_H3_prok <- bray_C_prok
meta_H3_prok <- meta_C

bray_H3_fun <- bray_C_fun
meta_H3_fun <- meta_C_fun


# ============================================================
# 4. H1 / H2 OMNIBUS PERMANOVA
#
# Site + ecological condition
#
# Condition tested with permutations restricted within Site.
#
# Site:
# R2 and pseudo-F reported descriptively.
# p-value intentionally set to NA.
# ============================================================

run_H12_omnibus <- function(
    bray,
    meta,
    hypothesis,
    dataset,
    kingdom
) {
  
  meta <- standardise_H12_meta(
    bray,
    meta
  )
  
  perms <- make_site_permutations(
    meta
  )
  
  set.seed(123)
  
  fit <- vegan::adonis2(
    bray ~
      Site_master +
      Condition_master,
    data = meta,
    permutations = perms,
    by = "terms"
  )
  
  site_row <- extract_row(
    fit,
    "Site_master"
  )
  
  condition_row <- extract_row(
    fit,
    "Condition_master"
  )
  
  # Whole model R2
  whole_R2 <-
    site_row$R2 +
    condition_row$R2
  
  # Adjusted whole-model R2
  db <- vegan::dbrda(
    bray ~
      Site_master +
      Condition_master,
    data = meta
  )
  
  whole_adj_R2 <-
    vegan::RsquareAdj(
      db
    )$adj.r.squared
  
  rows <- bind_rows(
    
    tibble(
      Hypothesis = hypothesis,
      Dataset = dataset,
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + ecological condition",
      Term = "Site",
      Comparison = NA_character_,
      n = nrow(meta),
      df = site_row$df,
      pseudo_F = site_row$pseudo_F,
      R2 = site_row$R2,
      Adjusted_R2 = NA_real_,
      p_raw = NA_real_,
      p_adjusted = NA_real_,
      Adjustment = NA_character_,
      Inferential_status =
        "Descriptive contribution; Site used as permutation block",
      Permutation_structure =
        "9999 permutations restricted within Site"
    ),
    
    tibble(
      Hypothesis = hypothesis,
      Dataset = dataset,
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + ecological condition",
      Term = "Ecological condition",
      Comparison = NA_character_,
      n = nrow(meta),
      df = condition_row$df,
      pseudo_F =
        condition_row$pseudo_F,
      R2 =
        condition_row$R2,
      Adjusted_R2 = NA_real_,
      p_raw =
        condition_row$p_raw,
      p_adjusted = NA_real_,
      Adjustment = "None",
      Inferential_status =
        "Primary inferential test",
      Permutation_structure =
        "9999 permutations restricted within Site"
    ),
    
    tibble(
      Hypothesis = hypothesis,
      Dataset = dataset,
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + ecological condition",
      Term = "Whole model",
      Comparison = NA_character_,
      n = nrow(meta),
      df =
        site_row$df +
        condition_row$df,
      pseudo_F = NA_real_,
      R2 = whole_R2,
      Adjusted_R2 =
        whole_adj_R2,
      p_raw = NA_real_,
      p_adjusted = NA_real_,
      Adjustment = NA_character_,
      Inferential_status =
        "Whole-model variance explained",
      Permutation_structure =
        "R2 from complete Site + ecological condition model"
    )
  )
  
  list(
    table = rows,
    fit = fit,
    meta = meta
  )
}


# ============================================================
# 5. H1 / H2 PAIRWISE PERMANOVA + PERMDISP
# ============================================================

run_H12_pairwise <- function(
    bray,
    meta,
    hypothesis,
    dataset,
    kingdom,
    pairs
) {
  
  meta <- standardise_H12_meta(
    bray,
    meta
  )
  
  bray_mat <- as.matrix(
    bray
  )
  
  perma_rows <- list()
  disp_rows <- list()
  
  for (
    i in seq_along(pairs)
  ) {
    
    pair <- pairs[[i]]
    
    keep <-
      meta$Condition_master %in%
      pair
    
    meta_sub <-
      droplevels(
        meta[
          keep,
          ,
          drop = FALSE
        ]
      )
    
    ids <-
      rownames(
        meta_sub
      )
    
    bray_sub <-
      as.dist(
        bray_mat[
          ids,
          ids
        ]
      )
    
    perms <-
      make_site_permutations(
        meta_sub
      )
    
    # -----------------------------
    # PERMANOVA
    # -----------------------------
    
    set.seed(123)
    
    ad <- vegan::adonis2(
      bray_sub ~
        Site_master +
        Condition_master,
      data = meta_sub,
      permutations = perms,
      by = "terms"
    )
    
    cond_row <-
      extract_row(
        ad,
        "Condition_master"
      )
    
    contrast_name <-
      paste(
        pair,
        collapse = " vs "
      )
    
    perma_rows[[i]] <- tibble(
      Hypothesis = hypothesis,
      Dataset = dataset,
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Pairwise ecological condition after Site",
      Term =
        "Ecological condition",
      Comparison =
        contrast_name,
      n =
        nrow(meta_sub),
      df =
        cond_row$df,
      pseudo_F =
        cond_row$pseudo_F,
      R2 =
        cond_row$R2,
      Adjusted_R2 =
        NA_real_,
      p_raw =
        cond_row$p_raw,
      p_adjusted =
        NA_real_,
      Adjustment =
        "Holm",
      Inferential_status =
        "Planned pairwise comparison",
      Permutation_structure =
        "9999 permutations restricted within Site"
    )
    
    # -----------------------------
    # PERMDISP
    # -----------------------------
    
    bd <- vegan::betadisper(
      bray_sub,
      group =
        meta_sub$Condition_master,
      type = "median",
      bias.adjust = TRUE
    )
    
    set.seed(123)
    
    pd <- vegan::permutest(
      bd,
      permutations = perms
    )
    
    pd_tab <-
      as.data.frame(
        pd$tab
      )
    
    disp_rows[[i]] <- tibble(
      Hypothesis = hypothesis,
      Dataset = dataset,
      Kingdom = kingdom,
      Test_type = "PERMDISP",
      Model =
        "Pairwise dispersion diagnostic",
      Term =
        "Ecological condition",
      Comparison =
        contrast_name,
      n =
        nrow(meta_sub),
      df =
        pd_tab$Df[1],
      pseudo_F =
        pd_tab$F[1],
      R2 =
        NA_real_,
      Adjusted_R2 =
        NA_real_,
      p_raw =
        pd_tab$`Pr(>F)`[1],
      p_adjusted =
        NA_real_,
      Adjustment =
        "Holm",
      Inferential_status =
        "Dispersion diagnostic",
      Permutation_structure =
        "9999 permutations restricted within Site"
    )
  }
  
  perma_tab <-
    bind_rows(
      perma_rows
    ) %>%
    mutate(
      p_adjusted =
        p.adjust(
          p_raw,
          method = "holm"
        )
    )
  
  disp_tab <-
    bind_rows(
      disp_rows
    ) %>%
    mutate(
      p_adjusted =
        p.adjust(
          p_raw,
          method = "holm"
        )
    )
  
  list(
    permanova = perma_tab,
    permdisp = disp_tab
  )
}


# ============================================================
# 6. GLOBAL H1 / H2 PERMDISP
# ============================================================

run_H12_global_disp <- function(
    bray,
    meta,
    hypothesis,
    dataset,
    kingdom
) {
  
  meta <- standardise_H12_meta(
    bray,
    meta
  )
  
  perms <-
    make_site_permutations(
      meta
    )
  
  bd <- vegan::betadisper(
    bray,
    group =
      meta$Condition_master,
    type = "median",
    bias.adjust = TRUE
  )
  
  set.seed(123)
  
  pd <- vegan::permutest(
    bd,
    permutations = perms
  )
  
  tab <-
    as.data.frame(
      pd$tab
    )
  
  tibble(
    Hypothesis = hypothesis,
    Dataset = dataset,
    Kingdom = kingdom,
    Test_type = "PERMDISP",
    Model =
      "Global dispersion diagnostic",
    Term =
      "Ecological condition",
    Comparison =
      "Overall",
    n =
      nrow(meta),
    df =
      tab$Df[1],
    pseudo_F =
      tab$F[1],
    R2 = NA_real_,
    Adjusted_R2 = NA_real_,
    p_raw =
      tab$`Pr(>F)`[1],
    p_adjusted = NA_real_,
    Adjustment = "None",
    Inferential_status =
      "Dispersion diagnostic",
    Permutation_structure =
      "9999 permutations restricted within Site"
  )
}


# ============================================================
# 7. H3 PAIRED PERMANOVA
#
# TWO MODELS:
#
# A. Overall BS vs RZ
#    Site + Condition + Habitat
#
# B. Interaction
#    Site + Condition * Habitat
#
# Permutations restricted within PlotID.
# ============================================================

run_H3 <- function(
    bray,
    meta,
    kingdom
) {
  
  meta <-
    standardise_H3_meta(
      bray,
      meta
    )
  
  perms <-
    make_plot_permutations(
      meta
    )
  
  # ----------------------------------------------------------
  # H3-A: overall habitat model
  # ----------------------------------------------------------
  
  set.seed(123)
  
  fit_overall <- vegan::adonis2(
    bray ~
      Site_master +
      Condition_master +
      Habitat_master,
    data = meta,
    permutations = perms,
    by = "terms"
  )
  
  site_o <-
    extract_row(
      fit_overall,
      "Site_master"
    )
  
  cond_o <-
    extract_row(
      fit_overall,
      "Condition_master"
    )
  
  habitat_o <-
    extract_row(
      fit_overall,
      "Habitat_master"
    )
  
  R2_overall <-
    site_o$R2 +
    cond_o$R2 +
    habitat_o$R2
  
  db_overall <- vegan::dbrda(
    bray ~
      Site_master +
      Condition_master +
      Habitat_master,
    data = meta
  )
  
  adj_overall <-
    vegan::RsquareAdj(
      db_overall
    )$adj.r.squared
  
  # ----------------------------------------------------------
  # H3-B: interaction model
  # ----------------------------------------------------------
  
  set.seed(123)
  
  fit_interaction <- vegan::adonis2(
    bray ~
      Site_master +
      Condition_master *
      Habitat_master,
    data = meta,
    permutations = perms,
    by = "terms"
  )
  
  site_i <-
    extract_row(
      fit_interaction,
      "Site_master"
    )
  
  cond_i <-
    extract_row(
      fit_interaction,
      "Condition_master"
    )
  
  hab_i <-
    extract_row(
      fit_interaction,
      "Habitat_master"
    )
  
  int_term <-
    grep(
      ":",
      rownames(
        as.data.frame(
          fit_interaction
        )
      ),
      value = TRUE
    )[1]
  
  int_i <-
    extract_row(
      fit_interaction,
      int_term
    )
  
  R2_interaction <-
    site_i$R2 +
    cond_i$R2 +
    hab_i$R2 +
    int_i$R2
  
  db_interaction <- vegan::dbrda(
    bray ~
      Site_master +
      Condition_master *
      Habitat_master,
    data = meta
  )
  
  adj_interaction <-
    vegan::RsquareAdj(
      db_interaction
    )$adj.r.squared
  
  # ----------------------------------------------------------
  # Master rows
  # ----------------------------------------------------------
  
  omnibus_rows <- bind_rows(
    
    # Descriptive Site contribution from full model
    tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + condition × soil habitat",
      Term = "Site",
      Comparison = NA_character_,
      n = nrow(meta),
      df = site_i$df,
      pseudo_F = site_i$pseudo_F,
      R2 = site_i$R2,
      Adjusted_R2 = NA_real_,
      p_raw = NA_real_,
      p_adjusted = NA_real_,
      Adjustment = NA_character_,
      Inferential_status =
        "Descriptive between-plot contribution",
      Permutation_structure =
        "9999 permutations restricted within PlotID"
    ),
    
    # Descriptive condition contribution
    tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + condition × soil habitat",
      Term =
        "Ecological condition",
      Comparison = NA_character_,
      n = nrow(meta),
      df = cond_i$df,
      pseudo_F = cond_i$pseudo_F,
      R2 = cond_i$R2,
      Adjusted_R2 = NA_real_,
      p_raw = NA_real_,
      p_adjusted = NA_real_,
      Adjustment = NA_character_,
      Inferential_status =
        "Descriptive between-plot contribution",
      Permutation_structure =
        "9999 permutations restricted within PlotID"
    ),
    
    # Primary overall habitat test
    tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + condition + soil habitat",
      Term =
        "Soil habitat",
      Comparison =
        "Overall BS vs RZ",
      n = nrow(meta),
      df = habitat_o$df,
      pseudo_F =
        habitat_o$pseudo_F,
      R2 =
        habitat_o$R2,
      Adjusted_R2 = NA_real_,
      p_raw =
        habitat_o$p_raw,
      p_adjusted = NA_real_,
      Adjustment = "None",
      Inferential_status =
        "Primary inferential test",
      Permutation_structure =
        "9999 permutations restricted within PlotID"
    ),
    
    # Primary interaction
    tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + condition × soil habitat",
      Term =
        "Ecological condition × soil habitat",
      Comparison = NA_character_,
      n = nrow(meta),
      df =
        int_i$df,
      pseudo_F =
        int_i$pseudo_F,
      R2 =
        int_i$R2,
      Adjusted_R2 = NA_real_,
      p_raw =
        int_i$p_raw,
      p_adjusted = NA_real_,
      Adjustment = "None",
      Inferential_status =
        "Primary H3 interaction test",
      Permutation_structure =
        "9999 permutations restricted within PlotID"
    ),
    
    # Reduced whole model
    tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + condition + soil habitat",
      Term = "Whole model",
      Comparison =
        "Overall habitat model",
      n = nrow(meta),
      df =
        site_o$df +
        cond_o$df +
        habitat_o$df,
      pseudo_F = NA_real_,
      R2 =
        R2_overall,
      Adjusted_R2 =
        adj_overall,
      p_raw = NA_real_,
      p_adjusted = NA_real_,
      Adjustment = NA_character_,
      Inferential_status =
        "Whole-model variance explained",
      Permutation_structure =
        "Complete reduced model"
    ),
    
    # Full interaction whole model
    tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Site + condition × soil habitat",
      Term = "Whole model",
      Comparison =
        "Interaction model",
      n = nrow(meta),
      df =
        site_i$df +
        cond_i$df +
        hab_i$df +
        int_i$df,
      pseudo_F = NA_real_,
      R2 =
        R2_interaction,
      Adjusted_R2 =
        adj_interaction,
      p_raw = NA_real_,
      p_adjusted = NA_real_,
      Adjustment = NA_character_,
      Inferential_status =
        "Whole-model variance explained",
      Permutation_structure =
        "Complete interaction model"
    )
  )
  
  
  # ==========================================================
  # H3 WITHIN-CONDITION BS vs RZ
  # ==========================================================
  
  conditions <-
    c(
      "Natural",
      "BESE",
      "BARE"
    )
  
  bray_mat <-
    as.matrix(
      bray
    )
  
  within_rows <- list()
  
  for (
    i in seq_along(conditions)
  ) {
    
    cond <- conditions[i]
    
    meta_sub <-
      meta %>%
      filter(
        Condition_master == cond
      ) %>%
      droplevels()
    
    ids <-
      rownames(
        meta_sub
      )
    
    bray_sub <-
      as.dist(
        bray_mat[
          ids,
          ids
        ]
      )
    
    perms_sub <-
      make_plot_permutations(
        meta_sub
      )
    
    set.seed(123)
    
    fit_sub <- vegan::adonis2(
      bray_sub ~
        Site_master +
        Habitat_master,
      data = meta_sub,
      permutations = perms_sub,
      by = "terms"
    )
    
    hab_row <-
      extract_row(
        fit_sub,
        "Habitat_master"
      )
    
    within_rows[[i]] <- tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMANOVA",
      Model =
        "Within-condition paired BS vs RZ",
      Term =
        "Soil habitat",
      Comparison =
        paste0(
          cond,
          ": BS vs RZ"
        ),
      n =
        nrow(meta_sub),
      df =
        hab_row$df,
      pseudo_F =
        hab_row$pseudo_F,
      R2 =
        hab_row$R2,
      Adjusted_R2 =
        NA_real_,
      p_raw =
        hab_row$p_raw,
      p_adjusted =
        NA_real_,
      Adjustment =
        "Holm",
      Inferential_status =
        "Planned within-condition test",
      Permutation_structure =
        "9999 permutations restricted within PlotID"
    )
  }
  
  within_table <-
    bind_rows(
      within_rows
    ) %>%
    mutate(
      p_adjusted =
        p.adjust(
          p_raw,
          method = "holm"
        )
    )
  
  list(
    omnibus =
      omnibus_rows,
    within =
      within_table,
    meta =
      meta,
    fit_overall =
      fit_overall,
    fit_interaction =
      fit_interaction
  )
}


# ============================================================
# 8. H3 PERMDISP
#
# Matches previous diagnostic strategy:
# unrestricted PERMDISP because PERMDISP does not reproduce
# the paired PlotID PERMANOVA design.
# ============================================================

run_H3_permdisp <- function(
    bray,
    meta,
    kingdom
) {
  
  meta <-
    standardise_H3_meta(
      bray,
      meta
    )
  
  meta$Condition_habitat <-
    interaction(
      meta$Condition_master,
      meta$Habitat_master,
      sep = " | ",
      drop = TRUE
    )
  
  # ------------------------------
  # Global six-group PERMDISP
  # ------------------------------
  
  bd_global <-
    vegan::betadisper(
      bray,
      group =
        meta$Condition_habitat,
      type = "median",
      bias.adjust = TRUE
    )
  
  set.seed(123)
  
  pd_global <-
    vegan::permutest(
      bd_global,
      permutations = 9999
    )
  
  global_tab <-
    as.data.frame(
      pd_global$tab
    )
  
  global_row <- tibble(
    Hypothesis = "H3",
    Dataset =
      "Paired bulk soil–rhizosphere",
    Kingdom = kingdom,
    Test_type = "PERMDISP",
    Model =
      "Six condition × habitat groups",
    Term =
      "Multivariate dispersion",
    Comparison =
      "Overall",
    n = nrow(meta),
    df =
      global_tab$Df[1],
    pseudo_F =
      global_tab$F[1],
    R2 = NA_real_,
    Adjusted_R2 = NA_real_,
    p_raw =
      global_tab$`Pr(>F)`[1],
    p_adjusted = NA_real_,
    Adjustment = "None",
    Inferential_status =
      "Dispersion diagnostic",
    Permutation_structure =
      "9999 unrestricted permutations"
  )
  
  # ------------------------------
  # Within condition
  # ------------------------------
  
  conditions <-
    c(
      "Natural",
      "BESE",
      "BARE"
    )
  
  bray_mat <-
    as.matrix(
      bray
    )
  
  rows <- list()
  
  for (
    i in seq_along(conditions)
  ) {
    
    cond <- conditions[i]
    
    meta_sub <-
      meta %>%
      filter(
        Condition_master == cond
      ) %>%
      droplevels()
    
    ids <-
      rownames(
        meta_sub
      )
    
    bray_sub <-
      as.dist(
        bray_mat[
          ids,
          ids
        ]
      )
    
    bd <- vegan::betadisper(
      bray_sub,
      group =
        meta_sub$Habitat_master,
      type = "median",
      bias.adjust = TRUE
    )
    
    set.seed(123)
    
    pd <- vegan::permutest(
      bd,
      permutations = 9999
    )
    
    tab <-
      as.data.frame(
        pd$tab
      )
    
    rows[[i]] <- tibble(
      Hypothesis = "H3",
      Dataset =
        "Paired bulk soil–rhizosphere",
      Kingdom = kingdom,
      Test_type = "PERMDISP",
      Model =
        "Within-condition BS vs RZ dispersion",
      Term =
        "Soil habitat",
      Comparison =
        paste0(
          cond,
          ": BS vs RZ"
        ),
      n =
        nrow(meta_sub),
      df =
        tab$Df[1],
      pseudo_F =
        tab$F[1],
      R2 = NA_real_,
      Adjusted_R2 = NA_real_,
      p_raw =
        tab$`Pr(>F)`[1],
      p_adjusted =
        NA_real_,
      Adjustment =
        "Holm",
      Inferential_status =
        "Dispersion diagnostic",
      Permutation_structure =
        "9999 unrestricted permutations"
    )
  }
  
  within <-
    bind_rows(
      rows
    ) %>%
    mutate(
      p_adjusted =
        p.adjust(
          p_raw,
          method = "holm"
        )
    )
  
  bind_rows(
    global_row,
    within
  )
}


# ============================================================
# 9. RUN H1
# ============================================================

H1_prok <-
  run_H12_omnibus(
    bray_H1_prok,
    meta_H1_prok,
    hypothesis = "H1",
    dataset =
      "Bulk soil: Natural, BESE, BARE, Degraded",
    kingdom =
      "Prokaryotic"
  )

H1_fun <-
  run_H12_omnibus(
    bray_H1_fun,
    meta_H1_fun,
    hypothesis = "H1",
    dataset =
      "Bulk soil: Natural, BESE, BARE, Degraded",
    kingdom =
      "Fungal"
  )


# Planned H1 pairwise family
H1_pairs <- list(
  c(
    "Natural",
    "Degraded"
  ),
  c(
    "Natural",
    "BESE"
  ),
  c(
    "Natural",
    "BARE"
  ),
  c(
    "BESE",
    "BARE"
  )
)

H1_pair_prok <-
  run_H12_pairwise(
    bray_H1_prok,
    meta_H1_prok,
    "H1",
    "Bulk soil: Natural, BESE, BARE, Degraded",
    "Prokaryotic",
    H1_pairs
  )

H1_pair_fun <-
  run_H12_pairwise(
    bray_H1_fun,
    meta_H1_fun,
    "H1",
    "Bulk soil: Natural, BESE, BARE, Degraded",
    "Fungal",
    H1_pairs
  )


H1_disp_prok <-
  run_H12_global_disp(
    bray_H1_prok,
    meta_H1_prok,
    "H1",
    "Bulk soil: Natural, BESE, BARE, Degraded",
    "Prokaryotic"
  )

H1_disp_fun <-
  run_H12_global_disp(
    bray_H1_fun,
    meta_H1_fun,
    "H1",
    "Bulk soil: Natural, BESE, BARE, Degraded",
    "Fungal"
  )


# ============================================================
# 10. RUN H2
# ============================================================

H2_prok <-
  run_H12_omnibus(
    bray_H2_prok,
    meta_H2_prok,
    hypothesis = "H2",
    dataset =
      "Rhizosphere: Natural, BESE, BARE",
    kingdom =
      "Prokaryotic"
  )

H2_fun <-
  run_H12_omnibus(
    bray_H2_fun,
    meta_H2_fun,
    hypothesis = "H2",
    dataset =
      "Rhizosphere: Natural, BESE, BARE",
    kingdom =
      "Fungal"
  )


H2_pairs <- list(
  c(
    "Natural",
    "BESE"
  ),
  c(
    "Natural",
    "BARE"
  ),
  c(
    "BESE",
    "BARE"
  )
)

H2_pair_prok <-
  run_H12_pairwise(
    bray_H2_prok,
    meta_H2_prok,
    "H2",
    "Rhizosphere: Natural, BESE, BARE",
    "Prokaryotic",
    H2_pairs
  )

H2_pair_fun <-
  run_H12_pairwise(
    bray_H2_fun,
    meta_H2_fun,
    "H2",
    "Rhizosphere: Natural, BESE, BARE",
    "Fungal",
    H2_pairs
  )


H2_disp_prok <-
  run_H12_global_disp(
    bray_H2_prok,
    meta_H2_prok,
    "H2",
    "Rhizosphere: Natural, BESE, BARE",
    "Prokaryotic"
  )

H2_disp_fun <-
  run_H12_global_disp(
    bray_H2_fun,
    meta_H2_fun,
    "H2",
    "Rhizosphere: Natural, BESE, BARE",
    "Fungal"
  )


# ============================================================
# 11. RUN H3
# ============================================================

H3_prok <-
  run_H3(
    bray_H3_prok,
    meta_H3_prok,
    "Prokaryotic"
  )

H3_fun <-
  run_H3(
    bray_H3_fun,
    meta_H3_fun,
    "Fungal"
  )


H3_disp_prok <-
  run_H3_permdisp(
    bray_H3_prok,
    meta_H3_prok,
    "Prokaryotic"
  )

H3_disp_fun <-
  run_H3_permdisp(
    bray_H3_fun,
    meta_H3_fun,
    "Fungal"
  )


# ============================================================
# 12. BUILD MASTER PERMANOVA TABLE
# ============================================================

beta_PERMANOVA_master <-
  bind_rows(
    
    # H1 omnibus
    H1_prok$table,
    H1_fun$table,
    
    # H1 pairwise
    H1_pair_prok$permanova,
    H1_pair_fun$permanova,
    
    # H2 omnibus
    H2_prok$table,
    H2_fun$table,
    
    # H2 pairwise
    H2_pair_prok$permanova,
    H2_pair_fun$permanova,
    
    # H3 omnibus / interaction
    H3_prok$omnibus,
    H3_fun$omnibus,
    
    # H3 within condition
    H3_prok$within,
    H3_fun$within
    
  ) %>%
  mutate(
    R2_percent =
      100 * R2,
    
    Adjusted_R2_percent =
      100 * Adjusted_R2
  )


# ============================================================
# 13. BUILD MASTER PERMDISP TABLE
# ============================================================

beta_PERMDISP_master <-
  bind_rows(
    
    # H1 global
    H1_disp_prok,
    H1_disp_fun,
    
    # H1 pairwise
    H1_pair_prok$permdisp,
    H1_pair_fun$permdisp,
    
    # H2 global
    H2_disp_prok,
    H2_disp_fun,
    
    # H2 pairwise
    H2_pair_prok$permdisp,
    H2_pair_fun$permdisp,
    
    # H3
    H3_disp_prok,
    H3_disp_fun
    
  ) %>%
  mutate(
    R2_percent =
      NA_real_,
    
    Adjusted_R2_percent =
      NA_real_
  )


# ============================================================
# 14. ONE COMPLETE BETA MASTER TABLE
# ============================================================

beta_H1_H2_H3_MASTER <-
  bind_rows(
    beta_PERMANOVA_master,
    beta_PERMDISP_master
  ) %>%
  select(
    Hypothesis,
    Dataset,
    Kingdom,
    Test_type,
    Model,
    Term,
    Comparison,
    n,
    df,
    pseudo_F,
    R2,
    R2_percent,
    Adjusted_R2,
    Adjusted_R2_percent,
    p_raw,
    p_adjusted,
    Adjustment,
    Inferential_status,
    Permutation_structure
  )


# ============================================================
# 15. MAIN-TEXT STATISTICS ONLY
#
# Includes:
# H1/H2:
#   Site contribution
#   Ecological condition
#   Whole model
#
# H3:
#   Site contribution
#   Ecological condition contribution
#   Overall soil habitat
#   Condition × habitat interaction
#   Full whole model
# ============================================================

beta_MAIN_TEXT_stats <-
  beta_PERMANOVA_master %>%
  filter(
    
    (
      Hypothesis %in%
        c(
          "H1",
          "H2"
        ) &
        is.na(Comparison) &
        Term %in%
        c(
          "Site",
          "Ecological condition",
          "Whole model"
        )
    ) |
      
      (
        Hypothesis == "H3" &
          (
            Term %in%
              c(
                "Site",
                "Ecological condition",
                "Ecological condition × soil habitat"
              ) |
              
              (
                Term == "Soil habitat" &
                  Comparison ==
                  "Overall BS vs RZ"
              ) |
              
              (
                Term == "Whole model" &
                  Comparison ==
                  "Interaction model"
              )
          )
      )
  ) %>%
  arrange(
    Hypothesis,
    Kingdom,
    factor(
      Term,
      levels = c(
        "Site",
        "Ecological condition",
        "Soil habitat",
        "Ecological condition × soil habitat",
        "Whole model"
      )
    )
  )


# ============================================================
# 16. SAVE
# ============================================================

dir.create(
  "outputs/beta_master",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  beta_PERMANOVA_master,
  "outputs/beta_master/Beta_H1_H2_H3_MASTER_PERMANOVA.csv",
  row.names = FALSE
)

write.csv(
  beta_PERMDISP_master,
  "outputs/beta_master/Beta_H1_H2_H3_MASTER_PERMDISP.csv",
  row.names = FALSE
)

write.csv(
  beta_H1_H2_H3_MASTER,
  "outputs/beta_master/Beta_H1_H2_H3_MASTER_ALL.csv",
  row.names = FALSE
)

write.csv(
  beta_MAIN_TEXT_stats,
  "outputs/beta_master/Beta_H1_H2_H3_MAIN_TEXT_STATS.csv",
  row.names = FALSE
)

saveRDS(
  list(
    PERMANOVA =
      beta_PERMANOVA_master,
    PERMDISP =
      beta_PERMDISP_master,
    MASTER =
      beta_H1_H2_H3_MASTER,
    MAIN_TEXT =
      beta_MAIN_TEXT_stats
  ),
  "outputs/beta_master/Beta_H1_H2_H3_MASTER_RESULTS.rds"
)


# ============================================================
# 17. PRINT ONLY WHAT WE NEED TO CHECK
# ============================================================

cat("\n========================================\n")
cat("MAIN-TEXT BETA STATISTICS\n")
cat("========================================\n\n")

print(
  beta_MAIN_TEXT_stats %>%
    select(
      Hypothesis,
      Kingdom,
      Model,
      Term,
      Comparison,
      n,
      df,
      pseudo_F,
      R2_percent,
      Adjusted_R2_percent,
      p_raw,
      Inferential_status
    ),
  n = Inf,
  width = Inf
)


cat("\n========================================\n")
cat("H1 PAIRWISE PERMANOVA\n")
cat("========================================\n\n")

print(
  beta_PERMANOVA_master %>%
    filter(
      Hypothesis == "H1",
      !is.na(Comparison),
      grepl(
        "Pairwise",
        Model
      )
    ) %>%
    select(
      Kingdom,
      Comparison,
      pseudo_F,
      R2_percent,
      p_raw,
      p_adjusted
    ),
  n = Inf,
  width = Inf
)


cat("\n========================================\n")
cat("H2 PAIRWISE PERMANOVA\n")
cat("========================================\n\n")

print(
  beta_PERMANOVA_master %>%
    filter(
      Hypothesis == "H2",
      !is.na(Comparison)
    ) %>%
    select(
      Kingdom,
      Comparison,
      pseudo_F,
      R2_percent,
      p_raw,
      p_adjusted
    ),
  n = Inf,
  width = Inf
)


cat("\n========================================\n")
cat("H3 WITHIN-CONDITION PERMANOVA\n")
cat("========================================\n\n")

print(
  beta_PERMANOVA_master %>%
    filter(
      Hypothesis == "H3",
      Model ==
        "Within-condition paired BS vs RZ"
    ) %>%
    select(
      Kingdom,
      Comparison,
      pseudo_F,
      R2_percent,
      p_raw,
      p_adjusted
    ),
  n = Inf,
  width = Inf
)


cat("\n========================================\n")
cat("FILES SAVED\n")
cat("========================================\n")

print(
  list.files(
    "outputs/beta_master",
    full.names = FALSE
  )
)


# ============================================================
# SUPPLEMENTARY TABLE S7
# ONE EXCEL WORKBOOK WITH THREE TABS
#
# S7a = Omnibus PERMANOVA
# S7b = Pairwise / within-condition PERMANOVA
# S7c = PERMDISP diagnostics
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(openxlsx)
})

# ============================================================
# 0. CHECK REQUIRED S7a OBJECT
# ============================================================

stopifnot(
  exists("S7a_beta_omnibus")
)


# ============================================================
# 1. S7a — OMNIBUS PERMANOVA
# ============================================================

S7a <- S7a_beta_omnibus %>%
  mutate(
    pseudo_F = round(pseudo_F, 4),
    R2 = round(R2, 5),
    R2_percent = round(R2_percent, 2)
  ) %>%
  select(
    Hypothesis,
    Kingdom,
    Dataset,
    Model,
    Term,
    n,
    df,
    pseudo_F,
    R2,
    R2_percent,
    p_display,
    Inference
  )


# ============================================================
# 2. S7b — PAIRWISE / WITHIN-CONDITION PERMANOVA
# ============================================================

S7b_parts <- list()


# ------------------------------------------------------------
# H1 — BULK SOIL PAIRWISE
# Existing final object: pairwise_results
# ------------------------------------------------------------

if (exists("pairwise_results")) {
  
  H1_pair <- as.data.frame(
    pairwise_results
  )
  
  # Flexible column finder
  get_col <- function(dat, candidates) {
    x <- candidates[candidates %in% names(dat)]
    if (length(x) == 0) return(rep(NA, nrow(dat)))
    dat[[x[1]]]
  }
  
  S7b_parts[["H1"]] <- tibble(
    Hypothesis = "H1",
    Kingdom = as.character(
      get_col(
        H1_pair,
        c("Kingdom")
      )
    ),
    Dataset =
      "Bulk soil",
    Comparison = as.character(
      get_col(
        H1_pair,
        c("Comparison", "Contrast")
      )
    ),
    n = as.numeric(
      get_col(
        H1_pair,
        c("n", "N")
      )
    ),
    df = 1,
    pseudo_F = as.numeric(
      get_col(
        H1_pair,
        c(
          "PERMANOVA_F",
          "Pseudo_F",
          "F"
        )
      )
    ),
    R2 = as.numeric(
      get_col(
        H1_pair,
        c(
          "PERMANOVA_R2",
          "R2"
        )
      )
    ),
    p_raw = as.numeric(
      get_col(
        H1_pair,
        c(
          "PERMANOVA_p",
          "p_value",
          "p"
        )
      )
    ),
    p_adjusted = as.numeric(
      get_col(
        H1_pair,
        c(
          "PERMANOVA_p_Holm",
          "p_Holm"
        )
      )
    ),
    Adjustment = "Holm",
    Permutation_structure =
      "9999 permutations restricted within Site"
  )
}


# ------------------------------------------------------------
# H2 — RHIZOSPHERE PAIRWISE
# Existing final object: pairwise_beta
# ------------------------------------------------------------

if (exists("pairwise_beta")) {
  
  H2_pair <- as.data.frame(
    pairwise_beta
  )
  
  S7b_parts[["H2"]] <- tibble(
    Hypothesis = "H2",
    Kingdom = case_when(
      H2_pair$Kingdom %in%
        c("Bacteria", "Prokaryotic") ~
        "Prokaryotic",
      
      H2_pair$Kingdom == "Fungi" ~
        "Fungal",
      
      TRUE ~
        as.character(
          H2_pair$Kingdom
        )
    ),
    Dataset =
      "Rhizosphere",
    Comparison =
      H2_pair$Contrast,
    n =
      H2_pair$n,
    df = 1,
    pseudo_F =
      H2_pair$PERMANOVA_F,
    R2 =
      H2_pair$PERMANOVA_R2,
    p_raw =
      H2_pair$PERMANOVA_p,
    p_adjusted =
      H2_pair$PERMANOVA_p_Holm,
    Adjustment =
      "Holm",
    Permutation_structure =
      "9999 permutations restricted within Site"
  )
}


# ------------------------------------------------------------
# H3 — PROKARYOTIC BS vs RZ WITHIN CONDITION
# Existing final object: within_table
# ------------------------------------------------------------

if (exists("within_table")) {
  
  H3_prok_within <-
    as.data.frame(
      within_table
    )
  
  S7b_parts[["H3_prok"]] <- tibble(
    Hypothesis = "H3",
    Kingdom =
      "Prokaryotic",
    Dataset =
      "Paired bulk soil–rhizosphere",
    Comparison =
      paste0(
        H3_prok_within$Ecological_condition,
        ": BS vs RZ"
      ),
    n = NA_real_,
    df = 1,
    pseudo_F =
      H3_prok_within$F,
    R2 =
      H3_prok_within$R2,
    p_raw =
      H3_prok_within$p_raw,
    p_adjusted =
      H3_prok_within$p_Holm,
    Adjustment =
      "Holm",
    Permutation_structure =
      "9999 permutations restricted within PlotID"
  )
}


# ------------------------------------------------------------
# H3 — FUNGAL BS vs RZ WITHIN CONDITION
# Existing final object: within_table_fun
# ------------------------------------------------------------

if (exists("within_table_fun")) {
  
  H3_fun_within <-
    as.data.frame(
      within_table_fun
    )
  
  S7b_parts[["H3_fun"]] <- tibble(
    Hypothesis = "H3",
    Kingdom =
      "Fungal",
    Dataset =
      "Paired bulk soil–rhizosphere",
    Comparison =
      paste0(
        H3_fun_within$Ecological_condition,
        ": BS vs RZ"
      ),
    n = NA_real_,
    df = 1,
    pseudo_F =
      H3_fun_within$F,
    R2 =
      H3_fun_within$R2,
    p_raw =
      H3_fun_within$p_raw,
    p_adjusted =
      H3_fun_within$p_Holm,
    Adjustment =
      "Holm",
    Permutation_structure =
      "9999 permutations restricted within PlotID"
  )
}


S7b <- bind_rows(
  S7b_parts
) %>%
  mutate(
    R2_percent =
      100 * R2,
    
    pseudo_F =
      round(
        pseudo_F,
        4
      ),
    
    R2 =
      round(
        R2,
        5
      ),
    
    R2_percent =
      round(
        R2_percent,
        2
      ),
    
    p_display = case_when(
      is.na(p_raw) ~ "—",
      p_raw < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.4f",
        p_raw
      )
    ),
    
    p_adjusted_display = case_when(
      is.na(p_adjusted) ~ "—",
      p_adjusted < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.4f",
        p_adjusted
      )
    )
  ) %>%
  select(
    Hypothesis,
    Kingdom,
    Dataset,
    Comparison,
    n,
    df,
    pseudo_F,
    R2,
    R2_percent,
    p_display,
    p_adjusted_display,
    Adjustment,
    Permutation_structure
  )


# ============================================================
# 3. S7c — PERMDISP
# ============================================================

S7c_parts <- list()


# ------------------------------------------------------------
# H1 PERMDISP
# pairwise_results already contains pairwise PERMDISP
# ------------------------------------------------------------

if (exists("pairwise_results")) {
  
  H1_disp <- as.data.frame(
    pairwise_results
  )
  
  S7c_parts[["H1_pairwise"]] <- tibble(
    Hypothesis = "H1",
    Kingdom = as.character(
      get_col(
        H1_disp,
        c("Kingdom")
      )
    ),
    Dataset =
      "Bulk soil",
    Analysis =
      "Pairwise PERMDISP",
    Comparison = as.character(
      get_col(
        H1_disp,
        c("Comparison", "Contrast")
      )
    ),
    pseudo_F = as.numeric(
      get_col(
        H1_disp,
        c("PERMDISP_F")
      )
    ),
    p_raw = as.numeric(
      get_col(
        H1_disp,
        c("PERMDISP_p")
      )
    ),
    p_adjusted = as.numeric(
      get_col(
        H1_disp,
        c("PERMDISP_p_Holm")
      )
    ),
    Adjustment =
      "Holm",
    Permutation_structure =
      "9999 permutations restricted within Site"
  )
}


# ------------------------------------------------------------
# H2 PERMDISP
# pairwise_beta contains pairwise PERMDISP
# ------------------------------------------------------------

if (exists("pairwise_beta")) {
  
  H2_disp <- as.data.frame(
    pairwise_beta
  )
  
  S7c_parts[["H2_pairwise"]] <- tibble(
    Hypothesis = "H2",
    Kingdom = case_when(
      H2_disp$Kingdom %in%
        c("Bacteria", "Prokaryotic") ~
        "Prokaryotic",
      
      H2_disp$Kingdom == "Fungi" ~
        "Fungal",
      
      TRUE ~
        as.character(
          H2_disp$Kingdom
        )
    ),
    Dataset =
      "Rhizosphere",
    Analysis =
      "Pairwise PERMDISP",
    Comparison =
      H2_disp$Contrast,
    pseudo_F =
      H2_disp$PERMDISP_F,
    p_raw =
      H2_disp$PERMDISP_p,
    p_adjusted =
      H2_disp$PERMDISP_p_Holm,
    Adjustment =
      "Holm",
    Permutation_structure =
      "9999 permutations restricted within Site"
  )
}


# ------------------------------------------------------------
# H3 PROKARYOTIC WITHIN-CONDITION PERMDISP
# ------------------------------------------------------------

if (exists("permdisp_summary")) {
  
  H3_prok_disp <- as.data.frame(
    permdisp_summary
  )
  
  S7c_parts[["H3_prok"]] <- tibble(
    Hypothesis = "H3",
    Kingdom =
      "Prokaryotic",
    Dataset =
      "Paired bulk soil–rhizosphere",
    Analysis =
      "Within-condition PERMDISP",
    Comparison =
      paste0(
        H3_prok_disp$Ecological_condition,
        ": BS vs RZ"
      ),
    pseudo_F =
      H3_prok_disp$F,
    p_raw =
      H3_prok_disp$p_raw,
    p_adjusted =
      if (
        "p_Holm" %in%
        names(H3_prok_disp)
      ) {
        H3_prok_disp$p_Holm
      } else {
        p.adjust(
          H3_prok_disp$p_raw,
          method = "holm"
        )
      },
    Adjustment =
      "Holm",
    Permutation_structure =
      "PERMDISP diagnostic; unrestricted permutations"
  )
}


# ------------------------------------------------------------
# H3 FUNGAL WITHIN-CONDITION PERMDISP
# ------------------------------------------------------------

if (exists("permdisp_table_fun")) {
  
  H3_fun_disp <- as.data.frame(
    permdisp_table_fun
  )
  
  S7c_parts[["H3_fun"]] <- tibble(
    Hypothesis = "H3",
    Kingdom =
      "Fungal",
    Dataset =
      "Paired bulk soil–rhizosphere",
    Analysis =
      "Within-condition PERMDISP",
    Comparison =
      paste0(
        H3_fun_disp$Ecological_condition,
        ": BS vs RZ"
      ),
    pseudo_F =
      H3_fun_disp$F,
    p_raw =
      H3_fun_disp$p_raw,
    p_adjusted =
      if (
        "p_Holm" %in%
        names(H3_fun_disp)
      ) {
        H3_fun_disp$p_Holm
      } else {
        p.adjust(
          H3_fun_disp$p_raw,
          method = "holm"
        )
      },
    Adjustment =
      "Holm",
    Permutation_structure =
      "PERMDISP diagnostic; unrestricted permutations"
  )
}


S7c <- bind_rows(
  S7c_parts
) %>%
  mutate(
    pseudo_F =
      round(
        pseudo_F,
        4
      ),
    
    p_display = case_when(
      is.na(p_raw) ~ "—",
      p_raw < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.4f",
        p_raw
      )
    ),
    
    p_adjusted_display = case_when(
      is.na(p_adjusted) ~ "—",
      p_adjusted < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.4f",
        p_adjusted
      )
    )
  ) %>%
  select(
    Hypothesis,
    Kingdom,
    Dataset,
    Analysis,
    Comparison,
    pseudo_F,
    p_display,
    p_adjusted_display,
    Adjustment,
    Permutation_structure
  )


# ============================================================
# 4. CREATE ONE EXCEL WORKBOOK
# ============================================================

dir.create(
  "outputs/supplementary_tables",
  recursive = TRUE,
  showWarnings = FALSE
)

wb <- createWorkbook()


# ------------------------------------------------------------
# SHEET S7a
# ------------------------------------------------------------

addWorksheet(
  wb,
  "S7a_Omnibus_PERMANOVA"
)

writeData(
  wb,
  "S7a_Omnibus_PERMANOVA",
  S7a,
  withFilter = TRUE
)


# ------------------------------------------------------------
# SHEET S7b
# ------------------------------------------------------------

addWorksheet(
  wb,
  "S7b_Pairwise_PERMANOVA"
)

writeData(
  wb,
  "S7b_Pairwise_PERMANOVA",
  S7b,
  withFilter = TRUE
)


# ------------------------------------------------------------
# SHEET S7c
# ------------------------------------------------------------

addWorksheet(
  wb,
  "S7c_PERMDISP"
)

writeData(
  wb,
  "S7c_PERMDISP",
  S7c,
  withFilter = TRUE
)


# ============================================================
# 5. BASIC FORMATTING
# ============================================================

header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  valign = "center",
  wrapText = TRUE,
  border = "Bottom"
)

for (
  sheet in c(
    "S7a_Omnibus_PERMANOVA",
    "S7b_Pairwise_PERMANOVA",
    "S7c_PERMDISP"
  )
) {
  
  addStyle(
    wb,
    sheet,
    style = header_style,
    rows = 1,
    cols = 1:ncol(
      if (
        sheet == "S7a_Omnibus_PERMANOVA"
      ) {
        S7a
      } else if (
        sheet == "S7b_Pairwise_PERMANOVA"
      ) {
        S7b
      } else {
        S7c
      }
    ),
    gridExpand = TRUE
  )
  
  freezePane(
    wb,
    sheet,
    firstRow = TRUE
  )
  
  setColWidths(
    wb,
    sheet,
    cols = 1:20,
    widths = "auto"
  )
}


# ============================================================
# 6. SAVE WORKBOOK
# ============================================================

S7_file <-
  "outputs/supplementary_tables/Supplementary_Table_S7_Beta_diversity.xlsx"

saveWorkbook(
  wb,
  S7_file,
  overwrite = TRUE
)


# ============================================================
# 7. FINAL CHECK
# ============================================================

cat("\n========================================\n")
cat("SUPPLEMENTARY TABLE S7 CREATED\n")
cat("========================================\n\n")

cat(
  "Workbook:",
  S7_file,
  "\n\n"
)

cat(
  "S7a rows:",
  nrow(S7a),
  "\n"
)

cat(
  "S7b rows:",
  nrow(S7b),
  "\n"
)

cat(
  "S7c rows:",
  nrow(S7c),
  "\n"
)

cat("\nSheets:\n")
print(
  names(wb)
)