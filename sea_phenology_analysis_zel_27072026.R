#### Southeast Asian masting dataset - Analysis Script ####
#### Zoë Lieb ####
#### Last updated: July 2026 ####
#### re-cleaned version (27-07-2026)####

#### ----------------------------------------------------------------- ####
####                        SETUP AND LIBRARIES                         ####
#### ----------------------------------------------------------------- ####

library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggrepel)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(circular)
library(patchwork)
library(scales)
library(colorspace)
library(lme4)
library(broom.mixed)


#### ----------------------------------------------------------------- ####
####                         READ DATA                                  ####
#### ----------------------------------------------------------------- ####

#### Step 1 - Read harmonised dataset and literature review events ####


all_harmonised <- read_csv(
  "exported data summaries/sea_phenology_combined_data.csv",
  show_col_types = FALSE
)

# pheno: one row per discrete literature review masting event
# with explicit onset_date and end_date — used for onset analysis
pheno <- read_csv(
  "exported data summaries/pheno_lit_review_events.csv",
  show_col_types = FALSE
)



#### ----------------------------------------------------------------- ####
####                    LOOKUP TABLES AND CONSTANTS                     ####
#### ----------------------------------------------------------------- ####

#### Step 2 - Site name recode lookup ####
# Applied to pheno to match all_harmonised site names

pheno_site_recode <- c(
  "Bukit Timah Nature Reserve"                    = "Singapore",
  "Central Catchment Nature Reserve (MacRitchie)" = "Singapore",
  "Bukit Timah"                                   = "Singapore",
  "Lambir Hills National Park"                    = "Lambir Hills",
  "Kepong"                                        = "FRIM Malaysia - Kepong - Ampang",
  "FRIM Malaysia"                                 = "FRIM Malaysia - Kepong - Ampang",
  "Ampang"                                        = "FRIM Malaysia - Kepong - Ampang",
  "Ulu Segama Forest Reserve"                     = "Ulu Segama Malua"
)

# Sites to exclude from onset analysis (no clear onset date available)
excluded_sites <- c(
  "Doi Suthep-Pui National Park",
  "Khao Ang Runai Wildlife Sanctuary",
  "Khun Wang Royal Agriculture Research Center",
  "Lobo"
)

# Covidence IDs to exclude from onset analysis
# (only peak month recorded, not true onset date)
exclude_cov_ids <- c(85, 359, 405, 108, 474, 74, 17)

#### Step 3 - Subregion grouping lookup ####
# Used for circular plots and bubble chart

subregion_recode <- c(
  "Singapore"           = "Peninsular Malaysia & Singapore",
  "Peninsular Malaysia" = "Peninsular Malaysia & Singapore",
  "Sabah"               = "Borneo",
  "Sarawak"             = "Borneo",
  "West Kalimantan"     = "Borneo",
  "Central Kalimantan"  = "Borneo",
  "East Kalimantan"     = "Borneo",
  "Brunei"              = "Borneo"
)

#### Step 4 - Source type lookup ####
# Maps dataset names to human-readable source types

source_type_lookup <- c(
  "FRIM"                       = "Long-term monitoring",
  "Lambir"                     = "Long-term monitoring",
  "Gunung_Palung"              = "Long-term monitoring",
  "Pasoh"                      = "Long-term monitoring",
  "Bukit_Timah"                = "Long-term monitoring",
  "Gunung_Tarak"               = "Unpublished long-term",
  "BBBR"                       = "Unpublished long-term",
  "additional_contributed_obs" = "Expert observations",
  "Literature_review"          = "Literature review"
)

#### Step 5 - Color constants ####
# All figure colors defined here for consistency across the paper

# Phenophase colors
col_FL      <- "#6C739D"   # Glaucous — Flowering
col_FR      <- "#7FBEAB"   # Muted Teal — Fruiting
col_FL_sig  <- "#6C739D"   # Significant flowering (full color)
col_FL_nonsig <- "#B8BCE0" # Non-significant flowering (lighter)
col_FR_sig  <- "#7FBEAB"   # Significant fruiting (full color)
col_FR_nonsig <- "#B8DDD4" # Non-significant fruiting (lighter)
col_FL_border <- "#3D4270" # Dark purple border for significant FL bars
col_FR_border <- "#4A8C7A" # Dark teal border for significant FR bars

# Data source colors
col_ltm        <- "#8ECCA5"  # Celadon — Long-term monitoring
col_unpub_ltm  <- "#F0BE19"  # Saffron — Unpublished long-term
col_expert     <- "#D37D56"  # Burnt Peach — Expert observations
col_lit        <- "#708AFF"  # Cornflower Blue — Literature review

# Timeline phenophase colors (different from analysis colors for contrast)
col_timeline_fl <- "#916392"  # Medium purple — Flowering in timeline
col_timeline_fr <- "#371B38"  # Dark purple — Fruiting in timeline

# Convenience named vectors
phase_cols <- c("FL" = col_FL, "FR" = col_FR)
phase_labels <- c("FL" = "Flowering", "FR" = "Fruiting")

source_cols <- c(
  "Long-term monitoring"  = col_ltm,
  "Unpublished long-term" = col_unpub_ltm,
  "Expert observations"   = col_expert,
  "Literature review"     = col_lit
)


#### Theme constants ####

#### Custom theme — consistent across all manuscript figures ####
# Target: Ecology Letters full page width (173mm = 6.8 inches)
# base_size = 11 calibrated for this width

theme_masting <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = "grey92", linewidth = 0.3),
      axis.text         = element_text(size = base_size * 0.9),
      axis.title        = element_text(size = base_size),
      legend.text       = element_text(size = base_size * 0.9),
      legend.title      = element_text(size = base_size, face = "bold"),
      strip.text        = element_text(size = base_size * 0.9, face = "bold"),
      plot.title        = element_text(size = base_size * 1.1, face = "bold"),
      plot.subtitle     = element_text(size = base_size * 0.85,
                                       colour = "grey40"),
      plot.margin       = margin(6, 6, 6, 6)
    )
}

# Standard export dimensions for Ecology Letters full page width
fig_width  <- 6.8   # inches (173mm)
fig_height_single <- 4    # for single panel figures
fig_height_double <- 6    # for two-panel figures  
fig_height_triple <- 8    # for three-panel or tall figures
fig_dpi    <- 300





########################################################################
###                   Descriptive Data Visualisation                 ###
########################################################################

#### Timeline figure - all sites, colored by source type ####

# Source type lookup — matches data integration diagram categories
source_type_lookup <- c(
  "FRIM"                       = "Long-term monitoring",
  "Lambir"                     = "Long-term monitoring",
  "Gunung_Palung"              = "Long-term monitoring",
  "Pasoh"                      = "Long-term monitoring",
  "Bukit_Timah"                = "Long-term monitoring",
  "Gunung_Tarak"               = "Unpublished long-term",
  "BBBR"                       = "Unpublished long-term",
  "additional_contributed_obs" = "Expert observations",
  "Literature_review"          = "Literature review"
)

source_cols <- c(
  "Long-term monitoring"  = "#8ECCA5",   # Celadon
  "Unpublished long-term" = "#F0BE19",   # Saffron
  "Expert observations"   = "#D37D56",   # Burnt Peach
  "Literature review"     = "#708AFF"    # Cornflower Blue
)

# FL/FR phase colors — consistent with rest of paper
phase_cols_timeline <- c(
  "FL" = "#916392",
  "FR" = "#371B38"
)

# Build timeline dataframe
timeline_df <- all_harmonised %>%
  filter(!is.na(month_date)) %>%
  mutate(
    source_type = dplyr::recode(dataset_name, !!!source_type_lookup)
  ) %>%
  distinct(source_type, dataset_name, site, month_date, phenophase) %>%
  mutate(
    site = factor(site, levels = all_harmonised %>%
                    filter(!is.na(month_date)) %>%
                    group_by(site) %>%
                    summarise(first_obs = min(month_date), .groups = "drop") %>%
                    arrange(first_obs) %>%
                    pull(site))
  )

# Activity bands for FL/FR masting months
activity_band <- timeline_df %>%
  filter(phenophase %in% c("FL", "FR")) %>%
  distinct(site, month_date, phenophase)

# Build plot


sea_sites_temporal <- ggplot() +
  geom_line(
    data      = timeline_df,
    aes(x = month_date, y = site, colour = source_type),
    linewidth = 2.5,
    lineend   = "square",
    alpha     = 0.8
  ) +
  geom_tile(
    data   = activity_band,
    aes(x = month_date, y = site, fill = phenophase),
    height = 0.85,
    width  = 40,
    alpha  = 0.9
  ) +
  scale_colour_manual(
    name   = "Data source",
    values = source_cols
  ) +
  scale_fill_manual(
    name   = "Masting event",
    values = c("FL" = col_timeline_fl, "FR" = col_timeline_fr),
    labels = c("FL" = "Flowering", "FR" = "Fruiting")
  ) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y",
    expand      = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x = "Year",
    y = NULL
  ) +
  theme_masting() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 8),
    legend.position    = "right"
  ) +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 3, alpha = 1)),
    fill   = guide_legend(override.aes = list(size = 4, alpha = 0.8))
  )

sea_sites_temporal

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/all_datasets_timeline_figure.png",
  plot   = sea_sites_temporal,
  width  = fig_width,
  height = fig_height_double,
  dpi    = fig_dpi
)




#### ----------------------------------------------------------------- ####
####          H1: SEASONALITY OF MASTING ONSET                         ####
#### ----------------------------------------------------------------- ####

#### Step 6 - Proportion of observed site-months masting by calendar month ####
# Divides masting site-months by total observed site-months per calendar month.
# Controls for uneven sampling effort across months.
# A site-month is "observed" if phenophase is FL, FR, or none (not NA).

observed_by_month <- all_harmonised %>%
  filter(phenophase %in% c("FL", "FR", "none")) %>%
  mutate(calendar_month = month(month_date, label = TRUE, abbr = TRUE)) %>%
  group_by(calendar_month, event_type) %>%
  summarise(n_observed = n(), .groups = "drop") %>%
  filter(event_type %in% c("FL", "FR")) %>%
  rename(phenophase = event_type)

masting_by_month <- all_harmonised %>%
  filter(phenophase %in% c("FL", "FR")) %>%
  mutate(calendar_month = month(month_date, label = TRUE, abbr = TRUE)) %>%
  dplyr::count(calendar_month, phenophase, name = "n_masting")

prop_counts <- masting_by_month %>%
  left_join(observed_by_month, by = c("calendar_month", "phenophase")) %>%
  mutate(
    prop_masting = n_masting / n_observed,
    month_num    = as.numeric(calendar_month)
  )

#### Step 7 - Binomial significance test for proportion by calendar month ####
# Tests whether each month's proportion of masting site-months is
# significantly above the overall baseline proportion.

baseline <- all_harmonised %>%
  filter(phenophase %in% c("FL", "FR", "none")) %>%
  group_by(event_type) %>%
  summarise(
    n_masting     = sum(phenophase %in% c("FL", "FR")),
    n_observed    = n(),
    baseline_prop = n_masting / n_observed,
    .groups       = "drop"
  ) %>%
  filter(event_type %in% c("FL", "FR")) %>%
  rename(phenophase = event_type)

prop_counts_sig <- prop_counts %>%
  left_join(
    baseline %>%
      rename(n_total = n_observed, n_total_masting = n_masting),
    by = "phenophase"
  ) %>%
  rowwise() %>%
  mutate(
    binom_p     = binom.test(
      x           = n_masting,
      n           = n_observed,
      p           = baseline_prop,
      alternative = "greater"
    )$p.value,
    significant = binom_p < 0.05
  ) %>%
  ungroup() %>%
  mutate(
    fill_group = case_when(
      phenophase == "FL" & significant  ~ "FL_sig",
      phenophase == "FL" & !significant ~ "FL_nonsig",
      phenophase == "FR" & significant  ~ "FR_sig",
      phenophase == "FR" & !significant ~ "FR_nonsig"
    ),
    border_col = case_when(
      phenophase == "FL" & significant  ~ col_FL_border,
      phenophase == "FL" & !significant ~ col_FL_nonsig,
      phenophase == "FR" & significant  ~ col_FR_border,
      phenophase == "FR" & !significant ~ col_FR_nonsig
    )
  )

# Print results
prop_counts_sig %>%
  dplyr::select(calendar_month, phenophase, prop_masting,
                baseline_prop, binom_p, significant) %>%
  arrange(phenophase, calendar_month) %>%
  print(n = 30)

#### Step 8/9 - Plots: proportion masting by calendar month (line, then bar) ####

#### Step 9 - Plot: proportion masting by calendar month (line version) ####
# Alternative version requested by advisor — solid lines connecting monthly
# proportions, no bars or trend lines.

p_prop_line <- ggplot(prop_counts,
                      aes(x     = month_num,
                          y     = prop_masting,
                          colour = phenophase,
                          group  = phenophase)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_x_continuous(
    breaks = 1:12,
    labels = month.abb
  ) +
  scale_colour_manual(
    name   = "Phenophase",
    values = phase_cols,
    labels = phase_labels
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "Month",
    y = "Proportion of observed site-months masting"
  ) +
  theme_masting()

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/proportion_masting_line.png",
  plot   = p_prop_line,
  width  = fig_width,
  height = fig_height_single,
  dpi    = fig_dpi
)

#### Proportional bar chart w/ significance ####

p_prop <- ggplot(prop_counts_sig,
                 aes(x = calendar_month, y = prop_masting)) +
  geom_col(
    aes(fill = fill_group, colour = fill_group),
    position = "dodge",
    width    = 0.8
  ) +
  scale_fill_manual(
    name   = NULL,
    values = c(
      "FL_sig"    = col_FL_sig,
      "FL_nonsig" = col_FL_nonsig,
      "FR_sig"    = col_FR_sig,
      "FR_nonsig" = col_FR_nonsig
    ),
    labels = c(
      "FL_sig"    = "Flowering (p < 0.05)",
      "FL_nonsig" = "Flowering (p ≥ 0.05)",
      "FR_sig"    = "Fruiting (p < 0.05)",
      "FR_nonsig" = "Fruiting (p ≥ 0.05)"
    )
  ) +
  scale_colour_manual(
    name   = NULL,
    values = c(
      "FL_sig"    = col_FL_border,
      "FL_nonsig" = col_FL_nonsig,
      "FR_sig"    = col_FR_border,
      "FR_nonsig" = col_FR_nonsig
    ),
    labels = c(
      "FL_sig"    = "Flowering (p < 0.05)",
      "FL_nonsig" = "Flowering (p ≥ 0.05)",
      "FR_sig"    = "Fruiting (p < 0.05)",
      "FR_nonsig" = "Fruiting (p ≥ 0.05)"
    )
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Month",
    y = "Proportion of observed site-months masting"
  ) +
  theme_masting()

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/proportion_masting_by_calendar_month_v2.png",
  plot   = p_prop,
  width  = fig_width,
  height = fig_height_single,
  dpi    = fig_dpi
)

#### ----------------------------------------------------------------- ####
####          H1: CIRCULAR ONSET MONTH ANALYSIS                        ####
#### ----------------------------------------------------------------- ####

#### Step 10 - Build all_onsets_full ####
# One row per masting onset event (FL or FR) per site.
# Long-term monitoring and expert obs: onset detected from is_fl_onset /
# is_fr_onset flags in all_harmonised.
# Literature review: onset taken directly from pheno which has explicit
# onset_date per event — more reliable than consecutive month detection.

# FL onsets — long-term monitoring and expert obs
fl_onsets_ltm <- all_harmonised %>%
  filter(
    dataset_name %in% c("FRIM", "Lambir", "Gunung_Palung", "Pasoh",
                        "Bukit_Timah", "additional_contributed_obs",
                        "BBBR", "Gunung_Tarak"),
    is_fl_onset == TRUE
  ) %>%
  mutate(phase = "FL") %>%
  dplyr::select(site, onset_date = month_date, phase, dataset_name, area)

# FR onsets — long-term monitoring and expert obs
fr_onsets_ltm <- all_harmonised %>%
  filter(
    dataset_name %in% c("FRIM", "Lambir", "Gunung_Palung", "Pasoh",
                        "Bukit_Timah", "additional_contributed_obs",
                        "BBBR", "Gunung_Tarak"),
    is_fr_onset == TRUE
  ) %>%
  mutate(phase = "FR") %>%
  dplyr::select(site, onset_date = month_date, phase, dataset_name, area)

# FL onsets — literature review (from pheno directly)
fl_onsets_lit <- pheno %>%
  filter(
    fruit_flower == "flower",
    !covidence_id %in% exclude_cov_ids,
    !site %in% excluded_sites
  ) %>%
  mutate(
    site         = dplyr::recode(site, !!!pheno_site_recode),
    phase        = "FL",
    dataset_name = "Literature_review"
  ) %>%
  dplyr::select(site, onset_date, phase, dataset_name, area)

# FR onsets — literature review (from pheno directly)
fr_onsets_lit <- pheno %>%
  filter(
    fruit_flower == "fruit",
    !covidence_id %in% exclude_cov_ids,
    !site %in% excluded_sites
  ) %>%
  mutate(
    site         = dplyr::recode(site, !!!pheno_site_recode),
    phase        = "FR",
    dataset_name = "Literature_review"
  ) %>%
  dplyr::select(site, onset_date, phase, dataset_name, area)

# Combine all onsets
# onset_month and year calculated here consistently for all datasets
all_onsets_full <- bind_rows(
  fl_onsets_ltm,
  fr_onsets_ltm,
  fl_onsets_lit,
  fr_onsets_lit
) %>%
  mutate(
    onset_month = month(onset_date),
    year        = year(onset_date)
  ) %>%
  arrange(year)

# Sanity checks
all_onsets_full %>% dplyr::count(phase, dataset_name)
all_onsets_full %>% dplyr::count(area)
range(all_onsets_full$year)


#### Step 11 - Rayleigh tests for onset month seasonality ####
# Tests whether FL and FR onset months are non-uniformly distributed
# around the calendar year.

# FL — all sites
fl_circular <- circular(
  all_onsets_full %>% filter(phase == "FL") %>%
    pull(onset_month) * (2 * pi / 12),
  units = "radians", template = "none"
)
rayleigh_fl   <- rayleigh.test(fl_circular)
mean_fl       <- mean(fl_circular)
mean_month_fl <- (mean_fl * 12 / (2 * pi)) %% 12

# FR — all sites
fr_circular <- circular(
  all_onsets_full %>% filter(phase == "FR") %>%
    pull(onset_month) * (2 * pi / 12),
  units = "radians", template = "none"
)
rayleigh_fr   <- rayleigh.test(fr_circular)
mean_fr       <- mean(fr_circular)
mean_month_fr <- (mean_fr * 12 / (2 * pi)) %% 12

# Summary table
tibble::tibble(
  phase      = c("FL", "FR"),
  n          = c(nrow(all_onsets_full %>% filter(phase == "FL")),
                 nrow(all_onsets_full %>% filter(phase == "FR"))),
  r          = c(round(rayleigh_fl$statistic, 3),
                 round(rayleigh_fr$statistic, 3)),
  p_value    = c(rayleigh_fl$p.value, rayleigh_fr$p.value),
  mean_onset = c(month.abb[round(mean_month_fl)],
                 month.abb[round(mean_month_fr)])
) %>% print()

#### Step 12 - Rayleigh tests by subregion ####

rayleigh_by_subregion <- all_onsets_full %>%
  mutate(subregion = dplyr::recode(area, !!!subregion_recode)) %>%
  filter(subregion %in% c("Peninsular Malaysia & Singapore", "Borneo")) %>%
  group_by(subregion, phase) %>%
  summarise(
    n       = n(),
    r       = rayleigh.test(
      circular(onset_month * (2 * pi / 12),
               units = "radians", template = "none")
    )$statistic,
    p_value = rayleigh.test(
      circular(onset_month * (2 * pi / 12),
               units = "radians", template = "none")
    )$p.value,
    mean_month = {
      circ <- circular(onset_month * (2 * pi / 12),
                       units = "radians", template = "none")
      (mean(circ) * 12 / (2 * pi)) %% 12
    },
    .groups = "drop"
  ) %>%
  mutate(mean_month_name = month.abb[round(mean_month)])

print(rayleigh_by_subregion)


#### Step 13 - Three-panel circular onset plot ####
# Panel A: all sites; Panel B: Borneo; Panel C: Peninsular Malaysia & Singapore
# Uses patchwork for layout, subregion-specific Rayleigh stats in subtitles.

# Helper function to build one circular plot
# (final theme_masting() version — an earlier draft used theme_minimal()
# with manual title/subtitle sizing, but this house-style version is what
# actually feeds the saved figure)
make_circ_plot <- function(data, title, subtitle) {
  counts <- data %>%
    dplyr::count(onset_month, phase) %>%
    mutate(
      month_label = factor(month.abb[onset_month], levels = month.abb),
      phase       = factor(phase, levels = c("FL", "FR"))
    )
  ggplot(counts, aes(x = month_label, y = n, fill = phase)) +
    geom_col(
      position = "dodge",
      width    = 0.8,
      alpha    = 0.85,
      colour   = "grey30"
    ) +
    coord_polar(start = -pi/12) +
    scale_fill_manual(
      name   = "Phenophase",
      values = phase_cols,
      labels = phase_labels
    ) +
    scale_y_continuous(breaks = seq(0, 20, 5)) +
    labs(
      title    = title,
      subtitle = subtitle,
      x        = NULL,
      y        = "Number of onset events"
    ) +
    theme_masting() +
    theme(
      axis.text.y     = element_blank(),
      panel.grid      = element_line(colour = "grey90"),
      legend.position = "none"
    )
}

# Helper function to build subregion-specific subtitle
make_subtitle <- function(subregion_name, rayleigh_df) {
  fl <- rayleigh_df %>% filter(subregion == subregion_name, phase == "FL")
  fr <- rayleigh_df %>% filter(subregion == subregion_name, phase == "FR")
  paste0(
    "FL: r = ", round(fl$r, 2),
    ", p = ", format(fl$p_value, scientific = TRUE, digits = 2),
    ", mean = ", fl$mean_month_name, "\n",
    "FR: r = ", round(fr$r, 2),
    ", p = ", format(fr$p_value, scientific = TRUE, digits = 2),
    ", mean = ", fr$mean_month_name
  )
}

# Add subregion column to all_onsets_full
all_onsets_subregion <- all_onsets_full %>%
  mutate(subregion = dplyr::recode(area, !!!subregion_recode))

# Panel A — all sites
p_all <- make_circ_plot(
  data     = all_onsets_full,
  title    = "All sites",
  subtitle = paste0(
    "FL: r = ", round(rayleigh_fl$statistic, 2),
    ", p = ", format(rayleigh_fl$p.value, scientific = TRUE, digits = 2),
    ", mean = ", month.abb[round(mean_month_fl)], "\n",
    "FR: r = ", round(rayleigh_fr$statistic, 2),
    ", p = ", format(rayleigh_fr$p.value, scientific = TRUE, digits = 2),
    ", mean = ", month.abb[round(mean_month_fr)]
  )
)

# Panel B — Borneo
p_borneo <- make_circ_plot(
  data     = all_onsets_subregion %>% filter(subregion == "Borneo"),
  title    = "Borneo",
  subtitle = make_subtitle("Borneo", rayleigh_by_subregion)
)

# Panel C — Peninsular Malaysia & Singapore
p_pennmal <- make_circ_plot(
  data     = all_onsets_subregion %>%
    filter(subregion == "Peninsular Malaysia & Singapore"),
  title    = "Peninsular Malaysia & Singapore",
  subtitle = make_subtitle("Peninsular Malaysia & Singapore",
                           rayleigh_by_subregion)
)

# Add legend to panel A, combine with patchwork
p_all_legend <- p_all +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(
    title        = "Phenophase",
    override.aes = list(size = 5)
  ))

combined_circ_final <- (p_all_legend | p_borneo | p_pennmal) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(plot.tag = element_text(size = 11, face = "bold"))

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/onset_circular_three_panel.png",
  plot   = combined_circ_final,
  width  = fig_width,
  height = fig_height_double,
  dpi    = fig_dpi
)





###### testing if onset month is further from the mean month closer to the equator ######

#### Latitude vs onset month deviation from regional mean ####
# Uses circular (folded) deviation: max possible distance from the seasonal
# mean month is 6 months. Regression is unweighted (weighting by n_events
# made CIs unrealistically narrow without changing significance - MSL).

# Site coordinates lookup — needed here and again later for the onset map
site_coords <- all_harmonised %>%
  filter(!is.na(lat), !is.na(long)) %>%
  group_by(site) %>%
  summarise(
    long = median(long, na.rm = TRUE),
    lat  = median(lat,  na.rm = TRUE),
    .groups = "drop"
  )

lat_deviation <- all_onsets_full %>%
  left_join(site_coords, by = "site") %>%
  mutate(
    mean_month = if_else(phase == "FL", 4, 7),
    deviation  = pmin(abs(onset_month - mean_month),
                      12 - abs(onset_month - mean_month))
  ) %>%
  filter(!is.na(lat))

# Per-site summary
lat_deviation_summary <- lat_deviation %>%
  group_by(site, phase, lat) %>%
  summarise(
    mean_deviation = mean(deviation, na.rm = TRUE),
    n_events       = n(),
    .groups        = "drop"
  ) %>%
  mutate(phase_label = if_else(phase == "FL", "Flowering", "Fruiting"))

# Linear models (unweighted) — one row per site, so no random effect needed
lm_fl <- lm(mean_deviation ~ lat,
            data = lat_deviation_summary %>% filter(phase_label == "Flowering"))
lm_fr <- lm(mean_deviation ~ lat,
            data = lat_deviation_summary %>% filter(phase_label == "Fruiting"))

summary(lm_fl)  # Flowering: R2 = 0.37, p = 0.006
summary(lm_fr)  # Fruiting:  R2 = 0.08, p = 0.192

# Annotation labels for each facet panel
stat_labels <- data.frame(
  phase_label = c("Flowering", "Fruiting"),
  label = c(
    paste0("R² = ", round(summary(lm_fl)$r.squared, 2),
           ", p = ", round(summary(lm_fl)$coefficients[2,4], 3)),
    paste0("R² = ", round(summary(lm_fr)$r.squared, 2),
           ", p = ", round(summary(lm_fr)$coefficients[2,4], 3))
  )
)

# Colors matched to timeline figure palette
col_timeline_fl <- "#916392"  # light purple — Flowering
col_timeline_fr <- "#371B38"  # dark purple — Fruiting
phase_cols_timeline <- c("FL" = col_timeline_fl, "FR" = col_timeline_fr)

sites_to_label <- c(
  "FRIM Malaysia - Kepong - Ampang", "Lambir Hills",
  "Singapore", "Danum Valley", "Barito Ulu",
  "BBBR", "Gunung Palung National Park",
  "Ketambe", "Royal Belum State Park",
  "Ulu Segama Malua", "Maliau Basin", "Gombak",
  "Pasoh Forest Reserve", "Gunung Tarak",
  "Deramakot Forest Reserve", "Sungai Wain Protection Forest",
  "Tangkulap Forest Reserve", "Anduki Forest Reserve", "Badas"
)

lat_deviation_plot <- lat_deviation_summary %>%
  ggplot(aes(x = lat, y = mean_deviation, colour = phase)) +
  geom_point(aes(size = n_events), alpha = 0.8) +
  geom_smooth(
    method      = "lm",
    se          = TRUE,
    linewidth   = 1.2,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data         = lat_deviation_summary %>%
      filter(site %in% sites_to_label),
    aes(label = site),
    size         = 2.5,
    colour       = "grey30",
    show.legend  = FALSE,
    max.overlaps = Inf
  ) +
  geom_text(
    data  = stat_labels,
    aes(label = label),
    x     = -Inf, y = -Inf,
    hjust = -0.1, vjust = -0.5,
    size  = 3,
    colour = "grey30",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ phase_label) +
  scale_colour_manual(values = phase_cols_timeline, guide = "none") +
  scale_size_continuous(
    name   = "Number of\nonset events",
    range  = c(2, 8),
    breaks = c(2, 5, 10, 20)
  ) +
  scale_x_continuous(
    name   = "Latitude (degrees from equator)",
    breaks = seq(-3, 6, by = 2)
  ) +
  scale_y_continuous(
    name   = "Mean onset deviation from regional mean (months)",
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  theme_masting() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11,
                                    margin = margin(b = 8)),
    legend.position  = "top",
    legend.direction = "horizontal"
  ) +
  guides(size = guide_legend(title = "Number of\nonset events", nrow = 1))

lat_deviation_plot

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/latitude_deviation_plot.png",
  plot   = lat_deviation_plot,
  width  = fig_width,
  height = fig_height_single,
  dpi    = fig_dpi
)


# NOTE FOR MSL: your deviation formula in the earlier draft had a bug —
# deviation = pmin(onset_month - mean_month), 12 - (onset_month - mean_month)
# The comma after pmin(...) closes it with only ONE argument, so pmin() just
# returns onset_month - mean_month unchanged (no abs(), no actual folding).
# "12 - (onset_month - mean_month)" becomes a separate, unused mutate() column
# instead of the second pmin() comparison. Result: deviation was an unbounded
# signed difference (~-8 to +8), not a proper circular distance (max 6 months).
# This likely explains why the original figure showed R2 = 0.37/0.08 instead
# of the corrected R2 = 0.62/0.28 - the noisier, unfolded deviation was
# washing out some of the real lat-onset relationship.
# Fixed version: deviation = pmin(abs(onset_month - mean_month),
#                                 12 - abs(onset_month - mean_month))

#### ----------------------------------------------------------------- ####
####          H2: MASTING FREQUENCY AND DURATION BY SITE               ####
#### ----------------------------------------------------------------- ####

#### Step 14 - Define sites for frequency analysis ####
# Restricted to sites with effectively continuous monitoring coverage.
# FL analysis: all sites listed.
# FR analysis: all sites except Pasoh (awaiting FR data).

fl_sites <- c(
  "FRIM Malaysia - Kepong - Ampang",
  "Singapore",
  "Lambir Hills",
  "Gunung Palung National Park",
  "Pasoh Forest Reserve",
  "Danum Valley",
  "Ulu Segama Malua",
  "Maliau Basin",
  "Barito Ulu",
  "Gombak"
)

fr_sites <- fl_sites[fl_sites != "Pasoh Forest Reserve"]


#### Step 15 - Observation period per site ####
# Used as denominator for frequency calculation.

obs_period <- all_harmonised %>%
  filter(site %in% fl_sites) %>%
  group_by(site) %>%
  summarise(
    obs_start  = min(month_date),
    obs_end    = max(month_date),
    obs_years  = as.numeric(interval(min(month_date),
                                     max(month_date)) / years(1)),
    obs_months = n_distinct(month_date),
    .groups    = "drop"
  )


#### Step 16 - FL frequency and duration ####

# FL events from long-term monitoring and expert obs
# Uses is_fl_onset flags and cumsum to identify discrete events
fl_events_ltm <- all_harmonised %>%
  filter(
    site %in% fl_sites,
    dataset_name %in% c("FRIM", "Lambir", "Gunung_Palung", "Pasoh",
                        "Bukit_Timah", "additional_contributed_obs",
                        "BBBR", "Gunung_Tarak"),
    event_type == "FL",
    phenophase == "FL"
  ) %>%
  arrange(site, month_date) %>%
  group_by(site) %>%
  mutate(
    is_onset = is_fl_onset == TRUE,
    run_id   = cumsum(is_onset)
  ) %>%
  group_by(site, run_id) %>%
  summarise(
    onset_date      = min(month_date),
    end_date        = max(month_date),
    duration_months = interval(min(month_date),
                               max(month_date)) %/% months(1) + 1,
    dataset_name    = first(dataset_name),
    .groups         = "drop"
  ) %>%
  dplyr::select(site, onset_date, end_date, duration_months, dataset_name)

# FL events from literature review sites (Barito Ulu and Gombak)
# Uses explicit onset/end dates from pheno
fl_events_lit <- pheno %>%
  filter(
    fruit_flower == "flower",
    site %in% c("Barito Ulu", "Gombak"),
    !covidence_id %in% exclude_cov_ids
  ) %>%
  mutate(
    duration_months = interval(onset_date, end_date) %/% months(1) + 1,
    dataset_name    = "Literature_review"
  ) %>%
  dplyr::select(site, onset_date, end_date, duration_months, dataset_name)

fl_events <- bind_rows(fl_events_ltm, fl_events_lit)

fl_summary <- fl_events %>%
  left_join(obs_period, by = "site") %>%
  group_by(site, obs_years) %>%
  summarise(
    n_fl_events      = n(),
    mean_fl_duration = mean(duration_months),
    sd_fl_duration   = sd(duration_months),
    min_fl_duration  = min(duration_months),
    max_fl_duration  = max(duration_months),
    .groups          = "drop"
  ) %>%
  mutate(fl_frequency = n_fl_events / obs_years) %>%
  arrange(desc(obs_years))

print(fl_summary, width = Inf)


#### Step 17 - FR frequency and duration ####

# FR events from long-term monitoring and expert obs
fr_events_ltm <- all_harmonised %>%
  filter(
    site %in% fr_sites,
    dataset_name %in% c("FRIM", "Lambir", "Gunung_Palung",
                        "Bukit_Timah", "additional_contributed_obs",
                        "BBBR", "Gunung_Tarak"),
    event_type == "FR",
    phenophase == "FR"
  ) %>%
  arrange(site, month_date) %>%
  group_by(site) %>%
  mutate(run_id = cumsum(is_fr_onset)) %>%
  group_by(site, run_id) %>%
  summarise(
    onset_date      = min(month_date),
    end_date        = max(month_date),
    duration_months = interval(min(month_date),
                               max(month_date)) %/% months(1) + 1,
    dataset_name    = first(dataset_name),
    .groups         = "drop"
  ) %>%
  filter(run_id > 0) %>%
  dplyr::select(site, onset_date, end_date, duration_months, dataset_name)

# FR events from literature review sites
fr_events_lit <- pheno %>%
  filter(
    fruit_flower == "fruit",
    site %in% c("Barito Ulu", "Gombak"),
    !covidence_id %in% exclude_cov_ids
  ) %>%
  mutate(
    duration_months = interval(onset_date, end_date) %/% months(1) + 1,
    dataset_name    = "Literature_review"
  ) %>%
  dplyr::select(site, onset_date, end_date, duration_months, dataset_name)

fr_events <- bind_rows(fr_events_ltm, fr_events_lit)

fr_summary <- fr_events %>%
  left_join(obs_period, by = "site") %>%
  group_by(site, obs_years) %>%
  summarise(
    n_fr_events      = n(),
    mean_fr_duration = mean(duration_months),
    sd_fr_duration   = sd(duration_months),
    min_fr_duration  = min(duration_months),
    max_fr_duration  = max(duration_months),
    .groups          = "drop"
  ) %>%
  mutate(fr_frequency = n_fr_events / obs_years) %>%
  arrange(desc(obs_years))

print(fr_summary, width = Inf)


#### Step 18 - Combined frequency and duration summary table ####

freq_duration_summary <- fl_summary %>%
  dplyr::select(site, obs_years,
                n_fl_events, fl_frequency,
                mean_fl_duration, sd_fl_duration) %>%
  full_join(
    fr_summary %>%
      dplyr::select(site,
                    n_fr_events, fr_frequency,
                    mean_fr_duration, sd_fr_duration),
    by = "site"
  ) %>%
  arrange(desc(obs_years))

print(freq_duration_summary)


#### Step 19 - Bubble chart: frequency vs duration by site ####
# Bubble size = years of monitoring; fill = subregion; shape = phenophase.

# Site groupings and short labels
site_subregion <- c(
  "FRIM Malaysia - Kepong - Ampang" = "Peninsular Malaysia & Singapore",
  "Singapore"                       = "Peninsular Malaysia & Singapore",
  "Pasoh Forest Reserve"            = "Peninsular Malaysia & Singapore",
  "Gombak"                          = "Peninsular Malaysia & Singapore",
  "Lambir Hills"                    = "Borneo",
  "Danum Valley"                    = "Borneo",
  "Ulu Segama Malua"                = "Borneo",
  "Maliau Basin"                    = "Borneo",
  "Barito Ulu"                      = "Borneo",
  "Gunung Palung National Park"     = "Borneo"
)

site_labels <- c(
  "FRIM Malaysia - Kepong - Ampang" = "FRIM",
  "Gunung Palung National Park"     = "Gunung Palung",
  "Pasoh Forest Reserve"            = "Pasoh",
  "Ulu Segama Malua"                = "Ulu Segama"
)

bubble_data <- bind_rows(
  fl_summary %>%
    mutate(phase = "FL") %>%
    dplyr::select(site, obs_years, phase,
                  frequency = fl_frequency,
                  duration  = mean_fl_duration),
  fr_summary %>%
    mutate(phase = "FR") %>%
    dplyr::select(site, obs_years, phase,
                  frequency = fr_frequency,
                  duration  = mean_fr_duration)
) %>%
  mutate(
    subregion  = dplyr::recode(site, !!!site_subregion),
    site_label = dplyr::recode(site, !!!site_labels, .default = site)
  )



#### new format and save

bubble_plot <- ggplot(
  bubble_data %>% filter(!is.na(frequency), !is.na(duration)),
  aes(x     = frequency,
      y     = duration,
      size  = obs_years,
      fill  = subregion,
      shape = phase)
) +
  geom_point(alpha = 0.85, colour = "grey30") +
  ggrepel::geom_text_repel(
    aes(label = site_label),
    size               = 3,
    max.overlaps       = Inf,
    show.legend        = FALSE,
    box.padding        = 1.2,
    point.padding      = 0.5,
    min.segment.length = 0.3,
    segment.colour     = "grey50",
    segment.size       = 0.3
  ) +
  scale_size_continuous(
    name   = "Years monitored",
    range  = c(3, 12),
    breaks = c(10, 20, 30, 40),
    limits = c(10, 45)
  ) +
  scale_fill_manual(
    name   = "Subregion",
    values = c(
      "Peninsular Malaysia & Singapore" = "#7FBEAB",
      "Borneo"                          = "#2C365E"
    ),
    labels = c(
      "Peninsular Malaysia & Singapore" = "Peninsular Malaysia\n& Singapore",
      "Borneo"                          = "Borneo"
    ),
    guide  = guide_legend(override.aes = list(shape = 21, size = 5))
  ) +
  scale_shape_manual(
    name   = "Phenophase",
    values = c("FL" = 21, "FR" = 23),
    labels = c("FL" = "Flowering", "FR" = "Fruiting"),
    guide  = guide_legend(override.aes = list(size = 5, fill = "grey50"))
  ) +
  scale_x_continuous(
    name   = "Masting frequency (events per year)",
    limits = c(0, NA),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  scale_y_continuous(
    name   = "Mean mast duration (months)",
    limits = c(0, NA),
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  theme_masting()

bubble_plot

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/masting_frequency_duration_bubble.png",
  plot   = bubble_plot,
  width  = fig_width,
  height = fig_height_single + 1,
  dpi    = fig_dpi
)


### getting summaries for results text #####
obs_period %>%
  arrange(desc(obs_years)) %>%
  print(width = Inf)

freq_duration_summary %>%
  summarise(
    fl_freq_min  = min(fl_frequency, na.rm = TRUE),
    fl_freq_max  = max(fl_frequency, na.rm = TRUE),
    fr_freq_min  = min(fr_frequency, na.rm = TRUE),
    fr_freq_max  = max(fr_frequency, na.rm = TRUE),
    fl_dur_min   = min(mean_fl_duration, na.rm = TRUE),
    fl_dur_max   = max(mean_fl_duration, na.rm = TRUE),
    fr_dur_min   = min(mean_fr_duration, na.rm = TRUE),
    fr_dur_max   = max(mean_fr_duration, na.rm = TRUE)
  )

freq_duration_summary %>%
  summarise(fr_dur_max = max(mean_fr_duration, na.rm = TRUE))


# Correlation between frequency and duration
# Run separately for FL and FR
cor_fl <- cor.test(fl_summary$fl_frequency, fl_summary$mean_fl_duration,
                   method = "spearman")
cor_fr <- cor.test(fr_summary$fr_frequency, fr_summary$mean_fr_duration,
                   method = "spearman")

cat("FL frequency vs duration: rho =", round(cor_fl$estimate, 3),
    ", p =", round(cor_fl$p.value, 3), "\n")
cat("FR frequency vs duration: rho =", round(cor_fr$estimate, 3),
    ", p =", round(cor_fr$p.value, 3), "\n")

#### ----------------------------------------------------------------- ####
####     H3: TEMPORAL CLUSTERING OF MASTING EVENTS ACROSS SITES        ####
#### ----------------------------------------------------------------- ####

#### Step 20 - Define eligible sites and build site-month matrix ####
# Sites included: long-term monitoring + expert obs +
# literature review sites with ≥3 years coverage and none periods.
# For sites in both monitoring and literature, monitoring takes priority.

sync_sites <- c(
  "Gunung Palung National Park",
  "FRIM Malaysia - Kepong - Ampang",
  "Singapore",
  "Danum Valley",
  "Ulu Segama Malua",
  "Lambir Hills",
  "Pasoh Forest Reserve",
  "Maliau Basin",
  "Barito Ulu",
  "Gombak",
  "Sungai Wain Protection Forest",
  "Deramakot Forest Reserve",
  "Lanjak Entimau Wildlife Sanctuary",
  "Tangkulap Forest Reserve",
  "Ketambe"
)

site_dataset_priority <- all_harmonised %>%
  filter(site %in% sync_sites) %>%
  mutate(
    priority = case_when(
      dataset_name %in% c("FRIM", "Lambir", "Gunung_Palung",
                          "Pasoh", "Bukit_Timah",
                          "additional_contributed_obs",
                          "BBBR", "Gunung_Tarak")  ~ 1,
      dataset_name == "Literature_review"           ~ 2
    )
  ) %>%
  group_by(site) %>%
  filter(priority == min(priority)) %>%
  ungroup()

# Observation window per site
site_windows <- site_dataset_priority %>%
  group_by(site) %>%
  summarise(
    obs_start = min(month_date),
    obs_end   = max(month_date),
    .groups   = "drop"
  )

# Masting indicator per site per month
# TRUE = masting (FL or FR), FALSE = observed not masting
site_month_matrix <- site_dataset_priority %>%
  filter(phenophase %in% c("FL", "FR", "none")) %>%
  group_by(site, month_date) %>%
  summarise(
    is_masting = any(phenophase %in% c("FL", "FR")),
    .groups    = "drop"
  )


#### Step 21 - Build model data ####
# One row per site-month.
# Response: is_masting (0/1).
# Predictor: n_co_masting = number of OTHER sites masting within ±3 months.
# Note: this step is slow due to rowwise calculation — allow a few minutes.

proximity_window <- 3

model_data <- site_month_matrix %>%
  filter(site %in% sync_sites) %>%
  mutate(is_masting = as.integer(is_masting)) %>%
  rowwise() %>%
  mutate(
    n_co_masting = {
      focal_site  <- site
      focal_month <- month_date
      site_month_matrix %>%
        filter(
          site       != focal_site,
          site       %in% sync_sites,
          month_date >= focal_month - months(proximity_window),
          month_date <= focal_month + months(proximity_window),
          is_masting == TRUE
        ) %>%
        distinct(site) %>%
        nrow()
    }
  ) %>%
  ungroup()

cat("Model data rows:", nrow(model_data), "\n")
cat("Masting months:", sum(model_data$is_masting), "\n")
cat("Non-masting months:", sum(model_data$is_masting == 0), "\n")


#### Step 22 - Fit mixed effects logistic regression ####
# Fixed effect: n_co_masting (temporal clustering predictor)
# Random effect: site (controls for baseline masting frequency per site)

glmer_model <- glmer(
  is_masting ~ n_co_masting + (1 | site),
  data   = model_data,
  family = binomial
)

summary(glmer_model)

# Coefficients and odds ratios
glmer_results <- broom.mixed::tidy(
  glmer_model,
  conf.int     = TRUE,
  exponentiate = TRUE
) %>%
  filter(effect == "fixed")

print(glmer_results)

# Compare to null model (intercept + random effect only)
null_model <- glmer(
  is_masting ~ 1 + (1 | site),
  data   = model_data,
  family = binomial
)

cat("\nModel comparison (null vs full):\n")
anova(null_model, glmer_model)


#### Step 23 - Visualise logistic regression results ####

# Predicted probability curve with confidence ribbon
pred_data <- expand.grid(n_co_masting = seq(0, 7, by = 0.1)) %>%
  mutate(
    lp        = fixef(glmer_model)[1] +
      fixef(glmer_model)[2] * n_co_masting,
    lp_low    = (fixef(glmer_model)[1] -
                   1.96 * sqrt(vcov(glmer_model)[1,1])) +
      (fixef(glmer_model)[2] -
         1.96 * sqrt(vcov(glmer_model)[2,2])) * n_co_masting,
    lp_high   = (fixef(glmer_model)[1] +
                   1.96 * sqrt(vcov(glmer_model)[1,1])) +
      (fixef(glmer_model)[2] +
         1.96 * sqrt(vcov(glmer_model)[2,2])) * n_co_masting,
    prob      = plogis(lp),
    prob_low  = plogis(lp_low),
    prob_high = plogis(lp_high)
  )
obs_props <- model_data %>%
  group_by(n_co_masting) %>%
  summarise(
    n_total   = n(),
    n_masting = sum(is_masting),
    prop      = n_masting / n_total,
    .groups   = "drop"
  )

pred_integer <- pred_data %>%
  filter(n_co_masting == round(n_co_masting), n_co_masting <= 7)

#### new format and theme for save

p1 <- ggplot(pred_data, aes(x = n_co_masting, y = prob)) +
  geom_ribbon(
    aes(ymin = prob_low, ymax = prob_high),
    fill  = col_FL,
    alpha = 0.2
  ) +
  geom_line(colour = col_FL, linewidth = 1.2) +
  scale_x_continuous(
    breaks = 0:7,
    name   = "Number of other sites masting within ±3 months"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name   = "Predicted probability of masting",
    limits = c(0, NA)
  ) +
  annotate(
    "text",
    x     = 1,
    y     = 0.75,
    hjust = 0,
    label = paste0(
      "OR = ", round(exp(fixef(glmer_model)[2]), 2),
      " (95% CI: ",
      round(glmer_results$conf.low[2], 2), "–",
      round(glmer_results$conf.high[2], 2),
      "), p < 0.001"
    ),
    size   = 3,
    colour = "grey30"
  ) +
  theme_masting()

p2 <- ggplot() +
  geom_col(
    data  = obs_props,
    aes(x = n_co_masting, y = prop),
    fill  = col_FL,
    alpha = 0.4,
    width = 0.6
  ) +
  geom_line(
    data      = pred_integer,
    aes(x = n_co_masting, y = prob),
    colour    = col_FL,
    linewidth = 1.2
  ) +
  geom_point(
    data   = pred_integer,
    aes(x = n_co_masting, y = prob),
    colour = col_FL,
    size   = 2.5
  ) +
  geom_text(
    data   = obs_props,
    aes(x = n_co_masting, y = prop + 0.01,
        label = paste0("n=", n_total)),
    size   = 3, colour = "grey40", vjust = 0
  ) +
  scale_x_continuous(
    breaks = 0:7,
    name   = "Number of other sites masting within ±3 months"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name   = "Observed proportion masting",
    limits = c(0, NA)
  ) +
  theme_masting()

### fixing logisitc plot and stacking vertically ###

# Shared x-axis label using patchwork
logistic_plot <- (p1 / p2) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.tag    = element_text(size = 11, face = "bold"),
    plot.margin = margin(5, 10, 5, 5)
  )

logistic_plot

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/temporal_clustering_logistic.png",
  plot   = logistic_plot,
  width  = fig_width,
  height = fig_height_double,
  dpi    = fig_dpi
)

##### SENSITIVITY ANALYSIS ######


#### Sensitivity analysis — site inclusion threshold ####
# Tests whether the 3-year minimum coverage threshold for site inclusion
# affects the clustering result. Compares 2-year, 3-year, and 5-year thresholds.

run_clustering_model <- function(min_years) {
  
  cat("\n#### Running model with min_years =", min_years, "####\n")
  
  # Identify eligible sites at this threshold
  eligible_sites <- all_harmonised %>%
    filter(site %in% sync_sites) %>%
    group_by(site) %>%
    summarise(
      span_years = as.numeric(interval(min(month_date),
                                       max(month_date)) / years(1)),
      has_none   = any(phenophase == "none"),
      .groups    = "drop"
    ) %>%
    filter(span_years >= min_years, has_none == TRUE) %>%
    pull(site)
  
  cat("Sites included:", length(eligible_sites), "\n")
  cat("Sites:", paste(eligible_sites, collapse = ", "), "\n")
  
  # Build site-month matrix for these sites
  site_data <- site_dataset_priority %>%
    filter(site %in% eligible_sites)
  
  site_matrix <- site_data %>%
    filter(phenophase %in% c("FL", "FR", "none")) %>%
    group_by(site, month_date) %>%
    summarise(
      is_masting = any(phenophase %in% c("FL", "FR")),
      .groups    = "drop"
    )
  
  # Build model data with co-masting predictor
  model_dat <- site_matrix %>%
    mutate(is_masting = as.integer(is_masting)) %>%
    rowwise() %>%
    mutate(
      n_co_masting = {
        focal_site  <- site
        focal_month <- month_date
        site_matrix %>%
          filter(
            site       != focal_site,
            month_date >= focal_month - months(3),
            month_date <= focal_month + months(3),
            is_masting == TRUE
          ) %>%
          distinct(site) %>%
          nrow()
      }
    ) %>%
    ungroup()
  
  # Fit models
  full_model <- glmer(
    is_masting ~ n_co_masting + (1 | site),
    data   = model_dat,
    family = binomial
  )
  
  null_model <- glmer(
    is_masting ~ 1 + (1 | site),
    data   = model_dat,
    family = binomial
  )
  
  # Extract results
  or        <- exp(fixef(full_model)["n_co_masting"])
  ci        <- exp(confint(full_model, parm = "beta_", method = "Wald"))
  lrt       <- anova(null_model, full_model)
  
  tibble::tibble(
    min_years    = min_years,
    n_sites      = length(eligible_sites),
    n_sitemonths = nrow(model_dat),
    OR           = round(or, 3),
    CI_low       = round(ci[2, 1], 3),
    CI_high      = round(ci[2, 2], 3),
    chi_sq       = round(lrt$Chisq[2], 2),
    p_value      = lrt$`Pr(>Chisq)`[2]
  )
}

# Run sensitivity analysis across three thresholds
# Note: this will take several minutes due to rowwise calculation
sensitivity_results <- bind_rows(
  run_clustering_model(min_years = 2),
  run_clustering_model(min_years = 3),
  run_clustering_model(min_years = 5)
)

print(sensitivity_results, width = Inf)

## make summary table

#### Supplementary Table - Sensitivity analysis results ####

sensitivity_table <- sensitivity_results %>%
  mutate(
    `Minimum coverage threshold (years)` = min_years,
    `Sites included (n)`                 = n_sites,
    `Site-months (n)`                    = n_sitemonths,
    `Odds ratio`                         = OR,
    `95% CI`                             = paste0(CI_low, "–", CI_high),
    `χ² (df = 1)`                        = chi_sq,
    `p-value`                            = format(p_value, scientific = TRUE, digits = 2)
  ) %>%
  dplyr::select(
    `Minimum coverage threshold (years)`,
    `Sites included (n)`,
    `Site-months (n)`,
    `Odds ratio`,
    `95% CI`,
    `χ² (df = 1)`,
    `p-value`
  )

print(sensitivity_table, width = Inf)

# Export as CSV for supplement
write.csv(
  sensitivity_table,
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/exported data summaries/supplementary_sensitivity_analysis.csv",
  row.names = FALSE
)

#### ----------------------------------------------------------------- ####
####     SPATIOTEMPORAL ANALYSIS: REGIONAL MASTING PERIOD IDENTIFICATION ####
#### ----------------------------------------------------------------- ####

#### Step 24 - Justify gap threshold using quiet period distribution ####
# Plots distribution of quiet period lengths to support the 6-month
# gap threshold used to split regional masting periods.

# First active_by_month pass: simple version for quiet run analysis only
active_by_month_simple <- all_harmonised %>%
  filter(!is.na(month_date)) %>%
  group_by(month_date) %>%
  summarise(
    n_sites_active   = n_distinct(site[phenophase %in% c("FL", "FR")]),
    n_sites_observed = n_distinct(site[!is.na(phenophase)]),
    .groups = "drop"
  ) %>%
  arrange(month_date) %>%
  complete(
    month_date = seq(min(month_date), max(month_date), by = "month"),
    fill = list(n_sites_active = 0, n_sites_observed = 0)
  )

quiet_runs <- active_by_month_simple %>%
  mutate(is_quiet = n_sites_active == 0) %>%
  mutate(
    run_id  = with(rle(is_quiet), rep(seq_along(lengths), lengths)),
    run_len = with(rle(is_quiet), rep(lengths, lengths))
  ) %>%
  filter(is_quiet) %>%
  distinct(run_id, run_len)

# Distribution of quiet run lengths
quiet_runs %>%
  ggplot(aes(x = run_len)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "grey40") +
  geom_vline(xintercept = 6, linetype = "dashed", colour = "red") +
  labs(
    x     = "Length of quiet period (months)",
    y     = "Count",
    title = "Distribution of quiet period lengths between active months",
    subtitle = "Red dashed line = 6-month threshold used to split regional masting periods"
  ) +
  theme_bw()

# Summary stats
quiet_runs %>%
  summarise(
    min    = min(run_len),
    median = median(run_len),
    mean   = mean(run_len),
    max    = max(run_len),
    n_runs = n()
  )

quiet_runs %>%
  dplyr::count(run_len) %>%
  arrange(run_len) %>%
  print(n = 40)


#### Step 25 - Build full active_by_month with FL/FR/none counts ####
# Full version needed for timeline plots and period identification.

active_by_month <- all_harmonised %>%
  filter(!is.na(month_date)) %>%
  group_by(month_date) %>%
  summarise(
    n_flowering      = n_distinct(site[phenophase == "FL"]),
    n_fruiting       = n_distinct(site[phenophase == "FR"]),
    n_none           = n_distinct(site[phenophase == "none"]),
    n_sites_active   = n_distinct(site[phenophase %in% c("FL", "FR")]),
    n_sites_observed = n_distinct(site[!is.na(phenophase)]),
    .groups = "drop"
  ) %>%
  arrange(month_date) %>%
  complete(
    month_date = seq(min(month_date), max(month_date), by = "month"),
    fill = list(
      n_flowering      = 0,
      n_fruiting       = 0,
      n_none           = 0,
      n_sites_active   = 0,
      n_sites_observed = 0
    )
  )


#### Step 26 - Identify regional masting periods ####
# Gap threshold: 6 consecutive months of no active sites splits periods.
# Minimum sites threshold: at least 3 sites active at some point.

mast_timeline <- active_by_month %>%
  mutate(
    is_active  = n_sites_active > 0,
    block_id   = with(rle(is_active), rep(seq_along(lengths), lengths)),
    block_len  = with(rle(is_active), rep(lengths, lengths)),
    block_type = if_else(is_active, "active", "quiet")
  )

splitting_gaps <- mast_timeline %>%
  filter(block_type == "quiet", block_len >= 6) %>%
  distinct(block_id)

mast_timeline <- mast_timeline %>%
  mutate(
    is_split   = block_id %in% splitting_gaps$block_id,
    new_period = block_type == "active" &
      (lag(is_split, default = TRUE) | row_number() == 1),
    period_id  = cumsum(new_period),
    period_id  = if_else(is_active, period_id, NA_integer_)
  )

period_summary <- mast_timeline %>%
  filter(!is.na(period_id)) %>%
  group_by(period_id) %>%
  summarise(
    period_start      = min(month_date),
    period_end        = max(month_date),
    n_months          = n(),
    max_sites_active  = max(n_sites_active),
    total_site_months = sum(n_sites_active),
    .groups           = "drop"
  ) %>%
  filter(max_sites_active >= 3)

period_summary %>%
  dplyr::select(period_id, period_start, period_end,
                n_months, max_sites_active) %>%
  print(n = 30)


#### Step 27 - Build plot data for timeline figures ####

# Long format FL/FR counts per month for stacked bars above x-axis
active_long <- active_by_month %>%
  dplyr::select(month_date, n_flowering, n_fruiting) %>%
  pivot_longer(
    cols      = c(n_flowering, n_fruiting),
    names_to  = "phenophase",
    values_to = "n_sites"
  ) %>%
  mutate(
    phenophase = recode(phenophase,
                        "n_flowering" = "Flowering",
                        "n_fruiting"  = "Fruiting")
  )

# None counts per month as negative values for bars below x-axis
none_data <- active_by_month %>%
  dplyr::select(month_date, n_none) %>%
  mutate(n_sites = -n_none)

# Shaded rectangles for identified regional masting periods
period_rects <- period_summary %>%
  mutate(
    xmin = period_start - days(15),
    xmax = period_end   + days(15)
  )


#### Step 29 - Regional masting timeline ####
# Two export versions: split (two rows) and single row for multi-panel figure.
# Masting period shading, observed non-masting bars, and FL/FR bars above axis.
# Consistent y-axis scale across rows. Legend includes all four data elements.

library(grid)

# Update active_long with full phenophase labels for legend
active_long <- active_by_month %>%
  dplyr::select(month_date, n_flowering, n_fruiting) %>%
  pivot_longer(
    cols      = c(n_flowering, n_fruiting),
    names_to  = "phenophase",
    values_to = "n_sites"
  ) %>%
  mutate(
    phenophase = recode(phenophase,
                        "n_flowering" = "Flowering (FL)",
                        "n_fruiting"  = "Fruiting (FR)")
  )

# Consistent y-axis scale across all rows
y_max <- max(
  max(active_long$n_sites, na.rm = TRUE),
  max(abs(none_data$n_sites), na.rm = TRUE)
)

# Fixed x-axis breaks every 5 years
x_breaks_fixed <- seq(
  as.Date("1970-01-01"),
  as.Date("2030-01-01"),
  by = "5 years"
)

# Helper function to build one row of the timeline
make_timeline_row <- function(active_data, none_data, rects,
                              date_from, date_to,
                              show_x_label = FALSE) {
  
  x_label <- if (show_x_label) "Year" else ""
  
  ggplot() +
    
    # masting period shading — in legend via fill aes
    geom_rect(
      data = rects %>% filter(xmin <= date_to, xmax >= date_from),
      aes(xmin = pmax(xmin, date_from),
          xmax = pmin(xmax, date_to),
          ymin = -Inf, ymax = Inf,
          fill = "Regional masting period"),
      colour = NA,
      alpha  = 1
    ) +
    
    # none bars below x-axis — in legend via fill aes
    geom_col(
      data  = none_data %>% filter(month_date >= date_from,
                                   month_date <  date_to),
      aes(x = month_date, y = n_sites,
          fill = "Observed non-masting"),
      colour = NA,
      width  = 31,
      alpha  = 0.9
    ) +
    
    # flowering and fruiting bars above x-axis
    geom_col(
      data     = active_data %>%
        filter(n_sites > 0,
               month_date >= date_from,
               month_date <  date_to),
      aes(x = month_date, y = n_sites, fill = phenophase),
      width    = 31,
      colour   = NA,
      alpha    = 0.85,
      position = "stack"
    ) +
    
    scale_fill_manual(
      name   = NULL,
      values = c(
        "Flowering (FL)"          = col_timeline_fl,
        "Fruiting (FR)"           = col_timeline_fr,
        "Observed non-masting"    = "grey55",
        "Regional masting period" = "#E8EEF4"
      ),
      breaks = c(
        "Flowering (FL)",
        "Fruiting (FR)",
        "Observed non-masting",
        "Regional masting period"
      )
    ) +
    
    geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
    
    scale_x_date(
      breaks      = x_breaks_fixed,
      date_labels = "%Y",
      limits      = c(date_from, date_to),
      expand      = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      labels = function(x) abs(x),
      limits = c(-y_max, y_max),
      breaks = c(-10, -5, 0, 5, 10)
    ) +
    labs(x = x_label, y = "") +
    theme_masting() +
    theme(
      panel.grid      = element_blank(),
      panel.border    = element_rect(colour = "grey60", fill = NA),
      legend.position = "none"
    )
}

#### Single row version for multi-panel figure ####

row_single <- make_timeline_row(
  active_data  = active_long, none_data = none_data, rects = period_rects,
  date_from    = as.Date("1972-01-01"), date_to = as.Date("2027-01-01"),
  show_x_label = TRUE
)

regional_timeline_single <- row_single +
  theme(
    legend.position  = "top",
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_legend(title = NULL, nrow = 1)) +
  labs(y = "Number of sites active")

regional_timeline_single

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/regional_masting_timeline_single.png",
  plot   = regional_timeline_single,
  width  = 6.8,
  height = 2.5,
  dpi    = fig_dpi
)


#### ----------------------------------------------------------------- ####
####     SPATIOTEMPORAL ANALYSIS: DATASET COVERAGE TIMELINE            ####
#### ----------------------------------------------------------------- ####

#### Step 30 - Site observation timeline colored by data source type ####
# Shows which datasets cover each site and when, with masting events
# overlaid as colored tiles.

timeline_df <- all_harmonised %>%
  filter(!is.na(month_date)) %>%
  mutate(
    source_type = dplyr::recode(dataset_name, !!!source_type_lookup)
  ) %>%
  distinct(source_type, dataset_name, site, month_date, phenophase) %>%
  mutate(
    site = factor(site, levels = all_harmonised %>%
                    filter(!is.na(month_date)) %>%
                    group_by(site) %>%
                    summarise(first_obs = min(month_date), .groups = "drop") %>%
                    arrange(first_obs) %>%
                    pull(site))
  )

activity_band <- timeline_df %>%
  filter(phenophase %in% c("FL", "FR")) %>%
  distinct(site, month_date, phenophase)

sea_sites_temporal <- ggplot() +
  geom_line(
    data      = timeline_df,
    aes(x = month_date, y = site, colour = source_type),
    linewidth = 2.5,
    lineend   = "square",
    alpha     = 0.8
  ) +
  geom_tile(
    data   = activity_band,
    aes(x = month_date, y = site, fill = phenophase),
    height = 0.85,
    width  = 40,
    alpha  = 0.9
  ) +
  scale_colour_manual(
    name   = "Data source",
    values = source_cols
  ) +
  scale_fill_manual(
    name   = "Masting event",
    values = c("FL" = col_timeline_fl, "FR" = col_timeline_fr),
    labels = c("FL" = "Flowering", "FR" = "Fruiting")
  ) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y",
    expand      = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    x        = "Year",
    y        = NULL,
    title    = "Temporal coverage of phenology datasets across Southeast Asia",
    subtitle = "Lines show observation coverage by data source type; bands show masting events"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9),
    legend.position    = "right",
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, colour = "grey40")
  ) +
  guides(
    colour = guide_legend(override.aes = list(linewidth = 3, alpha = 1)),
    fill   = guide_legend(override.aes = list(size = 4, alpha = 0.8))
  )

sea_sites_temporal

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/all_datasets_timeline_figure.png",
  plot   = sea_sites_temporal,
  width  = 10,
  height = 6,
  dpi    = 300
)


#### ----------------------------------------------------------------- ####
####     SPATIOTEMPORAL ANALYSIS: ONSET TIMING MAP BY MASTING PERIOD   ####
#### ----------------------------------------------------------------- ####

#### Step 31 - Identify candidate masting periods for onset map ####
# Additional filter beyond period_summary: ≥6 months AND ≥4 sites total.
# This ensures enough spatial coverage for meaningful onset patterns.

mast_period_data <- all_harmonised %>%
  filter(
    !is.na(month_date),
    phenophase %in% c("FL", "FR")
  ) %>%
  left_join(
    mast_timeline %>%
      filter(!is.na(period_id)) %>%
      dplyr::select(month_date, period_id),
    by = "month_date"
  ) %>%
  filter(!is.na(period_id))

period_eligibility <- mast_period_data %>%
  group_by(period_id) %>%
  summarise(
    n_months      = n_distinct(month_date),
    n_sites_total = n_distinct(site),
    .groups       = "drop"
  ) %>%
  filter(n_months >= 6, n_sites_total >= 4)

keep_periods   <- period_eligibility$period_id
period_ids_chr <- as.character(sort(keep_periods))

cat("Eligible masting periods:", length(keep_periods), "\n")
print(period_eligibility)


#### Step 32 - Calculate site onset per period ####
# (site_coords defined earlier, before its first use in the lat/deviation section)

site_onset <- mast_period_data %>%
  filter(period_id %in% keep_periods) %>%
  group_by(period_id, site) %>%
  summarise(
    onset_date = min(month_date),
    .groups    = "drop"
  ) %>%
  left_join(site_coords, by = "site") %>%
  filter(!is.na(lat), !is.na(long)) %>%
  group_by(period_id) %>%
  arrange(onset_date, .by_group = TRUE) %>%
  mutate(
    onset_rank   = dense_rank(onset_date),
    onset_scaled = scales::rescale(
      onset_rank,
      to   = c(0, 1),
      from = range(onset_rank, na.rm = TRUE)
    )
  ) %>%
  ungroup() %>%
  mutate(period_id = as.character(period_id))




#### Step 33 - Build facet labels with onset date range ####


facet_labs <- site_onset %>%
  group_by(period_id) %>%
  summarise(
    onset_min = format(min(onset_date), "%b %Y"),
    onset_max = format(max(onset_date), "%b %Y"),
    .groups   = "drop"
  ) %>%
  mutate(
    facet_lab = paste0(onset_min, " \u2192 ", onset_max)
  ) %>%
  { setNames(.$facet_lab, .$period_id) }

#### Step 34 - Assign colors per period ####

event_cols <- c(
  "1"  = "#FF4E00",
  "4"  = "#F5D547",
  "7"  = "#B1693B",
  "8"  = "#878844",
  "10" = "#33673B",
  "11" = "#638475",
  "13" = "#026B79",
  "16" = "#7B506F",
  "18" = "#1F1A38",
  "22" = "#6AD5CB"   # new period — added teal
)

site_onset <- site_onset %>%
  mutate(period_id = factor(period_id, levels = period_ids_chr))

#### Step 35 - Basemap ####

sea_coast <- rnaturalearth::ne_coastline(
  scale       = "medium",
  returnclass = "sf"
)


#### Step 36 - Faceted onset timing map ####
# Each panel = one regional masting period.
# Color = period identity; opacity = onset order (darker = earlier).

onset_map <- ggplot() +
  geom_sf(
    data      = sea_coast,
    color     = "grey60",
    linewidth = 0.3
  ) +
  geom_point(
    data = site_onset,
    aes(x     = long,
        y     = lat,
        color = period_id,
        alpha = onset_scaled),
    size = 3.5
  ) +
  facet_wrap(
    ~ period_id,
    labeller = labeller(period_id = facet_labs)
  ) +
  coord_sf(
    xlim   = c(95, 120),
    ylim   = c(-5, 10),
    expand = FALSE
  ) +
  scale_color_manual(
    values = event_cols,
    guide  = "none"  # remove mast event legend
  ) +
  scale_alpha_continuous(
    name     = "Onset timing (early → late)",
    range    = c(0.95, 0.35),
    breaks   = c(0, 0.25, 0.5, 0.75, 1.0),
    na.value = 0.25
  ) +
  scale_x_continuous(breaks = c(100, 110, 120)) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10)) +
  # update facet labels to show only date range, not event ID
  theme_masting() +
  theme(
    panel.grid         = element_blank(),
    panel.border       = element_rect(colour = "grey80", fill = NA,
                                      linewidth = 0.3),
    strip.background   = element_blank(),
    strip.text         = element_text(size = 8, colour = "grey30"),
    axis.text          = element_text(size = 7),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 9),
    legend.text        = element_text(size = 8),
    plot.margin        = margin(5, 15, 5, 5)  # extra right margin
  ) +
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  guides(
    alpha = guide_legend(
      title          = "Onset timing (early → late)",
      title.position = "top",
      direction      = "horizontal",
      nrow           = 1
    )
  )

onset_map

ggsave(
  "/Users/zoe/Library/CloudStorage/Dropbox/Zoë PhD projects/1_Chapters_ZL/Ch3 ZL - pheno map/phenology_r/output/masting_onset_map.png",
  plot   = onset_map,
  width  = fig_width,
  height = fig_height_triple,
  dpi    = fig_dpi
)



#### ----------------------------------------------------------------- ####
#### Site summary table (with thresholds) is now built separately in    ####
#### site_summary_table_long.R - see exported data summaries/           ####
#### site_summary_table_long.csv for the current version.               ####
#### ----------------------------------------------------------------- ####