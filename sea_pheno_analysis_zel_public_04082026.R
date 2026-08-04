####   Southeast Asian masting dataset - Analysis Script                ####
####   Zoë Evans Lieb                                                   ####
####   Last updated: 4 August 2026                                      ####
# ----------------------------------------------------------------- #
#
# This script runs the full analysis reported in the manuscript, operating
# on the published site-month dataset. Per-dataset processing scripts
# (thresholding of raw tree-level data for each contributing site) are not
# included here as these data are private; the published dataset is the
# output of those scripts.
#
# STRUCTURE
#   Setup                Libraries, paths, constants, theme, basemap
#   Steps 1-2            Read data
#   Steps 3-6            Lookup tables and figure constants
#   Step  7              Dataset coverage timeline (descriptive)
#   Steps 8-10           Masting activity by calendar month
#   Steps 11-19          Circular analysis of onset timing
#   Steps 20-24          Masting frequency and duration by site
#   Steps 25-29          Temporal clustering of masting across sites
#   Steps 30-33          Regional masting period identification
#   Steps 34-36          Onset timing maps
#
# TO RUN
#   1. Download this script and the two data files from the repository:
#        sea_phenology_combined_data.csv
#        sea_masting_literature_review_events.csv
#   2. Put all three in the same folder (if downloading from GitHub,
#       they will be in a nested folder that needs to be removed).
#   3. Open this script in RStudio, set that folder as your working
#      directory (Session > Set Working Directory > To Source File
#      Location), and run the script top to bottom.
#   4. Figures and tables are written to an "output" folder, created
#      automatically inside that same directory.
#
# No paths need to be edited for the default setup above. If your files
# are somewhere else, change data_dir below.

# ----------------------------------------------------------------- #
####                        SETUP                               ####
# ----------------------------------------------------------------- #

#### Required packages ####
# Install any of these you don't already have with:
#   install.packages(c("tidyverse", "lubridate", "ggrepel", "geosphere",
#                       "rnaturalearth", "rnaturalearthdata", "sf",
#                       "circular", "diptest", "patchwork", "scales",
#                       "colorspace", "lme4", "broom.mixed", "scico"))

library(tidyverse)
library(lubridate)
library(ggrepel)
library(geosphere)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(circular)
library(diptest)
library(patchwork)
library(scales)
library(colorspace)
library(lme4)
library(broom.mixed)
library(scico)
library(grid)


#### Paths ####
# By default, this script looks for the two data files in the same folder
# as itself (your working directory). Figures and tables are written to
# an "output" folder here, created automatically.
#
# If your files live somewhere else, change data_dir below to that path.

data_dir   <- getwd()
output_dir <- file.path(getwd(), "output")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

required_files <- c("sea_phenology_combined_data.csv",
                    "sea_masting_literature_review_events.csv")
missing_files  <- required_files[!file.exists(file.path(data_dir, required_files))]

if (length(missing_files) > 0) {
  stop(
    "Could not find the following file(s) in ", data_dir, ":\n  ",
    paste(missing_files, collapse = "\n  "),
    "\n\nDownload them from the repository and place them in the same ",
    "folder as this script, or set data_dir above to point at them."
  )
}


# ----------------------------------------------------------------- #
####                         READ DATA                          ####
# ----------------------------------------------------------------- #

#### Step 1 - Read harmonised site-month dataset ####
# One row per site per month per dataset, with phenophase classified as
# FL (flowering), FR (fruiting), or none (observed, not masting).

all_harmonised <- read_csv(
  file.path(data_dir, "sea_phenology_combined_data.csv"),
  show_col_types = FALSE
)


#### Step 2 - Read literature review event records ####
# One row per discrete literature-reported masting event, with explicit
# onset_date and end_date. Used for onset analysis only, where explicit dates are
# more reliable than detecting onsets from consecutive observed months.

pheno <- read_csv(
  file.path(data_dir, "sea_masting_literature_review_events.csv"),
  show_col_types = FALSE
)


# ----------------------------------------------------------------- #
####                 LOOKUP TABLES AND CONSTANTS                ####
# ----------------------------------------------------------------- #

#### Step 3 - Site name recode and exclusions ####
# pheno site names are recoded to match all_harmonised. Sites within 5 km of
# each other are lumped under a single site name.

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

# Sites excluded from onset analysis: no clear onset month available, and
# three of the four are in seasonal mainland Southeast Asia, outside the
# everwet dipterocarp scope of this study.
excluded_sites <- c(
  "Doi Suthep-Pui National Park",
  "Khao Ang Runai Wildlife Sanctuary",
  "Khun Wang Royal Agriculture Research Center",
  "Lobo"
)

# Literature records excluded from onset analysis: only a peak month was
# reported, not a true onset month.
exclude_cov_ids <- c(85, 359, 405, 108, 474, 74, 17)


#### Step 4 - Data source type lookup ####

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

# Any dataset_name not in the lookup would pass through dplyr::recode()
# unchanged and appear as its own category. Fail loudly instead.
stopifnot(
  length(setdiff(unique(all_harmonised$dataset_name),
                 names(source_type_lookup))) == 0
)

# Datasets contributing tree-level or expert-reported site-month states.
# Used to select onset flags and to prioritise datasets where a site appears
# in more than one source.
monitoring_datasets <- c("FRIM", "Lambir", "Gunung_Palung", "Pasoh",
                         "Bukit_Timah", "additional_contributed_obs",
                         "BBBR", "Gunung_Tarak")


#### Step 5 - Figure colours ####
# All figure colours defined here once, for consistency across the paper.

# Phenophase colours - analysis figures
col_FL        <- "#6C739D"   # Glaucous - flowering
col_FR        <- "#7FBEAB"   # Muted teal - fruiting
col_FL_sig    <- "#6C739D"
col_FL_nonsig <- "#B8BCE0"
col_FR_sig    <- "#7FBEAB"
col_FR_nonsig <- "#B8DDD4"
col_FL_border <- "#3D4270"
col_FR_border <- "#4A8C7A"

# Phenophase colours - timeline and rose figures (higher contrast)
col_timeline_fl <- "#916392"   # medium purple - flowering
col_timeline_fr <- "#371B38"   # dark purple - fruiting

# Data source colours
col_ltm       <- "#8ECCA5"   # celadon
col_unpub_ltm <- "#F0BE19"   # saffron
col_expert    <- "#D37D56"   # burnt peach
col_lit       <- "#708AFF"   # cornflower blue

# Named vectors
phase_cols          <- c("FL" = col_FL, "FR" = col_FR)
phase_cols_timeline <- c("FL" = col_timeline_fl, "FR" = col_timeline_fr)
phase_labels        <- c("FL" = "Flowering", "FR" = "Fruiting")

source_cols <- c(
  "Long-term monitoring"  = col_ltm,
  "Unpublished long-term" = col_unpub_ltm,
  "Expert observations"   = col_expert,
  "Literature review"     = col_lit
)


#### Step 6 - Theme and export dimensions ####
# Target: Ecology Letters full page width (173 mm = 6.8 inches).
# base_size = 11 calibrated for this width.

theme_masting <- function(base_size = 11) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      axis.text        = element_text(size = base_size * 0.9),
      axis.title       = element_text(size = base_size),
      legend.text      = element_text(size = base_size * 0.9),
      legend.title     = element_text(size = base_size, face = "bold"),
      strip.text       = element_text(size = base_size * 0.9, face = "bold"),
      plot.title       = element_text(size = base_size * 1.1, face = "bold"),
      plot.subtitle    = element_text(size = base_size * 0.85,
                                      colour = "grey40"),
      plot.margin      = margin(6, 6, 6, 6)
    )
}

fig_width         <- 6.8   # inches (173 mm)
fig_height_single <- 4     # single panel
fig_height_double <- 6     # two panels
fig_height_triple <- 8     # three panels or tall figures
fig_dpi           <- 300


#### Site coordinates and basemap ####
# Defined here because both are used by several sections below.

site_coords <- all_harmonised %>%
  filter(!is.na(lat), !is.na(long)) %>%
  group_by(site) %>%
  summarise(
    long    = median(long, na.rm = TRUE),
    lat     = median(lat,  na.rm = TRUE),
    .groups = "drop"
  )

sea_coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")


# ----------------------------------------------------------------- #
####                 DESCRIPTIVE: DATASET COVERAGE              ####
# ----------------------------------------------------------------- #

#### Step 7 - Site observation timeline by data source ####
# Lines show observation coverage per site, coloured by source type.
# Overlaid tiles show months classified as masting.

timeline_df <- all_harmonised %>%
  filter(!is.na(month_date)) %>%
  mutate(source_type = dplyr::recode(dataset_name, !!!source_type_lookup)) %>%
  distinct(source_type, dataset_name, site, month_date, phenophase) %>%
  mutate(
    site = factor(
      site,
      levels = all_harmonised %>%
        filter(!is.na(month_date)) %>%
        group_by(site) %>%
        summarise(first_obs = min(month_date), .groups = "drop") %>%
        arrange(first_obs) %>%
        pull(site)
    )
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
  scale_colour_manual(name = "Data source", values = source_cols) +
  scale_fill_manual(
    name   = "Masting event",
    values = phase_cols_timeline,
    labels = phase_labels
  ) +
  scale_x_date(
    date_breaks = "5 years",
    date_labels = "%Y",
    expand      = expansion(mult = c(0.01, 0.02))
  ) +
  labs(x = "Year", y = NULL) +
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
  file.path(output_dir, "all_datasets_timeline_figure.png"),
  plot   = sea_sites_temporal,
  width  = fig_width,
  height = fig_height_double,
  dpi    = fig_dpi
)

# ----------------------------------------------------------------- #

####         MASTING ACTIVITY BY CALENDAR MONTH                 ####

# ----------------------------------------------------------------- #

#### Step 8 - Proportion of observed site-months masting, by month ####
# Divides masting site-months by total observed site-months per calendar
# month, controlling for uneven sampling effort across the year.
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


#### Step 9 - Binomial tests against the overall baseline ####
# Tests whether each calendar month's proportion of masting site-months is
# significantly above the overall observed proportion.

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
    baseline %>% rename(n_total = n_observed, n_total_masting = n_masting),
    by = "phenophase"
  ) %>%
  rowwise() %>%
  mutate(
    binom_p = binom.test(
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
      phenophase == "FL" &  significant ~ "FL_sig",
      phenophase == "FL" & !significant ~ "FL_nonsig",
      phenophase == "FR" &  significant ~ "FR_sig",
      phenophase == "FR" & !significant ~ "FR_nonsig"
    )
  )

prop_counts_sig %>%
  dplyr::select(calendar_month, phenophase, prop_masting,
                baseline_prop, binom_p, significant) %>%
  arrange(phenophase, calendar_month) %>%
  print(n = 30)


#### Step 10 - Figure: proportion masting by calendar month ####
# Bars shaded by whether the month differs significantly from baseline.

prop_fill_vals <- c(
  "FL_sig"    = col_FL_sig,
  "FL_nonsig" = col_FL_nonsig,
  "FR_sig"    = col_FR_sig,
  "FR_nonsig" = col_FR_nonsig
)

prop_border_vals <- c(
  "FL_sig"    = col_FL_border,
  "FL_nonsig" = col_FL_nonsig,
  "FR_sig"    = col_FR_border,
  "FR_nonsig" = col_FR_nonsig
)

prop_fill_labs <- c(
  "FL_sig"    = "Flowering (p < 0.05)",
  "FL_nonsig" = "Flowering (p \u2265 0.05)",
  "FR_sig"    = "Fruiting (p < 0.05)",
  "FR_nonsig" = "Fruiting (p \u2265 0.05)"
)

p_prop <- ggplot(prop_counts_sig,
                 aes(x = calendar_month, y = prop_masting)) +
  geom_col(
    aes(fill = fill_group, colour = fill_group),
    position = "dodge",
    width    = 0.8
  ) +
  scale_fill_manual(name = NULL, values = prop_fill_vals,
                    labels = prop_fill_labs) +
  scale_colour_manual(name = NULL, values = prop_border_vals,
                      labels = prop_fill_labs) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Month", y = "Proportion of observed site-months masting") +
  theme_masting()

p_prop

ggsave(
  file.path(output_dir, "proportion_masting_by_calendar_month.png"),
  plot   = p_prop,
  width  = fig_width,
  height = fig_height_single,
  dpi    = fig_dpi
)


# ----------------------------------------------------------------- #
####          CIRCULAR ANALYSIS OF MASTING ONSET TIMING           ####
# ----------------------------------------------------------------- #
#
# Analysis order follows the order in which the finding emerged:
#   Step 11  Build onset dataset
#   Step 12  Pooled Rayleigh tests across all sites
#   Step 13  Test onset distributions for multimodality
#   Step 14  Split by latitude, repeat circular tests
#   Step 15  Test whether the two groups differ from each other
#   Step 16  Robustness checks
#   Steps 17-19  Figures


#### Step 11 - Build all_onsets_full ####
# One row per masting onset event (FL or FR) per site.
# Monitoring and expert obs: onset detected from is_fl_onset / is_fr_onset.
# Literature review: onset taken directly from pheno's explicit onset_date.

fl_onsets_ltm <- all_harmonised %>%
  filter(dataset_name %in% monitoring_datasets, is_fl_onset == TRUE) %>%
  mutate(phase = "FL") %>%
  dplyr::select(site, onset_date = month_date, phase, dataset_name, area)

fr_onsets_ltm <- all_harmonised %>%
  filter(dataset_name %in% monitoring_datasets, is_fr_onset == TRUE) %>%
  mutate(phase = "FR") %>%
  dplyr::select(site, onset_date = month_date, phase, dataset_name, area)

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

all_onsets_full <- bind_rows(
  fl_onsets_ltm, fr_onsets_ltm, fl_onsets_lit, fr_onsets_lit
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

# Every site with onsets must have coordinates, or it will be dropped
# silently by the joins below. Should return character(0).
stopifnot(
  length(setdiff(unique(all_onsets_full$site), unique(site_coords$site))) == 0
)


#### Step 12 - Pooled Rayleigh tests across all sites ####
# Tests whether onset months are non-uniformly distributed around the
# calendar year, pooling all sites.
#
# NOTE: reported in the manuscript as the starting point, not as a finding.
# It understates seasonality because it pools sites that reproduce roughly
# five months apart (see Step 14).

fl_circular <- circular(
  all_onsets_full %>% filter(phase == "FL") %>%
    pull(onset_month) * (2 * pi / 12),
  units = "radians", template = "none"
)
rayleigh_fl   <- rayleigh.test(fl_circular)
mean_month_fl <- (mean(fl_circular) * 12 / (2 * pi)) %% 12

fr_circular <- circular(
  all_onsets_full %>% filter(phase == "FR") %>%
    pull(onset_month) * (2 * pi / 12),
  units = "radians", template = "none"
)
rayleigh_fr   <- rayleigh.test(fr_circular)
mean_month_fr <- (mean(fr_circular) * 12 / (2 * pi)) %% 12

rayleigh_pooled <- tibble::tibble(
  phase      = c("FL", "FR"),
  n          = c(sum(all_onsets_full$phase == "FL"),
                 sum(all_onsets_full$phase == "FR")),
  r          = c(round(rayleigh_fl$statistic, 3),
                 round(rayleigh_fr$statistic, 3)),
  p_value    = c(rayleigh_fl$p.value, rayleigh_fr$p.value),
  mean_month = c(round(mean_month_fl, 1), round(mean_month_fr, 1))
)

print(rayleigh_pooled)


#### Step 13 - Test onset distributions for multimodality ####
# A circular mean is only interpretable if onsets are distributed around a
# single centre. Hartigans' dip test compares the observed distribution to
# the best-fitting unimodal distribution.
#
# Both phenophases depart significantly from unimodality, which is why the
# pooled means above are not treated as describing the regional distribution.

fl_onsets <- filter(all_onsets_full, phase == "FL")
fr_onsets <- filter(all_onsets_full, phase == "FR")

# FL: Mar-Apr peak, Jul trough, Aug-Oct second cluster
table(fl_onsets$onset_month)
table(fr_onsets$onset_month)

dip_fl <- diptest::dip.test(fl_onsets$onset_month)
dip_fr <- diptest::dip.test(fr_onsets$onset_month)

print(dip_fl)
print(dip_fr)


#### Step 14 - Split sites by latitude and repeat circular tests ####
# The two flowering modes correspond to site latitude. Sites are assigned to
# a northern group (lat >= 1 N) or a southern/equatorial group (lat < 1 N),
# and Rayleigh's test is repeated within each group and phenophase.

HEMI_CUTOFF <- 1   # degrees latitude; sensitivity-tested in Step 16

onsets_hemi <- all_onsets_full %>%
  left_join(site_coords, by = "site") %>%
  filter(!is.na(lat), !is.na(onset_month)) %>%
  mutate(
    hemisphere = if_else(lat >= HEMI_CUTOFF,
                         "Northern", "Southern/Equatorial")
  )

# Confirm the coordinate join dropped nothing
stopifnot(nrow(onsets_hemi) == nrow(all_onsets_full))
onsets_hemi %>% dplyr::count(hemisphere, area, phase) %>% print(n = 40)


# Helper: Rayleigh test and circular mean for any grouped data frame.
# Custom function, not from a package. Must be run before use.
rayleigh_by_group <- function(df) {
  df %>%
    group_modify(~ {
      theta <- circular::circular(.x$onset_month * 2 * pi / 12,
                                  type = "angles", units = "radians")
      rt <- circular::rayleigh.test(theta)
      mm <- as.numeric(circular::mean.circular(theta)) * 12 / (2 * pi)
      tibble::tibble(
        n          = nrow(.x),
        r          = rt$statistic,
        p          = rt$p.value,
        mean_month = ifelse(mm < 0, mm + 12, mm)   # unwrap negative angles
      )
    })
}

# MAIN RESULT: seasonality within each latitudinal group
rayleigh_hemisphere <- onsets_hemi %>%
  group_by(hemisphere, phase) %>%
  rayleigh_by_group()

print(rayleigh_hemisphere)

# Supporting descriptive: which onset window each site uses, by latitude
onset_window_by_site <- onsets_hemi %>%
  filter(phase == "FL") %>%
  mutate(window = case_when(
    onset_month %in% 2:6  ~ "Feb-Jun",
    onset_month %in% 8:11 ~ "Aug-Nov",
    TRUE                  ~ "other"
  )) %>%
  dplyr::count(site, lat, window) %>%
  pivot_wider(names_from = window, values_from = n, values_fill = 0) %>%
  arrange(lat)

print(onset_window_by_site, n = 30)


#### Step 15 - Test whether the two groups differ from each other ####
# Rayleigh (Step 14) establishes that each group is seasonal. It does NOT
# establish that the groups are seasonal at DIFFERENT times of year.
# Watson's two-sample test of homogeneity compares the two distributions
# directly. Non-parametric; makes no distributional assumption.
# Critical value at alpha = 0.001 is approximately 0.268.

to_circ <- function(x) {
  circular::circular(x * 2 * pi / 12, type = "angles", units = "radians")
}

fl_hemi <- filter(onsets_hemi, phase == "FL")
fr_hemi <- filter(onsets_hemi, phase == "FR")

watson_fl <- circular::watson.two.test(
  to_circ(fl_hemi$onset_month[fl_hemi$hemisphere == "Northern"]),
  to_circ(fl_hemi$onset_month[fl_hemi$hemisphere == "Southern/Equatorial"])
)

watson_fr <- circular::watson.two.test(
  to_circ(fr_hemi$onset_month[fr_hemi$hemisphere == "Northern"]),
  to_circ(fr_hemi$onset_month[fr_hemi$hemisphere == "Southern/Equatorial"])
)

print(watson_fl)
print(watson_fr)


#### Step 16 - Robustness checks ####
# Reported in Supplementary Materials.
#   16a. Does the result depend on where the dividing latitude is drawn?
#   16b. Does it survive removing the largest southern contributor?
#   16c. Is it detectable without any latitudinal grouping at all?

# 16a - Cutoff sensitivity.
# Only one site (Gunung Tarak, 0.4 N) falls between 0 and 1 degree, so the
# 0.5 and 1.0 cutoffs produce identical groupings.
cutoff_sensitivity <- purrr::map_dfr(c(0, 0.5, 1), function(cut) {
  onsets_hemi %>%
    mutate(hemisphere = if_else(lat >= cut,
                                "Northern", "Southern/Equatorial")) %>%
    group_by(hemisphere, phase) %>%
    rayleigh_by_group() %>%
    mutate(cutoff = cut, .before = 1)
})

print(cutoff_sensitivity, n = 20)

# 16b - Robustness to the largest southern contributor
drop_gp_sensitivity <- onsets_hemi %>%
  filter(site != "Gunung Palung National Park") %>%
  group_by(hemisphere, phase) %>%
  rayleigh_by_group()

print(drop_gp_sensitivity)

# 16c - By administrative area, with no latitudinal grouping applied.
# Restricted to areas with >= 8 onset events; a Rayleigh test on 2 events
# is not interpretable.
rayleigh_by_area <- onsets_hemi %>%
  group_by(area, phase) %>%
  filter(n() >= 8) %>%
  rayleigh_by_group() %>%
  arrange(phase, area)

print(rayleigh_by_area, n = 20)

write.csv(
  rayleigh_by_area,
  file.path(output_dir, "supplementary_rayleigh_by_area.csv"),
  row.names = FALSE
)


#### Step 17 - Figure: rose diagram of onset months by latitudinal group ####
# Main seasonality figure. Four panels: latitudinal group x phenophase.
# Radial scales are free, so wedge lengths are NOT comparable between panels.

rose_data <- onsets_hemi %>%
  mutate(
    phase_label = if_else(phase == "FL", "Flowering", "Fruiting"),
    hemi_label  = if_else(hemisphere == "Northern",
                          "North of 1\u00b0N", "South of 1\u00b0N")
  )

rose_plot <- ggplot(rose_data, aes(x = factor(onset_month, levels = 1:12))) +
  geom_bar(aes(fill = phase_label), width = 1,
           colour = "white", linewidth = 0.2) +
  coord_polar(start = 0) +
  facet_grid(hemi_label ~ phase_label, scales = "free_y") +
  scale_x_discrete(
    drop   = FALSE,
    labels = c("J","F","M","A","M","J","J","A","S","O","N","D")
  ) +
  scale_fill_manual(
    values = c(Flowering = col_timeline_fl, Fruiting = col_timeline_fr),
    guide  = "none"
  ) +
  labs(x = NULL, y = "Number of onset events") +
  theme_masting() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 10, colour = "grey20"),
    panel.spacing    = unit(1.2, "lines"),
    plot.margin      = margin(10, 10, 10, 10)
  )

rose_plot

ggsave(
  file.path(output_dir, "onset_rose_hemisphere.png"),
  plot   = rose_plot,
  width  = fig_width,
  height = fig_height_double,
  dpi    = fig_dpi
)

ggsave(
  file.path(output_dir, "onset_rose_hemisphere.svg"),
  plot   = rose_plot,
  width  = fig_width,
  height = fig_height_double,
  device = "svg"
)


#### Step 18 - Figure: circular plot, all sites pooled ####
# Shows the bimodal pooled distribution that motivated the latitudinal split.
# Supplementary figure.

circ_counts <- all_onsets_full %>%
  dplyr::count(onset_month, phase) %>%
  mutate(
    month_label = factor(month.abb[onset_month], levels = month.abb),
    phase       = factor(phase, levels = c("FL", "FR"))
  )

circ_all_sites <- ggplot(circ_counts,
                         aes(x = month_label, y = n, fill = phase)) +
  geom_col(position = "dodge", width = 0.8,
           alpha = 0.85, colour = "grey30") +
  coord_polar(start = -pi/12) +
  scale_fill_manual(name = "Phenophase",
                    values = phase_cols, labels = phase_labels) +
  scale_y_continuous(breaks = seq(0, 40, 10)) +
  labs(
    x        = NULL,
    y        = "Number of onset events",
    subtitle = paste0(
      "FL: r = ", round(rayleigh_fl$statistic, 2),
      ", p = ", format(rayleigh_fl$p.value, scientific = TRUE, digits = 2),
      "\nFR: r = ", round(rayleigh_fr$statistic, 2),
      ", p = ", format(rayleigh_fr$p.value, scientific = TRUE, digits = 2)
    )
  ) +
  theme_masting() +
  theme(
    axis.text.y     = element_blank(),
    panel.grid      = element_line(colour = "grey90"),
    legend.position = "bottom"
  )

circ_all_sites

ggsave(
  file.path(output_dir, "onset_circular_all_sites.png"),
  plot   = circ_all_sites,
  width  = fig_width,
  height = fig_height_single,
  dpi    = fig_dpi
)


#### Step 19 - Figure: onset events mapped ####
# One point per onset event, jittered around the site coordinate so that
# repeated events at the same site are visible. Jitter is +/- 0.7 degrees;
# state this in the figure caption. Supplementary figure.

set.seed(42)   # reproducible jitter

onset_events_map <- ggplot() +
  geom_sf(data = sea_coast, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  geom_jitter(
    data = onsets_hemi %>%
      mutate(phase_label = if_else(phase == "FL", "Flowering", "Fruiting")),
    aes(x = long, y = lat, colour = onset_month),
    width  = 0.7,
    height = 0.7,
    size   = 1.8,
    alpha  = 0.75
  ) +
  facet_wrap(~ phase_label) +
  coord_sf(xlim = c(95, 120), ylim = c(-5, 10), expand = FALSE) +
  scico::scale_colour_scico(
    palette = "romaO",            # cyclic palette; months wrap at Dec/Jan
    limits  = c(0.5, 12.5),
    breaks  = c(1, 4, 7, 10),
    labels  = c("Jan", "Apr", "Jul", "Oct"),
    name    = "Onset month"
  ) +
  scale_x_continuous(breaks = c(100, 110, 120)) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_masting() +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "grey80", fill = NA,
                                    linewidth = 0.3),
    strip.background = element_blank(),
    strip.text       = element_text(size = 9, colour = "grey30"),
    axis.text        = element_text(size = 7),
    legend.position  = "bottom",
    legend.direction = "horizontal"
  )

onset_events_map

ggsave(
  file.path(output_dir, "onset_events_map.png"),
  plot   = onset_events_map,
  width  = fig_width,
  height = fig_height_single,
  dpi    = fig_dpi
)


# ----------------------------------------------------------------- #
####          MASTING FREQUENCY AND DURATION BY SITE              ####
# ----------------------------------------------------------------- #

#### Step 20 - Define sites for frequency analysis ####
# Restricted to sites with effectively continuous monitoring coverage of at
# least 10 years, since masting can occur as infrequently as once per decade.
# FR analysis excludes Pasoh (FR data not available at time of analysis).

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

# Observation period per site - denominator for frequency
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


#### Step 21 - Flowering frequency and duration ####
# Discrete events identified by cumsum over onset flags for monitoring data,
# and from explicit onset/end dates for literature review sites.

fl_events_ltm <- all_harmonised %>%
  filter(
    site %in% fl_sites,
    dataset_name %in% monitoring_datasets,
    event_type == "FL",
    phenophase == "FL"
  ) %>%
  arrange(site, month_date) %>%
  group_by(site) %>%
  mutate(run_id = cumsum(is_fl_onset == TRUE)) %>%
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


#### Step 22 - Fruiting frequency and duration ####

fr_events_ltm <- all_harmonised %>%
  filter(
    site %in% fr_sites,
    dataset_name %in% monitoring_datasets,
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


#### Step 23 - Combined frequency and duration summary ####

freq_duration_summary <- fl_summary %>%
  dplyr::select(site, obs_years, n_fl_events, fl_frequency,
                mean_fl_duration, sd_fl_duration) %>%
  full_join(
    fr_summary %>%
      dplyr::select(site, n_fr_events, fr_frequency,
                    mean_fr_duration, sd_fr_duration),
    by = "site"
  ) %>%
  arrange(desc(obs_years))

print(freq_duration_summary)

# Ranges quoted in the results text
freq_duration_summary %>%
  summarise(
    fl_freq_min = min(fl_frequency,     na.rm = TRUE),
    fl_freq_max = max(fl_frequency,     na.rm = TRUE),
    fr_freq_min = min(fr_frequency,     na.rm = TRUE),
    fr_freq_max = max(fr_frequency,     na.rm = TRUE),
    fl_dur_min  = min(mean_fl_duration, na.rm = TRUE),
    fl_dur_max  = max(mean_fl_duration, na.rm = TRUE),
    fr_dur_min  = min(mean_fr_duration, na.rm = TRUE),
    fr_dur_max  = max(mean_fr_duration, na.rm = TRUE)
  )


#### Step 24 - Figure: frequency vs duration bubble chart ####
# Bubble size = years monitored; fill = subregion; shape = phenophase.

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
                  frequency = fl_frequency, duration = mean_fl_duration),
  fr_summary %>%
    mutate(phase = "FR") %>%
    dplyr::select(site, obs_years, phase,
                  frequency = fr_frequency, duration = mean_fr_duration)
) %>%
  mutate(
    subregion  = dplyr::recode(site, !!!site_subregion),
    site_label = dplyr::recode(site, !!!site_labels, .default = site)
  )

bubble_plot <- ggplot(
  bubble_data %>% filter(!is.na(frequency), !is.na(duration)),
  aes(x = frequency, y = duration, size = obs_years,
      fill = subregion, shape = phase)
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
    values = c("Peninsular Malaysia & Singapore" = "#7FBEAB",
               "Borneo"                          = "#2C365E"),
    labels = c("Peninsular Malaysia & Singapore" = "Peninsular Malaysia\n& Singapore",
               "Borneo"                          = "Borneo"),
    guide  = guide_legend(override.aes = list(shape = 21, size = 5))
  ) +
  scale_shape_manual(
    name   = "Phenophase",
    values = c("FL" = 21, "FR" = 23),
    labels = phase_labels,
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
  file.path(output_dir, "masting_frequency_duration_bubble.png"),
  plot   = bubble_plot,
  width  = fig_width,
  height = fig_height_single + 1,
  dpi    = fig_dpi
)


# ----------------------------------------------------------------- #
####     TEMPORAL CLUSTERING OF MASTING EVENTS ACROSS SITES        ####
# ----------------------------------------------------------------- #

#### Step 25 - Define eligible sites and build site-month matrix ####
# Sites with at least 3 years of coverage that includes observed non-masting
# months. Where a site appears in both monitoring and literature sources,
# monitoring takes priority.

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
      dataset_name %in% monitoring_datasets ~ 1,
      dataset_name == "Literature_review"   ~ 2
    )
  ) %>%
  group_by(site) %>%
  filter(priority == min(priority)) %>%
  ungroup()

# Masting indicator per site per month.
# TRUE = masting (FL or FR); FALSE = observed, not masting.
site_month_matrix <- site_dataset_priority %>%
  filter(phenophase %in% c("FL", "FR", "none")) %>%
  group_by(site, month_date) %>%
  summarise(
    is_masting = any(phenophase %in% c("FL", "FR")),
    .groups    = "drop"
  )


#### Step 26 - Build model data ####
# One row per site-month. Response: is_masting (0/1).
# Predictor: n_co_masting = number of OTHER sites masting within +/- 3 months.
# NOTE: slow, due to the rowwise calculation. Allow a few minutes.

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


#### Step 27 - Fit mixed effects logistic regression ####
# Fixed effect: n_co_masting. Random effect: site, controlling for
# differences in baseline masting frequency between sites.

glmer_model <- glmer(
  is_masting ~ n_co_masting + (1 | site),
  data   = model_data,
  family = binomial
)

summary(glmer_model)

glmer_results <- broom.mixed::tidy(
  glmer_model, conf.int = TRUE, exponentiate = TRUE
) %>%
  filter(effect == "fixed")

print(glmer_results)

null_model <- glmer(
  is_masting ~ 1 + (1 | site),
  data   = model_data,
  family = binomial
)

cat("\nModel comparison (null vs full):\n")
anova(null_model, glmer_model)


#### Step 28 - Figure: clustering model results ####
# Panel a: predicted probability curve with confidence ribbon.
# Panel b: observed proportions with fitted values at integer counts.

pred_data <- expand.grid(n_co_masting = seq(0, 7, by = 0.1)) %>%
  mutate(
    lp      = fixef(glmer_model)[1] + fixef(glmer_model)[2] * n_co_masting,
    lp_low  = (fixef(glmer_model)[1] - 1.96 * sqrt(vcov(glmer_model)[1,1])) +
      (fixef(glmer_model)[2] - 1.96 * sqrt(vcov(glmer_model)[2,2])) * n_co_masting,
    lp_high = (fixef(glmer_model)[1] + 1.96 * sqrt(vcov(glmer_model)[1,1])) +
      (fixef(glmer_model)[2] + 1.96 * sqrt(vcov(glmer_model)[2,2])) * n_co_masting,
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

p_clust_a <- ggplot(pred_data, aes(x = n_co_masting, y = prob)) +
  geom_ribbon(aes(ymin = prob_low, ymax = prob_high),
              fill = col_FL, alpha = 0.2) +
  geom_line(colour = col_FL, linewidth = 1.2) +
  scale_x_continuous(
    breaks = 0:7,
    name   = "Number of other sites masting within \u00b13 months"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name   = "Predicted probability of masting",
    limits = c(0, NA)
  ) +
  annotate(
    "text", x = 1, y = 0.75, hjust = 0, size = 3, colour = "grey30",
    label = paste0(
      "OR = ", round(exp(fixef(glmer_model)[2]), 2),
      " (95% CI: ", round(glmer_results$conf.low[2], 2), "\u2013",
      round(glmer_results$conf.high[2], 2), "), p < 0.001"
    )
  ) +
  theme_masting()

p_clust_b <- ggplot() +
  geom_col(data = obs_props, aes(x = n_co_masting, y = prop),
           fill = col_FL, alpha = 0.4, width = 0.6) +
  geom_line(data = pred_integer, aes(x = n_co_masting, y = prob),
            colour = col_FL, linewidth = 1.2) +
  geom_point(data = pred_integer, aes(x = n_co_masting, y = prob),
             colour = col_FL, size = 2.5) +
  geom_text(
    data = obs_props,
    aes(x = n_co_masting, y = prop + 0.01, label = paste0("n=", n_total)),
    size = 3, colour = "grey40", vjust = 0
  ) +
  scale_x_continuous(
    breaks = 0:7,
    name   = "Number of other sites masting within \u00b13 months"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name   = "Observed proportion masting",
    limits = c(0, NA)
  ) +
  theme_masting()

logistic_plot <- (p_clust_a / p_clust_b) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.tag    = element_text(size = 11, face = "bold"),
    plot.margin = margin(5, 10, 5, 5)
  )

logistic_plot

ggsave(
  file.path(output_dir, "temporal_clustering_logistic.png"),
  plot   = logistic_plot,
  width  = fig_width,
  height = fig_height_double,
  dpi    = fig_dpi
)


#### Step 29 - Clustering sensitivity analyses ####
# 29a. Minimum site coverage threshold (2, 3, 5 years).
# 29b. Co-masting window width (+/- 3 vs +/- 6 months).
# Both are slow. Each refits the model from scratch.

# 29a - Site inclusion threshold
run_clustering_model <- function(min_years) {

  cat("\n#### Running model with min_years =", min_years, "####\n")

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

  site_matrix <- site_dataset_priority %>%
    filter(site %in% eligible_sites,
           phenophase %in% c("FL", "FR", "none")) %>%
    group_by(site, month_date) %>%
    summarise(is_masting = any(phenophase %in% c("FL", "FR")),
              .groups = "drop")

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

  full_mod <- glmer(is_masting ~ n_co_masting + (1 | site),
                    data = model_dat, family = binomial)
  null_mod <- glmer(is_masting ~ 1 + (1 | site),
                    data = model_dat, family = binomial)

  or  <- exp(fixef(full_mod)["n_co_masting"])
  ci  <- exp(confint(full_mod, parm = "beta_", method = "Wald"))
  lrt <- anova(null_mod, full_mod)

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

coverage_sensitivity <- bind_rows(
  run_clustering_model(min_years = 2),
  run_clustering_model(min_years = 3),
  run_clustering_model(min_years = 5)
)

print(coverage_sensitivity, width = Inf)

write.csv(
  coverage_sensitivity,
  file.path(output_dir, "supplementary_coverage_sensitivity.csv"),
  row.names = FALSE
)


# 29b - Co-masting window width.
# Tests whether the +/- 3 month window under-detects co-masting between
# north and south, given the ~5 month offset in onset timing between them.
window_sensitivity <- purrr::map_dfr(c(3, 6), function(w) {

  md <- site_month_matrix %>%
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
            month_date >= focal_month - months(w),
            month_date <= focal_month + months(w),
            is_masting == TRUE
          ) %>%
          distinct(site) %>%
          nrow()
      }
    ) %>%
    ungroup()

  m  <- glmer(is_masting ~ n_co_masting + (1 | site),
              data = md, family = binomial)
  s  <- summary(m)$coefficients["n_co_masting", ]
  ci <- confint(m, parm = "n_co_masting", method = "Wald")

  tibble::tibble(
    window  = w,
    n_obs   = nrow(md),
    OR      = round(exp(s["Estimate"]), 3),
    ci_low  = round(exp(ci[1]), 3),
    ci_high = round(exp(ci[2]), 3),
    p       = s["Pr(>|z|)"]
  )
})

print(window_sensitivity, width = Inf)

write.csv(
  window_sensitivity,
  file.path(output_dir, "supplementary_window_sensitivity.csv"),
  row.names = FALSE
)


# ----------------------------------------------------------------- #
####          REGIONAL MASTING PERIOD IDENTIFICATION            ####
# ----------------------------------------------------------------- #

#### Step 30 - Justify the gap threshold ####
# Distribution of quiet period lengths, used to support the 6-month gap
# threshold that splits regional masting periods. Diagnostic only.

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

quiet_runs %>%
  summarise(min = min(run_len), median = median(run_len),
            mean = mean(run_len), max = max(run_len), n_runs = n())

quiet_runs %>% dplyr::count(run_len) %>% arrange(run_len) %>% print(n = 40)

ggplot(quiet_runs, aes(x = run_len)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "grey40") +
  geom_vline(xintercept = 6, linetype = "dashed", colour = "red") +
  labs(x = "Length of quiet period (months)", y = "Count") +
  theme_masting()


#### Step 31 - Monthly site activity counts ####

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
    fill = list(n_flowering = 0, n_fruiting = 0, n_none = 0,
                n_sites_active = 0, n_sites_observed = 0)
  )


#### Step 32 - Identify regional masting periods ####
# Consecutive months of activity are grouped into periods, split by gaps of
# 6 or more months with no site masting anywhere. Periods are retained if at
# least 3 sites were active at some point.
#
# These periods structure the timeline and onset map figures. They are NOT
# used to subset data for the seasonality, frequency, or clustering
# analyses, each of which uses its own site inclusion criteria.

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


#### Step 33 - Figure: regional masting timeline ####
# Bars above the axis: sites flowering and fruiting per month.
# Bars below the axis: sites observed but not masting.
# Shaded blocks: identified regional masting periods.

active_long <- active_by_month %>%
  dplyr::select(month_date, n_flowering, n_fruiting) %>%
  pivot_longer(cols = c(n_flowering, n_fruiting),
               names_to = "phenophase", values_to = "n_sites") %>%
  mutate(phenophase = recode(phenophase,
                             "n_flowering" = "Flowering (FL)",
                             "n_fruiting"  = "Fruiting (FR)"))

none_data <- active_by_month %>%
  dplyr::select(month_date, n_none) %>%
  mutate(n_sites = -n_none)

period_rects <- period_summary %>%
  mutate(xmin = period_start - days(15),
         xmax = period_end   + days(15))

y_max <- max(
  max(active_long$n_sites, na.rm = TRUE),
  max(abs(none_data$n_sites), na.rm = TRUE)
)

x_breaks_fixed <- seq(as.Date("1970-01-01"), as.Date("2030-01-01"),
                      by = "5 years")

date_from <- as.Date("1972-01-01")
date_to   <- as.Date("2027-01-01")

regional_timeline <- ggplot() +
  geom_rect(
    data = period_rects %>% filter(xmin <= date_to, xmax >= date_from),
    aes(xmin = pmax(xmin, date_from), xmax = pmin(xmax, date_to),
        ymin = -Inf, ymax = Inf, fill = "Regional masting period"),
    colour = NA
  ) +
  geom_col(
    data = none_data %>% filter(month_date >= date_from,
                                month_date <  date_to),
    aes(x = month_date, y = n_sites, fill = "Observed non-masting"),
    colour = NA, width = 31, alpha = 0.9
  ) +
  geom_col(
    data = active_long %>% filter(n_sites > 0,
                                  month_date >= date_from,
                                  month_date <  date_to),
    aes(x = month_date, y = n_sites, fill = phenophase),
    width = 31, colour = NA, alpha = 0.85, position = "stack"
  ) +
  geom_hline(yintercept = 0, colour = "black", linewidth = 0.4) +
  scale_fill_manual(
    name   = NULL,
    values = c(
      "Flowering (FL)"          = col_timeline_fl,
      "Fruiting (FR)"           = col_timeline_fr,
      "Observed non-masting"    = "grey55",
      "Regional masting period" = "#E8EEF4"
    ),
    breaks = c("Flowering (FL)", "Fruiting (FR)",
               "Observed non-masting", "Regional masting period")
  ) +
  scale_x_date(
    breaks = x_breaks_fixed, date_labels = "%Y",
    limits = c(date_from, date_to), expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    labels = function(x) abs(x),
    limits = c(-y_max, y_max),
    breaks = c(-10, -5, 0, 5, 10)
  ) +
  labs(x = "Year", y = "Number of sites active") +
  theme_masting() +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "grey60", fill = NA),
    legend.position  = "top",
    legend.direction = "horizontal"
  ) +
  guides(fill = guide_legend(title = NULL, nrow = 1))

regional_timeline

ggsave(
  file.path(output_dir, "regional_masting_timeline.png"),
  plot   = regional_timeline,
  width  = fig_width,
  height = 2.5,
  dpi    = fig_dpi
)

# ----------------------------------------------------------------- #
####        ONSET TIMING MAPS BY REGIONAL MASTING PERIOD         ####
# ----------------------------------------------------------------- #


#### Step 34 - Select periods to map ####
# Uses the same 11 regional masting periods identified in Step 32 and shown
# in Figure 3a, so the map (Figure 3b) and timeline are based on identical
# events. No additional duration/spatial-coverage filter is applied here.

mast_period_data <- all_harmonised %>%
  filter(!is.na(month_date), phenophase %in% c("FL", "FR")) %>%
  left_join(
    mast_timeline %>%
      filter(!is.na(period_id)) %>%
      dplyr::select(month_date, period_id),
    by = "month_date"
  ) %>%
  filter(!is.na(period_id))

keep_periods   <- period_summary$period_id
period_ids_chr <- as.character(sort(keep_periods))

cat("Masting periods mapped:", length(keep_periods), "\n")


#### Step 35 - Site onset per period ####

site_onset <- mast_period_data %>%
  filter(period_id %in% keep_periods) %>%
  group_by(period_id, site) %>%
  summarise(onset_date = min(month_date), .groups = "drop") %>%
  left_join(site_coords, by = "site") %>%
  filter(!is.na(lat), !is.na(long)) %>%
  group_by(period_id) %>%
  arrange(onset_date, .by_group = TRUE) %>%
  mutate(
    onset_rank   = dense_rank(onset_date),
    onset_scaled = scales::rescale(onset_rank, to = c(0, 1),
                                   from = range(onset_rank, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(period_id = factor(as.character(period_id), levels = period_ids_chr))

# Facet labels showing the date range of each period
facet_labs <- site_onset %>%
  group_by(period_id) %>%
  summarise(
    onset_min = format(min(onset_date), "%b %Y"),
    onset_max = format(max(onset_date), "%b %Y"),
    .groups   = "drop"
  ) %>%
  mutate(facet_lab = paste0(onset_min, " \u2192 ", onset_max)) %>%
  { setNames(.$facet_lab, .$period_id) }

# One colour per period
event_cols <- c(
  "1"  = "#FF4E00", "4"  = "#F5D547", "5"  = "#E8A33D",
  "7"  = "#B1693B", "9"  = "#8C5A2E", "10" = "#33673B",
  "11" = "#638475", "13" = "#026B79", "16" = "#7B506F",
  "18" = "#1F1A38", "22" = "#6AD5CB"
)

#### Step 36 - Figure: faceted onset timing map ####
# One panel per regional masting period.
# Colour = period identity; opacity = onset order (darker = earlier).

onset_map <- ggplot() +
  geom_sf(data = sea_coast, colour = "grey60", linewidth = 0.3) +
  geom_point(
    data = site_onset,
    aes(x = long, y = lat, colour = period_id, alpha = onset_scaled),
    size = 3.5
  ) +
  facet_wrap(~ period_id, labeller = labeller(period_id = facet_labs)) +
  coord_sf(xlim = c(95, 120), ylim = c(-5, 10), expand = FALSE) +
  scale_colour_manual(values = event_cols, guide = "none") +
  scale_alpha_continuous(
    name     = "Onset timing (early \u2192 late)",
    range    = c(0.95, 0.35),
    breaks   = c(0, 0.25, 0.5, 0.75, 1.0),
    na.value = 0.25
  ) +
  scale_x_continuous(breaks = c(100, 110, 120)) +
  scale_y_continuous(breaks = c(-5, 0, 5, 10)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_masting() +
  theme(
    panel.grid       = element_blank(),
    panel.border     = element_rect(colour = "grey80", fill = NA,
                                    linewidth = 0.3),
    strip.background = element_blank(),
    strip.text       = element_text(size = 8, colour = "grey30"),
    axis.text        = element_text(size = 7),
    legend.position  = "bottom",
    legend.direction = "horizontal",
    plot.margin      = margin(5, 15, 5, 5)
  ) +
  guides(alpha = guide_legend(title.position = "top",
                              direction = "horizontal", nrow = 1))

onset_map

ggsave(
  file.path(output_dir, "masting_onset_map.png"),
  plot   = onset_map,
  width  = fig_width,
  height = fig_height_triple,
  dpi    = fig_dpi
)


############################################### - #
# For questions or issues with this code and
# accompanying data, please contact Zoë Lieb at
# zoelieb1@gmail.com
############################################### - #
