# Southeast Asian Masting Phenology

Data and analysis code accompanying:

> Lieb, Z. et al. (in review). [Paper title]. *Nature Plants*.
>
> ****** Duplication error identified after journal submission.
> New datasheets and code to resolve are included in this repo,
> but this is not yet reflected in manuscript text.

## Contents

- `sea_phenology_combined_data.csv` — the full processed dataset used
  in all analyses: one row per site, month, and event type (flowering/
  fruiting), and data source, synthesised from long-term monitoring,
  expert observations, and a systematic literature review across 26
  sites in Southeast Asia. Where multiple sources report on the same
  site-month, all are retained as separate rows rather than collapsed;
  downstream analyses select among them as needed.

- `sea_masting_literature_review_events.csv` — masting events
  extracted from a systematic literature review of general flowering and
  mast fruiting research in Southeast Asia (Lieb et al. 2026). One row per
  reported flowering or fruiting event, with citation, site, and reported
  onset/end dates. Used directly in the analysis to identify onset months
  and event durations, since explicit reported dates are more reliable
  here than inferring onsets from consecutive monthly observations. A
  separate, site-month-expanded version of this same literature data is
  merged into `sea_phenology_combined_data.csv` alongside all other
  sources.

- `site_data_summary.csv` — a site-level summary table (one row per
  site x contributing data source) documenting observation periods,
  dipterocarp-only status, onset-detection threshold methods, and source
  citations for each site.

- `sea_phenology_combined_data_metadata.csv` — column definitions
  and units for the combined dataset.

- `site_data_summary_metadata.csv` — column definitions and units
  for the site summary table.

- `sea_pheno_analysis_zel_public_04082026.R` — full analysis script reproducing
  all figures (and additional visuals), tables, and statistical results reported in the paper.

## Reproducing the analysis

The analysis script reads `sea_phenology_combined_data.csv` and
`sea_masting_literature_review_events.csv` directly. Figures and
supplementary tables are written to an `output/` folder, created
automatically if it does not exist.

Required R packages: tidyverse, lubridate, ggrepel, geosphere,
rnaturalearth, rnaturalearthdata, sf, circular, diptest, patchwork, scales,
colorspace, lme4, broom.mixed, scico.

```r
# from the repo root
source("sea_pheno_analysis_public.R")
```
## Data sources and terms explanation

Defining terms in the dataset

- Data source types:
Lit review: the evidence for this site-month of a masting event was found in literature review (Lieb et al. 2026). Short citations are provided for each site in the site_data_summary.csv. Some sites have multiple citations. Full references of these citations are provided in the section below.

Full citations of lit review sources are included in sea_masting_literature_review_events.csv, column “citation”

- Published long term: refers to datasets sourced from literature that include raw tree reproduction proportions. Even though data may have been published, there was consultation with data providers to determine how to correctly interpret the dataset.

These include:
Numata_Satake_2022 – data published and additional material provided by Drs Shinya Numata and Akiko Satake for FRIM. (Numata et al. 2022)
Chen_published_2017 – published dataset by Dr Yu-Yun Chen for Pasoh Forest Reserve. (Chen et al. 2018)

- Provided data: This includes datasets that were made available by contributors but are otherwise unpublished. Or not publicly available

These include:

YIARI_provided_2025 – Data provided for Gunung Tarak and Bukit Baka Bukit Raya National Park by organization, Yayasan Inisiasi Alam Rehabilitasi Indonesia (YIARI), thanks to Gail Campbell-Smith and Karmele Llano Sánchez.
Marshall_provided_2025 – Data provided for Gunung Palung thanks to Andrew Marshall

- Additional contributed obs: refers to observations contributed by researchers based on their personal knowledge of a site, rather that from a dataset. These observations were provided via a monthly datasheet and contributors could indicate months of flowering, fruiting, no activity, or non-observed months.
Contributors for this section include: 
John Payne, Tetsuro Hosaka, Lisa Ong, Michael O’Brien, Shawn Lum, and Andrew Marshall


## Data availability note

The underlying per-dataset raw data for long-term monitoring sites and
expert observations are not included here, as much of this data is
privately held by third-party contributors. This repository contains the
final, processed site-month dataset used in the published analyses, along
with the literature review event records, consistent with Ecology
Letters' data-sharing requirements.

## Citation

If you use this data or code, please cite:

> [Full paper citation once available]
>
> Dataset DOI: [![DOI](https://zenodo.org/badge/1313426754.svg)](https://doi.org/10.5281/zenodo.22632580)

## License

- Code: MIT License (see `LICENSE`)
- Data: CC-BY 4.0

## Contact

Zoë Lieb — zoelieb1@gmail.com
