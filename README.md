# Southeast Asian Masting Phenology

Data and analysis code accompanying:

> Lieb, Z. et al. (in review). [Paper title]. *Ecology Letters*.

## Contents

- `data/sea_phenology_combined_data.csv` — the full processed dataset used
  in all analyses: one row per site, month, and event type (flowering/
  fruiting), synthesised from long-term monitoring, expert observations,
  and a systematic literature review across 26 sites in Southeast Asia.

- `data/sea_masting_literature_review_events.csv` — masting events
  extracted from a systematic literature review of general flowering and
  mast fruiting research in Southeast Asia (Lieb et al. 2026). One row per
  reported flowering or fruiting event, with citation, site, and reported
  onset/end dates. Used directly in the analysis to identify onset months
  and event durations, since explicit reported dates are more reliable
  here than inferring onsets from consecutive monthly observations. A
  separate, site-month-expanded version of this same literature data is
  merged into `sea_phenology_combined_data.csv` alongside all other
  sources.

- `data/site_data_summary.csv` — a site-level summary table (one row per
  site x contributing data source) documenting observation periods,
  dipterocarp-only status, onset-detection threshold methods, and source
  citations for each site.

- `metadata/sea_phenology_combined_data_metadata.csv` — column definitions
  and units for the combined dataset.

- `metadata/site_data_summary_metadata.csv` — column definitions and units
  for the site summary table.

- `analysis/sea_pheno_analysis_public.R` — full analysis script reproducing
  all figures, tables, and statistical results reported in the paper.

## Reproducing the analysis

The analysis script reads `data/sea_phenology_combined_data.csv` and
`data/sea_masting_literature_review_events.csv` directly. Figures and
supplementary tables are written to an `output/` folder, created
automatically if it does not exist.

Required R packages: tidyverse, lubridate, ggrepel, geosphere,
rnaturalearth, rnaturalearthdata, sf, circular, diptest, patchwork, scales,
colorspace, lme4, broom.mixed, scico.

```r
# from the repo root
source("analysis/sea_pheno_analysis_public.R")
```

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
> Dataset DOI: [Zenodo DOI once minted]

## License

- Code: MIT License (see `LICENSE`)
- Data: CC-BY 4.0

## Contact

Zoë Lieb — zoelieb1@gmail.com
