# Southeast Asian Masting Phenology

Data and analysis code accompanying:

> Lieb, Z. et al. (in review). [Paper title]. *Ecology Letters*.

## Contents

- `data/sea_phenology_combined_data.csv` — the full processed dataset used
  in all analyses: one row per site, month, and event type (flowering/
  fruiting), synthesised from long-term monitoring, expert observations,
  and a systematic literature review across 26 sites in Southeast Asia.
- `data/site_data_summary.csv` — a site-level summary table (one row per
  site x contributing data source) documenting observation periods,
  dipterocarp-only status, onset-detection threshold methods, and source
  citations for each site.
- `metadata/sea_phenology_combined_data_metadata.csv` — column definitions
  and units for the combined dataset.
- `metadata/site_data_summary_metadata.csv` — column definitions and units
  for the site summary table.
- `analysis/sea_pheno_analysis_cleaned.R` — full analysis script reproducing
  all figures, tables, and statistical results reported in the paper.

## Reproducing the analysis

The analysis script reads `data/sea_phenology_combined_data.csv` directly.
Required R packages: tidyverse, lubridate, ggrepel, geosphere,
rnaturalearth, rnaturalearthdata, sf, circular, patchwork, scales,
colorspace, lme4, broom.mixed.

```r
# from the repo root
source("analysis/sea_pheno_analysis_cleaned.R")
```

## Data availability note

The underlying per-dataset raw data (individual long-term monitoring
datasets, expert observations, and literature review records) are not
included here, as much of this data is privately held by third-party
contributors. This repository contains only the final, processed,
site-month dataset used in the published analyses, consistent with
Ecology Letters' data-sharing requirements.

## Citation

If you use this data or code, please cite:

> [Full paper citation once available]
>
> Dataset DOI: [Zenodo DOI once minted]

## License

- Code: MIT License (see `LICENSE`)
- Data: CC-BY 4.0

## Contact

Zoë Lieb — [email]
