# spagam

GAM-based spatial correction for field trials.

Fits [`mgcv`](https://cran.r-project.org/package=mgcv) tensor-product smooths to remove spatial trends from unreplicated greenhouse or field trial data, then returns corrected genotype estimates (BLUEs and/or BLUPs), spatial trend surfaces, diagnostic plots, and model summaries.

## Installation

```r
# install.packages("remotes")
remotes::install_github("alexd106/spagam")
```

## Usage

```r
library(spagam)

# Read your data
df <- read.csv("trial_data.csv")

# Run spatial correction
out <- correct_spatial(
  data        = df,
  pheno_cols  = c("yield", "height"),
  geno_col    = "Genotype",
  row_col     = "Row",
  col_col     = "Col",
  bench_col   = "Bench",   # NULL for single-bench
  output_dir  = "output/gam"
)

# Results
out$blues          # corrected genotype estimates
out$spatial_trends # per-plot spatial surface values
out$model_summary  # EDF, residual SD, convergence
```

Outputs (CSV and PNG) are written to `output_dir`.

## Simulation

```r
sim <- simulate_field_trial(n_bench = 3, n_rows = 10, n_cols = 20,
                             spatial_type = 5, seed = 42)
head(sim$data)
```

## License

MIT
