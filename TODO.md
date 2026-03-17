# spagam TODO

_Last updated: 2026-03-17

---

## Phase 1 — Package foundation

- [x] Scaffold package structure (`R/`, `tests/`, `vignettes/`, `inst/extdata/`)
- [x] Split `spatial_correct_gam.R` into modules: `input.R`, `diagnostics.R`,
      `fit-single.R`, `fit-joint.R`, `correct.R`
- [x] Migrate `simulate_spatial_data.R` → `simulate.R`
- [x] Rename `run_spatial_gam()` → `correct_spatial()`
- [x] Copy example data to `inst/extdata/BNI_simulation.csv`
- [x] Initialise git repo: `git init && git add . && git commit -m "initial scaffold"`
- [x] Create GitHub repo and push
- [x] Run ``devtools::document() to generate `NAMESPACE` and `man/` pages
- [x] Run `devtools::check()` and resolve any warnings or notes
- [x] Drop `.rda`/`.RData` input support; `correct_spatial()` now accepts a data frame directly
- [x] Remove `read_input()` and stale `man/read_input.Rd`
- [x] Fix unresolved roxygen links to internal functions (`adaptive_nseg`, `empirical_semivariogram`, `compute_spatial`)
- [x] Create `data/wheatdata.rda` from `inst/extdata/wheatdata.csv`

---

## Phase 2 — Tests

- [x] `test-input.R` — `adaptive_nseg()` bounds (n < 10, n = 40, n = 100)
- [x] `test-input.R` — `correct_spatial()` rejects non-data-frame and zero-row input
- [ ] `test-simulate.R` — `simulate_field_trial()` dimensions, column names,
      reproducibility via seed, `save_csv = FALSE` path
- [ ] `test-simulate.R` — all 7 spatial types run without error
- [ ] `test-fit-single.R` — `fit_mgcv_bench()` returns correct list structure
      for BLUEs, BLUPs, and both on small simulated data
- [ ] `test-fit-joint.R` — `fit_mgcv_joint()` returns correct list structure on
      two-bench simulated data
- [ ] `test-correct.R` — `correct_spatial()` single-bench: output CSVs exist,
      correct dimensions, no errors
- [ ] `test-correct.R` — `correct_spatial()` multi-bench: same checks
- [ ] `test-correct.R` — BLUPs + bench_col raises informative error

---

## Phase 3 — S3 classes

- [ ] Define `spatial_fit` S3 class for `correct_spatial()` return value
- [ ] `print.spatial_fit()` — compact one-line summary
- [ ] `summary.spatial_fit()` — model summary table
- [ ] `plot.spatial_fit()` — delegate to diagnostic plot functions
- [ ] Define `spatial_simulation` class for `simulate_field_trial()` return
- [ ] `print.spatial_simulation()`, `summary.spatial_simulation()`,
      `plot.spatial_simulation()`

---

## Phase 4 — Vignettes

- [ ] Getting-started vignette: simulate data → `correct_spatial()` → inspect outputs
- [ ] Adapt `spatial_correct_gam_guide.Rmd` into a full package vignette

---

## Phase 5 — GxE extension (future)

- [ ] Design API for multi-environment GxE modelling
- [ ] Implement joint GxE GAM model
- [ ] Extend `spatial_fit` class or define `gxe_fit`
