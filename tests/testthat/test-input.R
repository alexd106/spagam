test_that("adaptive_nseg() clamps to lower bound for small n", {
  # n < 10: floor(6/2) = 3, clamped up to 5
  expect_equal(spagam:::adaptive_nseg(6L), 5L)
  expect_equal(spagam:::adaptive_nseg(1L), 5L)
})

test_that("adaptive_nseg() is within bounds for moderate n", {
  # n = 40: floor(40/2) = 20, exactly at upper bound
  expect_equal(spagam:::adaptive_nseg(40L), 20L)
})

test_that("adaptive_nseg() clamps to upper bound for large n", {
  # n = 100: floor(100/2) = 50, clamped down to 20
  expect_equal(spagam:::adaptive_nseg(100L), 20L)
})

test_that("adaptive_nseg() returns an integer", {
  expect_type(spagam:::adaptive_nseg(20L), "integer")
})

# ---------------------------------------------------------------------------
# correct_spatial() input validation
# ---------------------------------------------------------------------------

test_that("correct_spatial() rejects non-data-frame input", {
  expect_error(
    correct_spatial(data = list(a = 1), pheno_cols = "BNI"),
    "`data` must be a data frame"
  )
  expect_error(
    correct_spatial(data = matrix(1:4, 2, 2), pheno_cols = "BNI"),
    "`data` must be a data frame"
  )
})

test_that("correct_spatial() rejects zero-row data frame", {
  empty <- simulate_field_trial(n_bench = 1, n_rows = 6, n_cols = 6,
                                save_csv = FALSE, seed = 1)$data[0, ]
  expect_error(
    correct_spatial(
      data       = empty,
      pheno_cols = "BNI",
      geno_col   = "Genotype",
      row_col    = "Row",
      col_col    = "Col",
      output_dir = tempfile()
    ),
    "`data` has zero rows"
  )
})
