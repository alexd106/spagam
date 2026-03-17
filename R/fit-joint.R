# ============================================================================
# fit-joint.R — joint multi-bench mgcv GAM fitter
# ============================================================================


#' Joint multi-bench mgcv::gam spatial correction
#'
#' Fits a single GAM across all benches:
#' `pheno ~ 0 + Genotype + bench_f + te(row, col, bs='ps', by=bench_f) + s(row_f, bench_f, bs='re') + s(col_f, bench_f, bs='re')`
#'
#' Only BLUEs are supported on the multi-bench path.
#'
#' @param data      Full data.frame (all benches)
#' @param pheno     Phenotype column name
#' @param geno_col  Genotype column name
#' @param row_col   Row coordinate column name
#' @param col_col   Column coordinate column name
#' @param bench_col Bench column name
#' @param k_row     Basis dimension for rows (NULL = auto)
#' @param k_col     Basis dimension for columns (NULL = auto)
#' @return list with elements: blues, fitted, residuals, spatial_smooth,
#'   spatial_total, row_re, col_re, converged, edf_spatial
#' @keywords internal
fit_mgcv_joint <- function(data, pheno, geno_col, row_col, col_col, bench_col,
                            k_row = NULL, k_col = NULL) {
  res <- list(blues = NULL, fitted = NULL, residuals = NULL,
              spatial_smooth = NULL, spatial_total = NULL,
              row_re = NULL, col_re = NULL,
              converged = NA, edf_spatial = NA)

  data[[geno_col]] <- as.factor(data[[geno_col]])
  data$bench_f     <- as.factor(data[[bench_col]])
  data$row_f       <- as.factor(data[[row_col]])
  data$col_f       <- as.factor(data[[col_col]])
  geno_levels      <- levels(data[[geno_col]])
  bench_levels     <- levels(data$bench_f)

  nr_vec <- sapply(bench_levels, function(b)
    length(unique(data[[row_col]][data$bench_f == b])))
  nc_vec <- sapply(bench_levels, function(b)
    length(unique(data[[col_col]][data$bench_f == b])))
  kr <- if (!is.null(k_row)) k_row else adaptive_nseg(min(nr_vec))
  kc <- if (!is.null(k_col)) k_col else adaptive_nseg(min(nc_vec))

  fm <- as.formula(paste0(
    pheno, " ~ 0 + ", geno_col, " + bench_f + ",
    "te(", row_col, ", ", col_col,
    ", bs=c('ps','ps'), k=c(", kr, ",", kc, "), by=bench_f) + ",
    "s(row_f, bench_f, bs='re') + s(col_f, bench_f, bs='re')"
  ))

  m <- tryCatch(
    gam(fm, data = data, method = "REML"),
    error = function(e) { warning("mgcv joint error: ", e$message); NULL }
  )
  if (is.null(m)) return(res)

  res$converged <- m$converged
  if (!isTRUE(m$converged)) warning("Joint GAM did not converge (", pheno, ")")

  coef_names  <- paste0(geno_col, geno_levels)
  blues_vals  <- coef(m)[coef_names]
  vcov_diag   <- diag(vcov(m))
  blues_se    <- sqrt(vcov_diag[coef_names])

  res$blues <- data.frame(
    Genotype = geno_levels,
    BLUE     = as.numeric(blues_vals),
    BLUE_SE  = as.numeric(blues_se),
    stringsAsFactors = FALSE
  )
  names(res$blues)[1] <- geno_col

  res$fitted    <- as.numeric(fitted(m))
  res$residuals <- as.numeric(residuals(m))

  terms_pred  <- predict(m, type = "terms")
  te_cols     <- grep("^te\\(", colnames(terms_pred), value = TRUE)
  row_re_cols <- grep("^s\\(row_f", colnames(terms_pred), value = TRUE)
  col_re_cols <- grep("^s\\(col_f", colnames(terms_pred), value = TRUE)

  smooth_val <- if (length(te_cols) > 0)
    as.numeric(rowSums(terms_pred[, te_cols, drop = FALSE]))
  else
    rep(0, nrow(data))
  row_re_val <- if (length(row_re_cols) > 0)
    as.numeric(rowSums(terms_pred[, row_re_cols, drop = FALSE]))
  else
    rep(0, nrow(data))
  col_re_val <- if (length(col_re_cols) > 0)
    as.numeric(rowSums(terms_pred[, col_re_cols, drop = FALSE]))
  else
    rep(0, nrow(data))

  res$spatial_smooth <- smooth_val
  res$spatial_total  <- smooth_val + row_re_val + col_re_val
  res$row_re         <- row_re_val
  res$col_re         <- col_re_val
  res$edf_spatial    <- sum(m$edf[grep("^te\\(", names(m$edf))])

  res
}
