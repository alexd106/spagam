# ============================================================================
# fit-single.R — single-bench mgcv GAM fitter
# ============================================================================


#' mgcv te_ps_re spatial correction for a single bench
#'
#' Formula (BLUEs):
#' `pheno ~ 0 + Genotype + te(row, col, bs='ps') + s(row_f, bs='re') + s(col_f, bs='re')`
#'
#' @param bench_data    data.frame for one bench
#' @param pheno         Phenotype column name
#' @param geno_col      Genotype column name
#' @param row_col       Row coordinate column name
#' @param col_col       Column coordinate column name
#' @param k_row         Basis dimension for rows (NULL = auto)
#' @param k_col         Basis dimension for columns (NULL = auto)
#' @param estimate_type `"BLUEs"`, `"BLUPs"`, or `"both"`
#' @return list with elements: blues, blups, fitted, residuals, spatial_smooth,
#'   spatial_total, row_re, col_re, converged, edf_spatial
#' @keywords internal
fit_mgcv_bench <- function(bench_data, pheno, geno_col, row_col, col_col,
                            k_row = NULL, k_col = NULL,
                            estimate_type = "BLUEs") {
  res <- list(blues = NULL, blups = NULL, fitted = NULL,
              residuals = NULL, spatial_smooth = NULL, spatial_total = NULL,
              row_re = NULL, col_re = NULL,
              converged = NA, edf_spatial = NA)

  bench_data[[geno_col]] <- as.factor(bench_data[[geno_col]])
  bench_data$row_f       <- as.factor(bench_data[[row_col]])
  bench_data$col_f       <- as.factor(bench_data[[col_col]])
  geno_levels            <- levels(bench_data[[geno_col]])

  n_row <- length(unique(bench_data[[row_col]]))
  n_col <- length(unique(bench_data[[col_col]]))
  kr    <- if (!is.null(k_row)) k_row else adaptive_nseg(n_row)
  kc    <- if (!is.null(k_col)) k_col else adaptive_nseg(n_col)

  te_term  <- paste0("te(", row_col, ", ", col_col,
                     ", bs=c('ps','ps'), k=c(", kr, ",", kc, "))")
  re_extra <- " + s(row_f, bs='re') + s(col_f, bs='re')"

  # ---- BLUEs -----------------------------------------------------------------
  if (estimate_type %in% c("BLUEs", "both")) {
    fm <- as.formula(paste0(pheno, " ~ 0 + ", geno_col, " + ", te_term, re_extra))
    m  <- tryCatch(
      gam(fm, data = bench_data, method = "REML"),
      error = function(e) { warning("mgcv BLUEs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      res$converged <- m$converged
      if (!isTRUE(m$converged))
        warning("GAM did not converge for BLUEs (", pheno, ")")

      coef_names <- paste0(geno_col, geno_levels)
      blues_vals <- coef(m)[coef_names]
      vcov_diag  <- diag(vcov(m))
      blues_se   <- sqrt(vcov_diag[coef_names])

      res$blues <- data.frame(
        Genotype = geno_levels,
        BLUE     = as.numeric(blues_vals),
        BLUE_SE  = as.numeric(blues_se),
        stringsAsFactors = FALSE
      )
      names(res$blues)[1] <- geno_col

      terms_pred <- predict(m, type = "terms")
      te_col     <- grep("^te\\(", colnames(terms_pred), value = TRUE)
      row_re_col <- grep("^s\\(row_f\\)", colnames(terms_pred), value = TRUE)
      col_re_col <- grep("^s\\(col_f\\)", colnames(terms_pred), value = TRUE)

      smooth_val <- as.numeric(terms_pred[, te_col[1]])
      row_re_val <- if (length(row_re_col)) as.numeric(terms_pred[, row_re_col[1]]) else rep(0, nrow(bench_data))
      col_re_val <- if (length(col_re_col)) as.numeric(terms_pred[, col_re_col[1]]) else rep(0, nrow(bench_data))

      res$spatial_smooth <- smooth_val
      res$spatial_total  <- smooth_val + row_re_val + col_re_val
      res$row_re         <- row_re_val
      res$col_re         <- col_re_val
      res$fitted         <- as.numeric(fitted(m))
      res$residuals      <- as.numeric(residuals(m))
      res$edf_spatial    <- sum(m$edf[grep("^te\\(", names(m$edf))])
    }
  }

  # ---- BLUPs -----------------------------------------------------------------
  if (estimate_type %in% c("BLUPs", "both")) {
    fm <- as.formula(paste0(pheno, " ~ s(", geno_col, ", bs='re') + ",
                            te_term, re_extra))
    m  <- tryCatch(
      gam(fm, data = bench_data, method = "REML"),
      error = function(e) { warning("mgcv BLUPs error: ", e$message); NULL }
    )
    if (!is.null(m)) {
      if (is.na(res$converged)) res$converged <- m$converged
      if (!isTRUE(m$converged))
        warning("GAM did not converge for BLUPs (", pheno, ")")

      re_pattern <- paste0("s\\(", geno_col, "\\)")
      re_idx     <- grep(re_pattern, names(coef(m)))
      blup_vals  <- coef(m)[re_idx]
      intercept  <- coef(m)["(Intercept)"]
      blup_total <- as.numeric(intercept) + as.numeric(blup_vals)

      mean_row <- mean(bench_data[[row_col]], na.rm = TRUE)
      mean_col <- mean(bench_data[[col_col]], na.rm = TRUE)
      newdat   <- setNames(
        data.frame(geno_levels, mean_row, mean_col),
        c(geno_col, row_col, col_col)
      )
      newdat[[geno_col]] <- factor(newdat[[geno_col]], levels = geno_levels)
      newdat$row_f <- factor(round(mean_row), levels = levels(bench_data$row_f))
      newdat$col_f <- factor(round(mean_col), levels = levels(bench_data$col_f))
      tp <- tryCatch(
        predict(m, newdata = newdat, type = "terms", se.fit = TRUE),
        error = function(e) NULL
      )
      re_col_name <- if (!is.null(tp))
        grep(re_pattern, colnames(tp$fit), value = TRUE)[1] else NA
      blup_se <- if (!is.null(tp) && !is.na(re_col_name))
        as.numeric(tp$se.fit[, re_col_name]) else rep(NA_real_, length(geno_levels))

      res$blups <- data.frame(
        Genotype = geno_levels,
        BLUP     = blup_total,
        BLUP_SE  = blup_se,
        stringsAsFactors = FALSE
      )
      names(res$blups)[1] <- geno_col

      if (is.null(res$fitted)) {
        terms_pred <- predict(m, type = "terms")
        te_col     <- grep("^te\\(", colnames(terms_pred), value = TRUE)
        row_re_col <- grep("^s\\(row_f\\)", colnames(terms_pred), value = TRUE)
        col_re_col <- grep("^s\\(col_f\\)", colnames(terms_pred), value = TRUE)

        smooth_val <- as.numeric(terms_pred[, te_col[1]])
        row_re_val <- if (length(row_re_col)) as.numeric(terms_pred[, row_re_col[1]]) else rep(0, nrow(bench_data))
        col_re_val <- if (length(col_re_col)) as.numeric(terms_pred[, col_re_col[1]]) else rep(0, nrow(bench_data))

        res$spatial_smooth <- smooth_val
        res$spatial_total  <- smooth_val + row_re_val + col_re_val
        res$row_re         <- row_re_val
        res$col_re         <- col_re_val
        res$fitted         <- as.numeric(fitted(m))
        res$residuals      <- as.numeric(residuals(m))
        if (is.na(res$edf_spatial))
          res$edf_spatial <- sum(m$edf[grep("^te\\(", names(m$edf))])
      }
    }
  }

  res
}
