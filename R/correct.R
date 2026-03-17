# ============================================================================
# correct.R — main user-facing entry point
# ============================================================================


#' Run mgcv-based spatial correction end to end
#'
#' Fits the spatial GAM (single- or multi-bench), writes CSV outputs and
#' diagnostic plots to `output_dir`, and returns the key results invisibly.
#'
#' Single-bench mode fits one GAM per trait and supports both BLUEs and BLUPs.
#' Multi-bench mode fits a joint GAM across all benches and supports BLUEs only.
#'
#' @param data          A data frame containing phenotype, genotype, and
#'   coordinate columns.
#' @param pheno_cols    Character vector of phenotype column names
#' @param output_dir    Output directory; created automatically if absent
#' @param geno_col      Genotype identifier column name
#' @param row_col       Row coordinate column name (must be numeric)
#' @param col_col       Column coordinate column name (must be numeric)
#' @param bench_col     Bench column name; `NULL` triggers the single-bench model
#' @param estimate_type One of `"BLUEs"` (default), `"BLUPs"`, or `"both"`.
#'   BLUPs and `"both"` require single-bench mode (`bench_col = NULL`).
#' @param k_row         Basis dimension for rows; `NULL` = auto (roughly half the unique row positions, bounded to \[5, 20\])
#' @param k_col         Basis dimension for columns; `NULL` = auto
#'
#' @return Invisibly: list with elements `blues`, `blups`, `spatial_trends`,
#'   `model_summary`
#'
#' @examples
#' \dontrun{
#' df  <- read.csv(system.file("extdata", "BNI_simulation.csv", package = "spagam"))
#' out <- correct_spatial(
#'   data        = df,
#'   pheno_cols  = "BNI",
#'   geno_col    = "Genotype",
#'   row_col     = "Row",
#'   col_col     = "Col",
#'   bench_col   = "Bench",
#'   output_dir  = tempdir()
#' )
#' }
#'
#' @export
correct_spatial <- function(
  data,
  pheno_cols,
  output_dir    = "output/gam",
  geno_col      = "geno",
  row_col       = "row",
  col_col       = "col",
  bench_col     = NULL,
  estimate_type = "BLUEs",
  k_row         = NULL,
  k_col         = NULL
) {

  # --------------------------------------------------------------------------
  # Input validation
  # --------------------------------------------------------------------------
  if (!is.data.frame(data))
    stop("`data` must be a data frame, not ", class(data)[1], ".", call. = FALSE)
  if (nrow(data) == 0)
    stop("`data` has zero rows.", call. = FALSE)
  required_cols <- c(geno_col, row_col, col_col, pheno_cols)
  if (!is.null(bench_col)) required_cols <- c(required_cols, bench_col)
  missing_cols  <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop("Required column(s) not found in data: ",
         paste(missing_cols, collapse = ", "), "\n",
         "  Available columns: ", paste(names(data), collapse = ", "))

  for (coord in c(row_col, col_col)) {
    if (!is.numeric(data[[coord]]))
      stop("Column '", coord, "' must be numeric (found: ",
           class(data[[coord]]), ")")
    if (anyNA(data[[coord]]))
      stop("Column '", coord, "' contains NA values; remove or impute before running.")
  }

  for (ph in pheno_cols) {
    if (!is.numeric(data[[ph]]))
      stop("Phenotype column '", ph, "' must be numeric.")
  }

  n_unique_row <- length(unique(data[[row_col]]))
  n_unique_col <- length(unique(data[[col_col]]))
  if (n_unique_row < 5)
    warning("Only ", n_unique_row, " unique row values -- spatial surface may not be identifiable.")
  if (n_unique_col < 5)
    warning("Only ", n_unique_col, " unique column values -- spatial surface may not be identifiable.")

  is_multibench <- !is.null(bench_col)
  if (is_multibench) {
    bench_ids <- sort(unique(data[[bench_col]]))
    n_bench   <- length(bench_ids)
  } else {
    n_bench   <- 1L
    bench_ids <- "single"
  }

  if (!is.null(bench_col) && estimate_type %in% c("BLUPs", "both")) {
    stop("Multi-bench BLUPs are not yet supported. ",
         "Use estimate_type = \"BLUEs\" with bench_col, ",
         "or use single-bench mode (bench_col = NULL) for BLUPs.",
         call. = FALSE)
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # --------------------------------------------------------------------------
  # Model fitting
  # --------------------------------------------------------------------------
  all_blues_list <- list()
  all_blups_list <- list()
  spatial_rows   <- list()
  summary_rows   <- list()

  for (pheno in pheno_cols) {
    if (!is_multibench) {
      # Single-bench path
      result <- fit_mgcv_bench(
        bench_data    = data,
        pheno         = pheno,
        geno_col      = geno_col,
        row_col       = row_col,
        col_col       = col_col,
        k_row         = k_row,
        k_col         = k_col,
        estimate_type = estimate_type
      )

      if (!is.null(result$blues)) {
        df <- result$blues
        names(df)[names(df) == "BLUE"]    <- pheno
        names(df)[names(df) == "BLUE_SE"] <- paste0(pheno, "_SE")
        all_blues_list[[pheno]] <- df
      }
      if (!is.null(result$blups)) {
        df <- result$blups
        names(df)[names(df) == "BLUP"]    <- pheno
        names(df)[names(df) == "BLUP_SE"] <- paste0(pheno, "_SE")
        all_blups_list[[pheno]] <- df
      }

      if (!is.null(result$spatial_smooth)) {
        spatial_rows[[length(spatial_rows) + 1]] <- data.frame(
          data[, c(geno_col, row_col, col_col), drop = FALSE],
          bench          = "single",
          pheno          = pheno,
          spatial_smooth = result$spatial_smooth,
          spatial_total  = result$spatial_total,
          stringsAsFactors = FALSE
        )
      }

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        pheno       = pheno,
        bench       = "single",
        resid_sd    = if (!is.null(result$residuals)) sd(result$residuals, na.rm = TRUE) else NA,
        edf_spatial = result$edf_spatial,
        converged   = result$converged,
        stringsAsFactors = FALSE
      )

      diag_plot <- tryCatch(
        make_diagnostic_plots(data, result, pheno, row_col, col_col, ""),
        error = function(e) { warning("Diagnostic plot failed: ", e$message); NULL }
      )
      if (!is.null(diag_plot)) {
        out_png <- file.path(output_dir, paste0("diagnostics_", pheno, "_single.png"))
        ggplot2::ggsave(out_png, diag_plot, width = 12, height = 8, dpi = 150)
      }

      re_plot <- tryCatch(
        plot_row_col_re(data, result, row_col, col_col, ""),
        error = function(e) { warning("RE plot failed: ", e$message); NULL }
      )
      if (!is.null(re_plot)) {
        out_png <- file.path(output_dir, paste0("row_col_re_", pheno, "_single.png"))
        ggplot2::ggsave(out_png, re_plot, width = 8, height = 4, dpi = 150)
      }

    } else {
      # Multi-bench path
      result <- fit_mgcv_joint(
        data      = data,
        pheno     = pheno,
        geno_col  = geno_col,
        row_col   = row_col,
        col_col   = col_col,
        bench_col = bench_col,
        k_row     = k_row,
        k_col     = k_col
      )

      if (!is.null(result$blues)) {
        df <- result$blues
        names(df)[names(df) == "BLUE"]    <- pheno
        names(df)[names(df) == "BLUE_SE"] <- paste0(pheno, "_SE")
        all_blues_list[[pheno]] <- df
      }

      if (!is.null(result$spatial_smooth)) {
        spatial_rows[[length(spatial_rows) + 1]] <- data.frame(
          data[, c(geno_col, row_col, col_col, bench_col), drop = FALSE],
          pheno          = pheno,
          spatial_smooth = result$spatial_smooth,
          spatial_total  = result$spatial_total,
          stringsAsFactors = FALSE
        )
      }

      for (b in bench_ids) {
        b_mask  <- data[[bench_col]] == b
        b_resid <- if (!is.null(result$residuals)) result$residuals[b_mask] else NULL
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          pheno       = pheno,
          bench       = as.character(b),
          resid_sd    = if (!is.null(b_resid)) sd(b_resid, na.rm = TRUE) else NA,
          edf_spatial = result$edf_spatial,
          converged   = result$converged,
          stringsAsFactors = FALSE
        )
      }

      for (b in bench_ids) {
        b_mask  <- data[[bench_col]] == b
        b_data  <- data[b_mask, , drop = FALSE]
        b_result <- list(
          spatial_smooth = if (!is.null(result$spatial_smooth)) result$spatial_smooth[b_mask] else NULL,
          residuals      = if (!is.null(result$residuals))      result$residuals[b_mask]      else NULL,
          fitted         = if (!is.null(result$fitted))         result$fitted[b_mask]         else NULL,
          row_re         = if (!is.null(result$row_re))         result$row_re[b_mask]         else NULL,
          col_re         = if (!is.null(result$col_re))         result$col_re[b_mask]         else NULL
        )
        diag_plot <- tryCatch(
          make_diagnostic_plots(b_data, b_result, pheno, row_col, col_col,
                                paste0("Bench: ", b)),
          error = function(e) { warning("Diagnostic plot failed (bench ", b, "): ",
                                        e$message); NULL }
        )
        if (!is.null(diag_plot)) {
          safe_bench <- gsub("[^A-Za-z0-9_.-]", "_", as.character(b))
          out_png <- file.path(output_dir,
                               paste0("diagnostics_", pheno, "_", safe_bench, ".png"))
          ggplot2::ggsave(out_png, diag_plot, width = 12, height = 8, dpi = 150)
        }

        re_plot <- tryCatch(
          plot_row_col_re(b_data, b_result, row_col, col_col,
                          paste0("Bench: ", b)),
          error = function(e) { warning("RE plot failed (bench ", b, "): ",
                                        e$message); NULL }
        )
        if (!is.null(re_plot)) {
          out_png <- file.path(output_dir,
                               paste0("row_col_re_", pheno, "_", safe_bench, ".png"))
          ggplot2::ggsave(out_png, re_plot, width = 8, height = 4, dpi = 150)
        }
      }
    }
  }

  # --------------------------------------------------------------------------
  # CSV outputs
  # --------------------------------------------------------------------------
  blues_df <- NULL
  blups_df <- NULL
  sp_df    <- NULL
  sum_df   <- NULL

  if (length(all_blues_list) > 0) {
    blues_df <- Reduce(function(a, b) merge(a, b, by = geno_col, all = TRUE),
                       all_blues_list)
    out_path <- file.path(output_dir, "BLUEs.csv")
    write.csv(blues_df, out_path, row.names = FALSE)
  }

  if (length(all_blups_list) > 0) {
    blups_df <- Reduce(function(a, b) merge(a, b, by = geno_col, all = TRUE),
                       all_blups_list)
    out_path <- file.path(output_dir, "BLUPs.csv")
    write.csv(blups_df, out_path, row.names = FALSE)
  }

  if (length(spatial_rows) > 0) {
    sp_df    <- do.call(rbind, spatial_rows)
    out_path <- file.path(output_dir, "spatial_trends.csv")
    write.csv(sp_df, out_path, row.names = FALSE)
  }

  if (length(summary_rows) > 0) {
    sum_df   <- do.call(rbind, summary_rows)
    out_path <- file.path(output_dir, "model_summary.csv")
    write.csv(sum_df, out_path, row.names = FALSE)
  }

  # --------------------------------------------------------------------------
  # Summary plots
  # --------------------------------------------------------------------------
  if (length(all_blues_list) > 0) {
    blues_long <- do.call(rbind, lapply(names(all_blues_list), function(ph) {
      df <- all_blues_list[[ph]]
      data.frame(trait = ph, value = df[[ph]], stringsAsFactors = FALSE)
    }))
    blues_long <- blues_long[!is.na(blues_long$value), ]

    p_dens <- ggplot(blues_long, aes(x = value)) +
      geom_histogram(aes(y = after_stat(density)),
                     bins = 30, fill = "#2c7bb6", alpha = 0.6) +
      geom_density(colour = "red", linewidth = 0.7) +
      facet_wrap(~ trait, scales = "free") +
      labs(title = "BLUE distribution per trait",
           x = "BLUE value", y = "Density") +
      theme_bw(base_size = 10) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))

    out_png <- file.path(output_dir, "blues_distribution.png")
    ggplot2::ggsave(out_png, p_dens,
                    width = max(6, 4 * length(all_blues_list)), height = 5, dpi = 150)
  }

  if (!is.null(sp_df) && !is.null(sum_df)) {
    sp_plots <- lapply(seq_len(nrow(sum_df)), function(i) {
      ph <- sum_df$pheno[i]
      b  <- sum_df$bench[i]
      if (is_multibench) {
        sub <- sp_df[sp_df$pheno == ph & sp_df[[bench_col]] == b, ]
      } else {
        sub <- sp_df[sp_df$pheno == ph, ]
      }
      if (nrow(sub) == 0) return(NULL)
      sub_plot <- sub[, c(row_col, col_col, "spatial_smooth"), drop = FALSE]
      plot_heatmap(sub_plot, row_col, col_col, "spatial_smooth",
                   title = paste0(ph, " | ", b))
    })
    sp_plots <- Filter(Negate(is.null), sp_plots)

    if (length(sp_plots) > 0) {
      n_col_grid <- min(length(sp_plots), 3L)
      composite  <- wrap_plots(sp_plots, ncol = n_col_grid) +
        plot_annotation(
          title = "Spatial smooth surfaces (te only)",
          theme = theme(plot.title = element_text(size = 13, face = "bold"))
        )
      out_png     <- file.path(output_dir, "spatial_surfaces.png")
      n_rows_grid <- ceiling(length(sp_plots) / n_col_grid)
      ggplot2::ggsave(out_png, composite,
                      width  = 5 * n_col_grid,
                      height = 4 * n_rows_grid,
                      dpi = 150)
    }
  }

  invisible(list(
    blues          = blues_df,
    blups          = blups_df,
    spatial_trends = sp_df,
    model_summary  = sum_df
  ))
}
