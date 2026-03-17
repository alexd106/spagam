# ============================================================================
# diagnostics.R — semivariogram, heatmap, and diagnostic plot helpers
# ============================================================================


#' Compute empirical semivariogram
#'
#' @param data      data.frame with row, col, and value columns
#' @param row_col   Row coordinate column name
#' @param col_col   Column coordinate column name
#' @param value_col Value column name
#' @param n_bins    Number of distance bins
#' @return data.frame(dist, sv) or NULL when < 10 non-NA observations
#' @keywords internal
empirical_semivariogram <- function(data, row_col, col_col, value_col,
                                    n_bins = 15) {
  ok  <- !is.na(data[[value_col]])
  d   <- data[ok, ]
  if (nrow(d) < 10) return(NULL)

  rows <- d[[row_col]]
  cols <- d[[col_col]]
  vals <- d[[value_col]]

  row_diff <- outer(rows, rows, "-")
  col_diff <- outer(cols, cols, "-")
  dists    <- sqrt(row_diff^2 + col_diff^2)
  sv_mat   <- (outer(vals, vals, "-"))^2 / 2

  up     <- upper.tri(dists)
  d_vec  <- dists[up]
  sv_vec <- sv_mat[up]

  max_d  <- quantile(d_vec, 0.6)
  keep   <- d_vec <= max_d & d_vec > 0
  d_vec  <- d_vec[keep]
  sv_vec <- sv_vec[keep]
  if (length(d_vec) == 0) return(NULL)

  breaks <- seq(0, max(d_vec), length.out = n_bins + 1)
  bin    <- cut(d_vec, breaks = breaks, include.lowest = TRUE)
  agg    <- aggregate(cbind(sv = sv_vec, dist = d_vec) ~ bin, FUN = mean)
  agg[order(agg$dist), c("dist", "sv")]
}


#' Plot a field heatmap
#'
#' @param data      data.frame with coordinates and fill column
#' @param row_col   Row coordinate column name
#' @param col_col   Column coordinate column name
#' @param value_col Fill column name
#' @param title     Plot title
#' @param limits    c(min, max) for fill scale; NULL = auto
#' @return ggplot object
#' @keywords internal
plot_heatmap <- function(data, row_col, col_col, value_col,
                          title = "", limits = NULL) {
  ok_rows <- !is.na(data[[row_col]]) & !is.na(data[[col_col]])
  d       <- data[ok_rows, ]
  if (is.null(limits)) limits <- range(d[[value_col]], na.rm = TRUE)

  ggplot(d, aes(x = .data[[col_col]], y = .data[[row_col]],
                fill = .data[[value_col]])) +
    geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_viridis_c(
      name     = value_col,
      limits   = limits,
      option   = "viridis",
      na.value = "grey70"
    ) +
    scale_y_reverse(breaks = sort(unique(d[[row_col]]))) +
    scale_x_continuous(breaks = sort(unique(d[[col_col]]))) +
    labs(title = title, x = col_col, y = row_col) +
    theme_bw(base_size = 8) +
    theme(
      plot.title   = element_text(size = 9, face = "bold", hjust = 0.5),
      panel.grid   = element_blank(),
      axis.text    = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.text  = element_text(size = 6)
    )
}


#' Plot empirical semivariogram
#'
#' @param vgram data.frame(dist, sv) from [empirical_semivariogram()]
#' @param title Plot title
#' @return ggplot object
#' @keywords internal
plot_variogram <- function(vgram, title = "") {
  if (is.null(vgram) || nrow(vgram) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "Insufficient data",
                 size = 3, colour = "grey50") +
        theme_void() +
        labs(title = title) +
        theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    )
  }
  ggplot(vgram, aes(x = dist, y = sv)) +
    geom_point(colour = "#2c7bb6", size = 1.5) +
    geom_line(colour = "#2c7bb6", linewidth = 0.5) +
    labs(title = title, x = "Distance (plot units)", y = "Semivariance") +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 7),
      axis.text  = element_text(size = 6)
    )
}


#' Plot row and column random effect magnitudes
#'
#' Returns a two-panel (row | col) strip plot showing the estimated random
#' effect for each row and column level. Useful for diagnosing systematic
#' edge or lane effects that the smooth surface cannot capture.
#'
#' @param bench_data  data.frame for the bench
#' @param result      Return value from fit_mgcv_bench() or fit_mgcv_joint()
#' @param row_col     Row coordinate column name
#' @param col_col     Column coordinate column name
#' @param bench_label Character label for title
#' @return patchwork ggplot or NULL if no RE data available
#' @keywords internal
plot_row_col_re <- function(bench_data, result, row_col, col_col,
                            bench_label = "") {
  if (is.null(result$row_re) && is.null(result$col_re)) return(NULL)

  plots <- list()

  if (!is.null(result$row_re)) {
    rd <- data.frame(level = bench_data[[row_col]], effect = result$row_re)
    rd <- aggregate(effect ~ level, data = rd, FUN = mean)
    plots$row <- ggplot(rd, aes(x = factor(level), y = effect)) +
      geom_col(fill = "#2c7bb6", alpha = 0.7, width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      labs(title = "Row random effects", x = row_col, y = "Effect") +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 7),
            axis.text  = element_text(size = 6))
  }

  if (!is.null(result$col_re)) {
    cd <- data.frame(level = bench_data[[col_col]], effect = result$col_re)
    cd <- aggregate(effect ~ level, data = cd, FUN = mean)
    plots$col <- ggplot(cd, aes(x = factor(level), y = effect)) +
      geom_col(fill = "#d7191c", alpha = 0.7, width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      labs(title = "Column random effects", x = col_col, y = "Effect") +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
            axis.title = element_text(size = 7),
            axis.text  = element_text(size = 6))
  }

  if (length(plots) == 0) return(NULL)

  page_title <- paste0("Row / column random effects",
                       if (bench_label != "") paste0("  |  ", bench_label) else "")

  wrap_plots(plots, ncol = length(plots)) +
    plot_annotation(
      title = page_title,
      theme = theme(plot.title = element_text(size = 11, face = "bold"))
    )
}


#' Compose 6-panel diagnostic plot for one bench x trait
#'
#' Panels: observed | spatial trend | residuals (top row);
#' variogram raw | variogram residuals | QQ residuals (bottom row).
#'
#' @param bench_data  data.frame for the bench (includes observed pheno)
#' @param result      Return value from fit_mgcv_bench() or fit_mgcv_joint()
#' @param pheno       Phenotype column name in bench_data
#' @param row_col     Row coordinate column name
#' @param col_col     Column coordinate column name
#' @param bench_label Character label for title
#' @return patchwork ggplot
#' @keywords internal
make_diagnostic_plots <- function(bench_data, result, pheno,
                                   row_col, col_col, bench_label = "") {
  pd <- bench_data[, c(row_col, col_col), drop = FALSE]
  pd$obs     <- bench_data[[pheno]]
  pd$spatial <- if (!is.null(result$spatial_smooth)) result$spatial_smooth
                else if (!is.null(result$spatial))   result$spatial
                else NA_real_
  pd$resid   <- if (!is.null(result$residuals)) result$residuals else NA_real_

  safe_range <- function(...) {
    r <- range(..., na.rm = TRUE)
    if (any(!is.finite(r))) c(0, 1) else r
  }

  p1 <- plot_heatmap(pd, row_col, col_col, "obs",
                     title = "Observed", limits = safe_range(pd$obs))
  p2 <- plot_heatmap(pd, row_col, col_col, "spatial",
                     title = "Spatial trend", limits = safe_range(pd$spatial))
  p3 <- plot_heatmap(pd, row_col, col_col, "resid",
                     title = "Residuals", limits = safe_range(pd$resid))

  vg_obs  <- empirical_semivariogram(pd, row_col, col_col, "obs")
  vg_res  <- empirical_semivariogram(pd, row_col, col_col, "resid")
  p4 <- plot_variogram(vg_obs, title = "Variogram: observed")
  p5 <- plot_variogram(vg_res, title = "Variogram: residuals")

  qq_df <- data.frame(resid = pd$resid[!is.na(pd$resid)])
  p6 <- if (nrow(qq_df) > 2) {
    ggplot(qq_df, aes(sample = resid)) +
      stat_qq(size = 1, colour = "#2c7bb6") +
      stat_qq_line(colour = "red", linewidth = 0.5) +
      labs(title = "QQ: corrected residuals", x = "Theoretical", y = "Sample") +
      theme_bw(base_size = 8) +
      theme(
        plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 7),
        axis.text  = element_text(size = 6)
      )
  } else {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient data",
               size = 3, colour = "grey50") +
      theme_void() +
      labs(title = "QQ: corrected residuals") +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  }

  page_title <- paste0(pheno,
                       if (bench_label != "") paste0("  |  ", bench_label) else "")

  (p1 | p2 | p3) / (p4 | p5 | p6) +
    plot_annotation(
      title = page_title,
      theme = theme(plot.title = element_text(size = 11, face = "bold"))
    )
}
