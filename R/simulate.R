# ============================================================================
# simulate.R — field trial simulation and raw data visualisation
# ============================================================================


#' Compute spatial effect for a bench grid
#'
#' @param grid              data.frame with Row and Col columns
#' @param n_rows            number of rows in bench
#' @param n_cols            number of columns in bench
#' @param spatial_type      integer 1-7:
#'   1 = no spatial effect;
#'   2 = row gradient (half-sine);
#'   3 = column gradient (half-sine);
#'   4 = row + column gradient (additive);
#'   5 = central patch (Gaussian spotlight);
#'   6 = edge effects (elevated at borders with diagonal asymmetry);
#'   7 = localised hotspot (Gaussian bell at an off-centre position)
#' @param spatial_intensity peak magnitude of spatial effect in phenotype units
#' @param spatial_scale     bandwidth multiplier (default 1)
#' @return numeric vector of length nrow(grid)
#' @keywords internal
compute_spatial <- function(grid, n_rows, n_cols, spatial_type,
                             spatial_intensity, spatial_scale = 1) {
  cr <- (n_rows + 1) / 2
  cc <- (n_cols + 1) / 2

  switch(as.character(spatial_type),

    "1" = rep(0, nrow(grid)),

    "2" = {
      raw <- spatial_intensity * sin(pi / 2 * (grid$Row - 1) * spatial_scale / (n_rows - 1))
      raw - mean(raw)
    },

    "3" = {
      raw <- spatial_intensity * sin(pi / 2 * (grid$Col - 1) * spatial_scale / (n_cols - 1))
      raw - mean(raw)
    },

    "4" = {
      row_effect <- sin(pi / 2 * (grid$Row - 1) * spatial_scale / (n_rows - 1))
      col_effect <- sin(pi / 2 * (grid$Col - 1) * spatial_scale / (n_cols - 1))
      raw <- spatial_intensity * (row_effect + col_effect) / 2
      raw - mean(raw)
    },

    "5" = {
      sigma_r <- (n_rows / 4) * spatial_scale
      sigma_c <- (n_cols / 4) * spatial_scale
      raw <- spatial_intensity * exp(-(
        (grid$Row - cr)^2 / (2 * sigma_r^2) +
        (grid$Col - cc)^2 / (2 * sigma_c^2)
      ))
      raw - mean(raw)
    },

    "6" = {
      d_row  <- pmin(grid$Row - 1, n_rows - grid$Row) / ((n_rows - 1) / 2)
      d_col  <- pmin(grid$Col - 1, n_cols - grid$Col) / ((n_cols - 1) / 2)
      d_edge <- pmin(d_row, d_col)

      edge_effect <- (1 - d_edge^spatial_scale) * 0.6
      diag_effect <- ((n_rows - grid$Row) / (n_rows - 1) +
                      (n_cols - grid$Col) / (n_cols - 1)) / 2 * 0.4

      raw <- spatial_intensity * (edge_effect + diag_effect)
      raw - mean(raw)
    },

    "7" = {
      cr_hot <- 0.35 * (n_rows - 1) + 1
      cc_hot <- 0.65 * (n_cols - 1) + 1

      d_row  <- (grid$Row - cr_hot) / n_rows
      d_col  <- (grid$Col - cc_hot) / n_cols
      sigma  <- spatial_scale * 0.17

      raw <- spatial_intensity * exp(-(d_row^2 + d_col^2) / (2 * sigma^2))
      raw - mean(raw)
    },

    stop("spatial_type must be 1-7: 1 (none), 2 (row gradient), 3 (col gradient), ",
         "4 (row + col gradient), 5 (central patch), 6 (edge effects), 7 (hotspot)")
  )
}


#' Simulate a multi-bench unreplicated field trial dataset
#'
#' Generates `n_rows * n_cols` genotypes, each appearing exactly once per bench.
#' Genotypes are randomly assigned to plots within each bench independently.
#'
#' The observed phenotype for each plot is:
#' `observed = true_genotype_effect + spatial_effect + N(0, sd_error^2)`
#'
#' When `save_csv = TRUE`, two CSV files are written to `output_dir`:
#' `<pheno_name>_simulation.csv` and `<pheno_name>_true_effects.csv`.
#'
#' @param n_bench           number of benches/fields
#' @param n_rows            rows per bench
#' @param n_cols            columns per bench; `n_rows * n_cols` genotypes simulated
#' @param mean_pheno        grand mean of the phenotype
#' @param sd_geno           SD of true genotype effects
#' @param spatial_type      integer 1-7; see [compute_spatial()] for details
#' @param spatial_intensity peak magnitude of the spatial gradient in phenotype units
#' @param spatial_scale     bandwidth multiplier (default 1)
#' @param sd_error          plot-level random noise SD
#' @param pheno_min         lower clip applied to true genotype effects
#' @param pheno_max         upper clip applied to true genotype effects
#' @param pheno_name        name of the phenotype column in output
#' @param spatial_intensity_per_bench numeric vector of length `n_bench` (optional);
#'   per-bench spatial intensity. `NULL` = use scalar `spatial_intensity`
#' @param spatial_type_per_bench integer vector of length `n_bench` (optional);
#'   per-bench spatial type. `NULL` = use scalar `spatial_type`
#' @param save_csv          if `TRUE` (default), write CSV files to `output_dir`
#' @param output_dir        directory for CSV output; created if absent
#' @param seed              random seed for reproducibility
#'
#' @return list with elements: `$data`, `$true_effects`, `$pheno_name`
#'
#' @examples
#' sim <- simulate_field_trial(n_bench = 2, n_rows = 8, n_cols = 10,
#'                              save_csv = FALSE, seed = 1)
#' head(sim$data)
#'
#' @export
simulate_field_trial <- function(n_bench           = 4,
                                  n_rows            = 10,
                                  n_cols            = 30,
                                  mean_pheno        = 55,
                                  sd_geno           = 12,
                                  spatial_type      = 5,
                                  spatial_intensity = 30,
                                  spatial_scale     = 1,
                                  sd_error          = 4,
                                  pheno_min         = 0,
                                  pheno_max         = 100,
                                  pheno_name        = "BNI",
                                  spatial_intensity_per_bench = NULL,
                                  spatial_type_per_bench      = NULL,
                                  save_csv          = TRUE,
                                  output_dir        = "data",
                                  seed              = 42) {
  set.seed(seed)

  n_geno    <- n_rows * n_cols
  genotypes <- paste0("G", sprintf(paste0("%0", nchar(n_geno), "d"), 1:n_geno))

  geno_true <- pmax(pheno_min,
                    pmin(pheno_max,
                         rnorm(n_geno, mean = mean_pheno, sd = sd_geno)))

  bench_list <- lapply(seq_len(n_bench), function(b) {
    grid          <- expand.grid(Row = 1:n_rows, Col = 1:n_cols)
    grid$Bench    <- paste0("Bench", b)
    grid$Genotype <- sample(genotypes)

    b_type      <- if (!is.null(spatial_type_per_bench))
                     spatial_type_per_bench[b] else spatial_type
    b_intensity <- if (!is.null(spatial_intensity_per_bench))
                     spatial_intensity_per_bench[b] else spatial_intensity

    grid$True_Spatial <- compute_spatial(grid, n_rows, n_cols,
                                         b_type, b_intensity, spatial_scale)

    grid[[pheno_name]] <- geno_true[match(grid$Genotype, genotypes)] +
                          grid$True_Spatial +
                          rnorm(nrow(grid), 0, sd_error)

    grid[, c("Bench", "Genotype", "Row", "Col", "True_Spatial", pheno_name)]
  })

  sim_df           <- do.call(rbind, bench_list)
  rownames(sim_df) <- NULL

  true_effects <- data.frame(Genotype = genotypes, True = geno_true,
                              stringsAsFactors = FALSE)
  colnames(true_effects)[2] <- paste0("True_", pheno_name)

  if (save_csv) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    data_file    <- file.path(output_dir, paste0(pheno_name, "_simulation.csv"))
    effects_file <- file.path(output_dir, paste0(pheno_name, "_true_effects.csv"))
    write.csv(sim_df,       data_file,    row.names = FALSE)
    write.csv(true_effects, effects_file, row.names = FALSE)
    cat("Data written to:         ", data_file, "\n")
    cat("True effects written to: ", effects_file, "\n")
  }

  list(data = sim_df, true_effects = true_effects, pheno_name = pheno_name)
}


#' Plot a field heatmap for a single bench
#'
#' @param bench_dfr   data.frame for one bench, or the full dataset when `bench`
#'   is specified
#' @param pheno       name of the response column
#' @param rn          name of the row column
#' @param cn          name of the column column
#' @param gt          name of the genotype column (used when `show_labels = TRUE`)
#' @param bench_col   name of the bench column; used only when `bench` is not NULL
#' @param bench       bench identifier to subset from `bench_dfr`; if NULL (default)
#'   `bench_dfr` is used as-is (must already be a single bench)
#' @param bench_label character label for the plot title
#' @param show_labels logical; if TRUE, genotype names are printed inside tiles
#' @param shared_limits numeric c(min, max) for fill scale; NULL = auto from bench
#' @param label_size  text size for genotype labels
#' @return ggplot object
#'
#' @export
plot_bench_heatmap <- function(bench_dfr,
                                pheno,
                                rn            = "Row",
                                cn            = "Col",
                                gt            = "Genotype",
                                bench_col     = "Bench",
                                bench         = NULL,
                                bench_label   = "",
                                show_labels   = FALSE,
                                shared_limits = NULL,
                                label_size    = 2.5) {
  if (!is.null(bench)) {
    bench_dfr   <- bench_dfr[bench_dfr[[bench_col]] == bench, ]
    if (bench_label == "") bench_label <- as.character(bench)
  }

  fill_limits <- if (!is.null(shared_limits)) shared_limits else
    range(bench_dfr[[pheno]], na.rm = TRUE)

  p <- ggplot(bench_dfr,
              aes(x = .data[[cn]], y = .data[[rn]], fill = .data[[pheno]])) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(
      name     = pheno,
      limits   = fill_limits,
      option   = "viridis",
      na.value = "grey70"
    ) +
    scale_y_reverse(breaks = sort(unique(bench_dfr[[rn]]))) +
    scale_x_continuous(breaks = sort(unique(bench_dfr[[cn]]))) +
    labs(title = bench_label, x = cn, y = rn) +
    theme_bw() +
    theme(
      plot.title   = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 8),
      axis.text    = element_text(size = 7),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7),
      panel.grid   = element_blank()
    )

  if (show_labels && gt %in% colnames(bench_dfr)) {
    p <- p + geom_text(aes(label = .data[[gt]]),
                       size     = label_size,
                       colour   = "white",
                       fontface = "bold")
  }

  p
}


#' Plot heatmaps for all benches and save as PNG
#'
#' @param dfr          data.frame with all benches
#' @param pheno        name of the response column
#' @param rn           name of the row column
#' @param cn           name of the column column
#' @param gt           name of the genotype column
#' @param bench_col    name of the bench column
#' @param show_labels  logical; show genotype names inside tiles
#' @param shared_scale logical; shared fill scale across all benches
#' @param output_png   file path for PNG output
#' @param ncol_panels  number of bench panels per row
#' @param panel_width  width per panel in inches
#' @param panel_height height per panel in inches
#' @param label_size   text size for genotype labels
#' @return invisibly, the patchwork ggplot object
#'
#' @export
plot_all_bench_heatmaps <- function(dfr,
                                     pheno,
                                     rn            = "Row",
                                     cn            = "Col",
                                     gt            = "Genotype",
                                     bench_col     = "Bench",
                                     show_labels   = FALSE,
                                     shared_scale  = TRUE,
                                     output_png    = "heatmaps.png",
                                     ncol_panels   = 2,
                                     panel_width   = 7,
                                     panel_height  = 4,
                                     label_size    = 2.5) {
  benches       <- unique(dfr[[bench_col]])
  shared_limits <- if (shared_scale) range(dfr[[pheno]], na.rm = TRUE) else NULL

  plots <- lapply(benches, function(b) {
    plot_bench_heatmap(
      bench_dfr     = dfr[dfr[[bench_col]] == b, ],
      pheno         = pheno,
      rn            = rn,
      cn            = cn,
      gt            = gt,
      bench_label   = as.character(b),
      show_labels   = show_labels,
      shared_limits = shared_limits,
      label_size    = label_size
    )
  })

  n_plots     <- length(plots)
  nrow_panels <- ceiling(n_plots / ncol_panels)

  composite <- wrap_plots(plots, ncol = ncol_panels) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title    = paste0("Raw observed values: ", pheno),
      subtitle = if (shared_scale) "Shared colour scale across benches" else
                   "Independent colour scale per bench",
      theme    = theme(
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey40")
      )
    )

  total_w <- panel_width  * min(n_plots, ncol_panels)
  total_h <- panel_height * nrow_panels + 0.5

  ggplot2::ggsave(output_png, composite, width = total_w, height = total_h, dpi = 150)
  cat("Heatmap saved to:", output_png, "\n")

  invisible(composite)
}
