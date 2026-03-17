#' spagam: GAM-Based Spatial Correction for Field Trials
#'
#' @description
#' Fits mgcv-based GAM spatial correction models for unreplicated field trials.
#' Supports single-bench designs (BLUEs and BLUPs) and joint multi-bench designs
#' (BLUEs only). Includes simulation utilities and diagnostic plotting.
#'
#' Main entry points:
#' - [correct_spatial()] — fit model, write outputs, return results
#' - [simulate_field_trial()] — generate a realistic simulated trial dataset
#' - [plot_all_bench_heatmaps()] — visualise raw observations across benches
#'
#' @keywords internal
"_PACKAGE"

#' @import ggplot2
#' @import patchwork
#' @importFrom mgcv gam
NULL
