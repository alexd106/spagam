# ============================================================================
# input.R — data reading and grid utilities
# ============================================================================


#' Choose number of B-spline basis functions adaptively
#'
#' Rule: roughly half the unique positions, bounded to \[5, 20\].
#'
#' @param n_unique Number of unique positions along one axis
#' @return integer
#' @keywords internal
adaptive_nseg <- function(n_unique) {
  as.integer(min(max(5L, floor(n_unique / 2L)), 20L))
}
