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


#' Read CSV or R binary (.rda/.RData) file
#'
#' @param fn         File path
#' @param rda_object Object name when fn is .rda; NULL = first data frame found
#' @return data.frame
#' @keywords internal
read_input <- function(fn, rda_object = NULL) {
  if (!file.exists(fn)) stop("File not found: ", fn)
  ext <- tolower(tools::file_ext(fn))
  if (ext == "csv") {
    return(read.csv(fn, stringsAsFactors = FALSE))
  }
  if (ext %in% c("rda", "rdata")) {
    e <- new.env(parent = emptyenv())
    load(fn, envir = e)
    if (!is.null(rda_object)) {
      if (!exists(rda_object, envir = e))
        stop("Object '", rda_object, "' not found in ", fn)
      return(as.data.frame(get(rda_object, envir = e)))
    }
    objs <- ls(e)
    dfs  <- Filter(function(nm) is.data.frame(get(nm, envir = e)), objs)
    if (length(dfs) == 0) stop("No data frames found in ", fn)
    if (length(dfs) > 1)
      message("Multiple data frames in ", fn, "; using '", dfs[1],
              "'. Set rda_object= to choose.")
    return(as.data.frame(get(dfs[1], envir = e)))
  }
  stop("Unsupported file type '", ext, "'. Use .csv, .rda, or .RData.")
}
