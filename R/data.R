#' Wheat yield trial data (Gilmour, Cullis & Verbyla 1997)
#'
#' Grain yield from a replicated wheat trial in South Australia. Originally
#' distributed in the **SpATS** package (GPL licence).
#'
#' @format A data frame with 330 rows and 7 columns:
#' \describe{
#'   \item{yield}{Grain yield (numeric)}
#'   \item{geno}{Genotype identifier}
#'   \item{rep}{Replicate block}
#'   \item{row}{Row coordinate (numeric)}
#'   \item{col}{Column coordinate (numeric)}
#'   \item{rowcode}{Row blocking code}
#'   \item{colcode}{Column blocking code}
#' }
#'
#' @details
#' The field layout is 22 rows x 15 columns, with 107 genotypes and 3 replicates.
#' The file is available via:
#' ```r
#' system.file("extdata", "wheatdata.csv", package = "spagam")
#' ```
#'
#' @source
#' Gilmour, A.R., Cullis, B.R., and Verbyla, A.P. (1997). Accounting for
#' natural and extraneous variation in the analysis of field experiments.
#' *Journal of Agricultural, Biological, and Environmental Statistics*,
#' **2**, 269–293.
#'
#' Originally distributed in the SpATS R package
#' (<https://cran.r-project.org/package=SpATS>) under the GPL licence.
"wheatdata"
