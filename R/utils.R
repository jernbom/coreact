#' Resolve X/Y Pipeline Arguments (Internal)
#'
#' Handles the logic for arguments that can be provided as a single vector (applied to both)
#' or a list of length 2 (specific to X and Y).
#'
#' @param val The input argument (list or vector).
#' @param arg_name String. Name of the argument for error reporting.
#' @return A list with elements `$x` and `$y`.
#' @keywords internal
resolve_xy_args <- function(val, arg_name) {
  if (is.list(val)) {
    if (length(val) == 1) return(list(x = val[[1]], y = val[[1]]))
    if (length(val) == 2) return(list(x = val[[1]], y = val[[2]]))
    stop(sprintf("Argument '%s' provided as a list must be length 1 or 2.", arg_name))
  } else {
    # If atomic vector (or NULL), recycle the same value for both X and Y
    return(list(x = val, y = val))
  }
}

#' Filter Coreact Data by Prevalence (Internal)
#'
#' Filters a coreact_data object based on row prevalence.
#' Automatically detects if the threshold is a percentage (< 1) or absolute count (>= 1).
#'
#' @param obj A `coreact_data` object.
#' @param thresh Numeric. The filter threshold.
#' @return A subsetted `coreact_data` object.
#' @keywords internal
filter_by_prevalence <- function(obj, thresh) {
  if (thresh <= 0) return(obj)

  # Determine cutoff: Absolute vs Percentage
  # If 0 < thresh < 1: Treat as percentage of total samples.
  # If thresh >= 1: Treat as absolute count.
  if (thresh < 1) {
    cutoff <- ceiling(thresh * ncol(obj$mat))
  } else {
    cutoff <- thresh
  }

  prev <- Matrix::rowSums(obj$mat > 0)
  keep <- prev >= cutoff

  if (sum(keep) == 0) {
    stop(sprintf("Filter removed all features from %s (Cutoff: %s).", obj$name, cutoff))
  }

  # Uses the S3 subsetting method '[.coreact_data' defined in class_coreact.R
  return(obj[keep, ])
}
