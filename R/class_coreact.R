#' @importFrom methods new
#' @importFrom Matrix Matrix rowSums
NULL

#' Constructor for coreact data
#'
#' Creates a structured object for antibody co-occurrence analysis.
#' Pre-calculates prevalence (row sums) to save compute time later.
#'
#' @param mat A matrix (binary/integer). Rows are features, columns are samples.
#' @param meta A data.frame containing feature annotations. Rows must match mat.
#' @param name A string identifier for this dataset (e.g., "viral", "human").
#'
#' @return An object of class 'coreact_data'
#' @export
new_coreact_data <- function(mat, meta, name = "unknown") {

  if (!is.matrix(mat) && !inherits(mat, "Matrix")) {
    stop("Input 'mat' must be a matrix or sparse Matrix.")
  }
  if (nrow(mat) != nrow(meta)) {
    stop(sprintf("Dimension mismatch: Matrix has %d rows but Metadata has %d rows.",
                 nrow(mat), nrow(meta)))
  }

  # Pre-calculate prevalence (rowSums on sparse matrix is efficient)
  prev <- Matrix::rowSums(mat, na.rm = TRUE)

  structure(
    list(
      mat = mat,
      meta = meta,
      prevalence = prev,
      n_samples = ncol(mat),
      n_features = nrow(mat),
      name = name
    ),
    class = "coreact_data"
  )
}

#' Filter features by prevalence
#'
#' Removes features that do not meet the specified prevalence threshold.
#'
#' @param data A 'coreact_data' object.
#' @param threshold Numeric. Values less than 1 are interpreted as a fraction (percentage). Values greater than or equal to 1 are interpreted as absolute counts.
#'
#' @return A filtered 'coreact_data' object, or NULL if all features are removed.
#' @export
filter_prevalence <- function(data, threshold) {
  stopifnot(inherits(data, "coreact_data"))
  stopifnot(is.numeric(threshold) && length(threshold) == 1)

  cut_off <- if (threshold < 1 && threshold > 0) {
    ceiling(threshold * data$n_samples)
  } else {
    threshold
  }

  keep_idx <- which(data$prevalence >= cut_off)

  if (length(keep_idx) == 0) {
    warning(sprintf("[%s] Prevalence filter removed all features.", data$name))
    return(NULL)
  }

  message(sprintf("[%s] Filtering: Keeping %d / %d features (Threshold >= %s)",
                  data$name, length(keep_idx), data$n_features, cut_off))

  new_coreact_data(
    mat = data$mat[keep_idx, , drop = FALSE],
    meta = data$meta[keep_idx, , drop = FALSE],
    name = data$name
  )
}
