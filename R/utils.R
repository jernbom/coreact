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

#' Write Sidecar Metadata (Internal)
#'
#' Writes a metadata TSV for features present in the final results.
#' Appends `_metadata_\{suffix\}` to the base filename.
#'
#' @param obj A `coreact_data` object (containing `$meta` and `$mat`).
#' @param keep_ids Character vector. Feature IDs appearing in the final interaction list.
#' @param out_path String. The path of the main results file.
#' @param suffix String. Suffix to distinguish X and Y (e.g., "x" or "y").
#' @return NULL (Writes file to disk).
#' @keywords internal
write_metadata_sidecar <- function(obj, keep_ids, out_path, suffix) {
  # Construct filename: results.tsv -> results_metadata_x.tsv
  # This regex inserts the suffix before the last extension
  if (grepl("\\.", basename(out_path))) {
    # If extension exists, insert suffix before it
    meta_file <- sub("(\\.[^.]+)$", sprintf("_metadata_%s\\1", suffix), out_path)
  } else {
    # If no extension, append suffix + .tsv default
    meta_file <- sprintf("%s_metadata_%s.tsv", out_path, suffix)
  }

  # Extract IDs from the object (assuming rownames of matrix are the master IDs)
  all_ids <- rownames(obj$mat)

  # Check which features are in the final results
  to_keep <- all_ids %in% keep_ids

  # Create the export table: Feature ID first
  export_df <- tibble::tibble(feature_id = all_ids[to_keep])

  # Bind metadata columns if they exist
  if (!is.null(obj$meta) && ncol(obj$meta) > 0) {
    # Note: obj$meta is expected to be row-aligned with obj$mat
    meta_subset <- obj$meta[to_keep, , drop = FALSE]
    export_df <- dplyr::bind_cols(export_df, meta_subset)
  }

  message(sprintf("Writing %s metadata (%d features) to %s ...", suffix, nrow(export_df), meta_file))

  # Ensure directory exists (redundant check but safe)
  out_dir <- dirname(meta_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data.table::fwrite(export_df, file = meta_file, sep = "\t", quote = FALSE)
}

#' Check if a Matrix is Binary (0/1)
#'
#' @description
#' Efficiently determines whether a matrix (dense or sparse) contains only
#' values of 0 and 1.
#'
#' @details
#' This function is optimized for performance depending on the input type:
#' \itemize{
#'   \item **Sparse Matrices (`dgCMatrix`):** It inspects only the stored values
#'   in the `@x` slot. Since implicit (non-stored) values are always 0, verifying
#'   that explicit values are either 0 or 1 is sufficient. This avoids expanding
#'   the matrix.
#'   \item **Dense Matrices:** It uses `any()` with short-circuiting logic. The
#'   scan stops immediately upon finding a value that is not 0 or 1.
#' }
#'
#' @param x A matrix-like object (base `matrix` or `Matrix::sparseMatrix`).
#'
#' @return `TRUE` if all elements are 0 or 1; otherwise `FALSE`.
#'
#' @keywords internal
#' # @export
is_binary_matrix <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    # Check only stored values (@x slot)
    # If a value is stored, it must be 0 or 1.
    # (Implicit values are always 0, which is valid).
    return(all(x@x %in% c(0, 1)))
  } else {
    # Short-circuit check for dense matrices
    # Returns FALSE immediately if a non-0/1 value is found
    return(!any(x != 0 & x != 1))
  }
}
