#' Import TSV as Coreact Data Object
#'
#' Reads a TSV, separates metadata columns from sample data, and converts
#' the sample data to a sparse matrix.
#'
#' @param path String. Path to the .tsv file.
#' @param meta_cols Vector (integer or character). Indices or names of columns
#'   containing metadata.
#' @param feature_id Vector (integer or character). The columns to use to construct
#'   unique row names (Feature IDs). Must be a subset of `meta_cols`.
#'   If multiple columns are provided, they are concatenated using `feature_id_sep`.
#'   If NULL, row numbers are used.
#' @param feature_id_sep String. Separator used when concatenating multiple
#'   feature_id columns. Default is "|".
#' @param name String. Identifier for the dataset.
#' @param n_cores Integer. Number of threads for fast reading.
#' @param binarize Logical. Should the data matrix be converted to binary?
#'
#' @return A 'coreact_data' object.
#' @importFrom data.table fread setDF
#' @export
import_coreact_tsv <- function(path,
                               meta_cols = 1,
                               feature_id = NULL,
                               feature_id_sep = "|",
                               name = "unknown",
                               n_cores = 1,
                               binarize = FALSE) {

  if (!file.exists(path)) stop(sprintf("File not found: %s", path))

  message(sprintf("[%s] Reading file: %s using %d threads...", name, basename(path), n_cores))

  # 1. Fast Read & Stabilize
  dt <- data.table::fread(path, header = TRUE, nThread = n_cores)
  data.table::setDF(dt)

  all_cols <- names(dt)

  # 2. Resolve Metadata Indices
  if (is.character(meta_cols)) {
    missing <- setdiff(meta_cols, all_cols)
    if (length(missing) > 0) stop(sprintf("Metadata columns not found: %s", paste(missing, collapse=", ")))
    meta_idx <- match(meta_cols, all_cols)
  } else if (is.numeric(meta_cols)) {
    if (any(meta_cols > ncol(dt)) || any(meta_cols < 1)) stop("Metadata indices out of bounds.")
    meta_idx <- meta_cols
  } else {
    stop("meta_cols must be a vector of column names or indices.")
  }

  # 3. Construct Row Names (Feature IDs)
  if (!is.null(feature_id)) {

    # Resolve feature_id to indices
    if (is.character(feature_id)) {
      missing_id <- setdiff(feature_id, all_cols)
      if (length(missing_id) > 0) stop(sprintf("Feature ID columns not found: %s", paste(missing_id, collapse=", ")))
      feat_idx <- match(feature_id, all_cols)
    } else {
      if (any(feature_id > ncol(dt)) || any(feature_id < 1)) stop("Feature ID indices out of bounds.")
      feat_idx <- feature_id
    }

    # feature_id must be within meta_cols
    if (!all(feat_idx %in% meta_idx)) {
      stop(sprintf("[%s] Error: 'feature_id' columns must be present in 'meta_cols'.", name))
    }

    # Construct the ID vector
    if (length(feat_idx) == 1) {
      row_names <- as.character(dt[[feat_idx]])
    } else {
      cols_to_paste <- dt[, feat_idx, drop = FALSE]
      row_names <- do.call(paste, c(cols_to_paste, sep = feature_id_sep))
    }

    # Check Uniqueness
    if (anyDuplicated(row_names)) {
      stop(sprintf("[%s] Error: Constructed Feature IDs are not unique.", name))
    }

  } else {
    row_names <- as.character(seq_len(nrow(dt)))
  }

  # 4. Separate Metadata
  meta_df <- dt[, meta_idx, drop = FALSE]
  rownames(meta_df) <- row_names

  # 5. Extract Matrix Data
  # Remove all metadata columns
  mat_raw <- as.matrix(dt[, -meta_idx, drop = FALSE])
  rownames(mat_raw) <- row_names

  if (binarize == TRUE) mat_raw <- +(mat_raw >= 1)

  # Cleanup
  rm(dt); gc()

  # 6. Create Sparse Matrix
  mat_sparse <- Matrix::Matrix(mat_raw, sparse = TRUE)
  rm(mat_raw); gc()

  new_coreact_data(mat_sparse, meta_df, name = name)
}
