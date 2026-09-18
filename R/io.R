#' Import TSV as Coreact Data Object
#'
#' Reads a TSV, separates metadata columns from sample data, and converts
#' the sample data to a sparse matrix.
#'
#' @param path String. Path to the .tsv file.
#' @param meta_cols Vector (integer or character). Indices or names of columns
#'    containing metadata.
#' @param sample_cols Vector (integer or character) or NULL. Indices or names of
#'    columns to be used as sample data.
#'    If `NULL`, all columns NOT in `meta_cols` are used.
#' @param feature_id Vector (integer or character). The columns to use to construct
#'    unique row names (Feature IDs). Must be a subset of `meta_cols`.
#'    If multiple columns are provided, they are concatenated using `feature_id_sep`.
#'    If NULL, row numbers are used.
#' @param feature_id_sep String. Separator used when concatenating multiple
#'    feature_id columns. Default is "|".
#' @param name String. Identifier for the dataset.
#' @param n_cores Integer. Number of threads for fast reading.
#' @param binarize Logical. Should the data matrix be converted to binary?
#'   If `TRUE`, every value `>= 1` becomes 1 and everything else 0. If `FALSE`
#'   (default), the sample columns must already contain only 0 and 1.
#'
#' @details
#' Sample columns must be numeric and complete: text columns and missing
#' values are rejected with an error naming the offending columns (use
#' [preview_cols()] to inspect headers and choose `meta_cols` /
#' `sample_cols`). The resulting matrix is always a general sparse
#' `dgCMatrix`.
#'
#' @return A 'coreact_data' object.
#' @importFrom data.table fread setDF
#' @importFrom utils head
#' @export
#' @examples
#' path_x <- system.file("extdata", "input_x.tsv", package = "coreact")
#' x <- import_coreact_tsv(path_x, meta_cols = 1:4, feature_id = "u_pep_id", name = "X")
#' x
#' x$mat[1:3, 1:5]
import_coreact_tsv <- function(path,
                               meta_cols = 1,
                               sample_cols = NULL,
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
  n_total_cols <- ncol(dt)

  # 2. Resolve Metadata Indices
  if (is.character(meta_cols)) {
    missing <- setdiff(meta_cols, all_cols)
    if (length(missing) > 0) stop(sprintf("[%s] Metadata columns not found: %s", name, paste(missing, collapse=", ")))
    meta_idx <- match(meta_cols, all_cols)
  } else if (is.numeric(meta_cols)) {
    if (any(meta_cols > n_total_cols) || any(meta_cols < 1)) stop(sprintf("[%s] Metadata indices out of bounds.", name))
    meta_idx <- meta_cols
  } else {
    stop("meta_cols must be a vector of column names or indices.")
  }

  # 3. Construct Row Names (Feature IDs)
  if (!is.null(feature_id)) {
    # Resolve feature_id to indices
    if (is.character(feature_id)) {
      missing_id <- setdiff(feature_id, all_cols)
      if (length(missing_id) > 0) stop(sprintf("[%s] Feature ID columns not found: %s", name, paste(missing_id, collapse=", ")))
      feat_idx <- match(feature_id, all_cols)
    } else {
      if (any(feature_id > n_total_cols) || any(feature_id < 1)) stop(sprintf("[%s] Feature ID indices out of bounds.", name))
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
  if (is.null(sample_cols)) {
    # Default: All columns except metadata
    mat_cols_idx <- setdiff(seq_len(n_total_cols), meta_idx)
  } else {
    # Specific Columns Requested
    if (is.character(sample_cols)) {
      missing_samp <- setdiff(sample_cols, all_cols)
      if (length(missing_samp) > 0) {
        stop(sprintf("[%s] Sample columns not found: %s", name, paste(head(missing_samp, 5), collapse=", ")))
      }
      mat_cols_idx <- match(sample_cols, all_cols)
    } else if (is.numeric(sample_cols)) {
      # Interpret indices relative to original file
      if (any(sample_cols > n_total_cols) || any(sample_cols < 1)) {
        stop(sprintf("[%s] Sample column indices out of bounds relative to original TSV.", name))
      }
      mat_cols_idx <- sample_cols
    } else {
      stop("sample_cols must be a vector of column names or indices.")
    }
  }

  if (length(mat_cols_idx) == 0) stop(sprintf("[%s] No sample columns selected for matrix.", name))

  if (anyDuplicated(mat_cols_idx)) {
    stop(sprintf("[%s] Duplicate sample columns requested: %s", name,
                 paste(unique(all_cols[mat_cols_idx[duplicated(mat_cols_idx)]]), collapse = ", ")))
  }
  overlap_cols <- intersect(mat_cols_idx, meta_idx)
  if (length(overlap_cols) > 0) {
    stop(sprintf("[%s] Columns cannot be both metadata and sample columns: %s", name,
                 paste(all_cols[overlap_cols], collapse = ", ")))
  }

  # Sample columns must be numeric; fread types a column as character/logical
  # when it contains text, or all-empty fields, respectively.
  col_is_num <- vapply(dt[mat_cols_idx], is.numeric, logical(1))
  if (!all(col_is_num)) {
    bad <- all_cols[mat_cols_idx][!col_is_num]
    stop(sprintf(
      "[%s] Sample columns must be numeric (0/1). Non-numeric columns: %s%s",
      name, paste(head(bad, 5), collapse = ", "),
      if (length(bad) > 5) sprintf(", ... (+%d more)", length(bad) - 5) else ""
    ))
  }

  # Create Matrix
  mat_raw <- as.matrix(dt[, mat_cols_idx, drop = FALSE])
  rownames(mat_raw) <- row_names

  # Cleanup
  rm(dt)

  n_na <- sum(is.na(mat_raw))
  if (n_na > 0) {
    stop(sprintf("[%s] Sample data contains %d missing (NA) values; coreact requires complete 0/1 data.",
                 name, n_na))
  }

  if (isTRUE(binarize)) {
    mat_raw <- +(mat_raw >= 1)
  } else if (!is_binary_matrix(mat_raw)) {
    stop(sprintf(paste0("[%s] Sample data must be a binary matrix (0/1), but other values were found. ",
                        "Use binarize = TRUE to threshold counts at >= 1."), name))
  }

  # 6. Create Sparse Matrix
  mat_sparse <- as_binary_sparse(mat_raw)
  rm(mat_raw)

  new_coreact_data(mat_sparse, meta_df, name = name)
}

#' Convert a Dense 0/1 Matrix to dgCMatrix (Internal)
#'
#' `Matrix::Matrix(x, sparse = TRUE)` picks a class based on the *shape* of
#' the data (diagonal, symmetric or triangular square inputs come back as
#' `ddiMatrix`, `dsCMatrix` or `dtCMatrix`). Building the sparse matrix from
#' the non-zero positions instead guarantees a general `dgCMatrix`, which is
#' what the rest of the package assumes, and never materialises anything
#' larger than the index of non-zero entries.
#'
#' @param x A base numeric matrix containing only 0 and 1, with dimnames.
#' @return A `dgCMatrix` with the same dimensions and dimnames as `x`.
#' @keywords internal
as_binary_sparse <- function(x) {
  nz <- which(x != 0)
  Matrix::sparseMatrix(
    i = ((nz - 1L) %% nrow(x)) + 1L,
    j = ((nz - 1L) %/% nrow(x)) + 1L,
    x = rep(1, length(nz)),
    dims = dim(x),
    dimnames = dimnames(x)
  )
}
