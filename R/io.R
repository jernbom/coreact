#' Import TSV as Coreact Data Object
#'
#' Reads a TSV, separates metadata columns from sample data, and converts
#' the sample data to a sparse matrix.
#'
#' @param path String. Path to the .tsv file.
#' @param meta_cols Vector (integer or character). Indices or names of columns
#'   containing metadata.
#' @param name String. Identifier for the dataset.
#' @param n_cores Integer. Number of threads for fast reading.
#'
#' @return A 'coreact_data' object.
#' @importFrom data.table fread setDF
#' @export
import_coreact_tsv <- function(path, meta_cols = 1, name = "unknown", n_cores = 1) {

  if (!file.exists(path)) stop(sprintf("File not found: %s", path))

  message(sprintf("[%s] Reading file: %s using %d threads...", name, basename(path), n_cores))

  # 1. Fast Read (Keep the speed)
  dt <- data.table::fread(path, header = TRUE, nThread = n_cores)

  # 2. Stabilize (Convert to standard data.frame immediately)
  # This avoids all data.table scoping/syntax issues while keeping the data in memory.
  data.table::setDF(dt)

  # 3. Identify Metadata vs Data Columns
  all_cols <- names(dt)

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

  # 4. Separate using Standard Base R
  # `drop = FALSE` guarantees we get a data.frame, never a vector.
  meta_df <- dt[, meta_idx, drop = FALSE]

  # Extract Data (Remove metadata columns)
  mat_raw <- as.matrix(dt[, -meta_idx, drop = FALSE])

  # Cleanup
  rm(dt); gc()

  # 5. Create Sparse Matrix
  mat_sparse <- Matrix::Matrix(mat_raw, sparse = TRUE)
  rm(mat_raw); gc()

  new_coreact_data(mat_sparse, meta_df, name = name)
}
