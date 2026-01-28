#' Import TSV as Coreact Data Object
#'
#' @param path String. Path to the .tsv file.
#' @param meta_cols Vector (names or indices) of metadata columns.
#' @param name String. Identifier for the dataset.
#' @param n_cores Integer. Number of threads for fast reading.
#' @export
import_coreact_tsv <- function(path, meta_cols = 1, name = "unknown", n_cores = 1) {

  if (!file.exists(path)) stop(sprintf("File not found: %s", path))

  message(sprintf("[%s] Reading file: %s using %d threads...", name, basename(path), n_cores))

  # 1. Fast Read
  dt <- data.table::fread(path, header = TRUE, nThread = n_cores)

  # 2. Identify Metadata
  all_cols <- names(dt)
  if (is.character(meta_cols)) {
    missing <- setdiff(meta_cols, all_cols)
    if (length(missing) > 0) stop(sprintf("Metadata columns not found: %s", paste(missing, collapse=", ")))
    meta_idx <- match(meta_cols, all_cols)
  } else {
    if (any(meta_cols > ncol(dt)) || any(meta_cols < 1)) stop("Metadata indices out of bounds.")
    meta_idx <- meta_cols
  }

  # 3. Separate & Convert
  meta_df <- as.data.frame(dt[, ..meta_idx])
  mat_raw <- as.matrix(dt[, -..meta_idx])

  rm(dt); gc()

  # 4. Create Sparse Matrix
  mat_sparse <- Matrix::Matrix(mat_raw, sparse = TRUE)
  rm(mat_raw); gc()

  new_coreact_data(mat_sparse, meta_df, name = name)
}
