# Shared helpers for the test suite. testthat sources helper-*.R before tests.

# Build a coreact_data object from a random binary matrix.
# Uses rand.x = rep(1, n) so every stored value is exactly 1.
make_obj <- function(n_rows, n_cols, prefix = "G", density = 0.5, name = NULL) {
  mat <- Matrix::rsparsematrix(n_rows, n_cols, density = density,
                               rand.x = function(n) rep(1, n))

  row_ids <- paste0(prefix, seq_len(n_rows))
  col_ids <- paste0("S", seq_len(n_cols))

  rownames(mat) <- row_ids
  colnames(mat) <- col_ids

  meta <- data.frame(
    id = row_ids,
    desc = paste("Desc", row_ids),
    row.names = row_ids,
    stringsAsFactors = FALSE
  )

  new_coreact_data(mat, meta, name = if (is.null(name)) paste0("Obj_", prefix) else name)
}

# Write a temporary TSV with ID + Desc metadata columns followed by binary samples.
# `cols_arg` is either a number of samples or an explicit vector of sample names.
write_pipeline_tsv <- function(rows, cols_arg, prefix, dir = tempdir()) {
  if (is.numeric(cols_arg) && length(cols_arg) == 1) {
    col_names <- paste0("S", seq_len(cols_arg))
  } else {
    col_names <- cols_arg
  }
  n_cols <- length(col_names)

  mat <- matrix(sample(0:1, rows * n_cols, replace = TRUE), nrow = rows)
  df <- as.data.frame(mat)
  colnames(df) <- col_names

  df$FeatureID <- paste0(prefix, seq_len(rows))
  df$Desc <- paste0("Desc_", df$FeatureID)
  df <- df[, c("FeatureID", "Desc", col_names)]

  tmp <- tempfile(fileext = ".tsv", tmpdir = dir)
  write.table(df, file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  tmp
}

# Deliberately naive reference implementation: loops over every pair and
# computes each statistic from the raw vectors. Used as an oracle for the
# vectorised engine. Returns a data.frame sorted by (idx_x, idx_y).
brute_force_coreact <- function(mat_x, mat_y) {
  mx <- as.matrix(mat_x)
  my <- as.matrix(mat_y)
  n <- ncol(mx)
  out <- list()
  for (i in seq_len(nrow(mx))) {
    for (j in seq_len(nrow(my))) {
      a <- mx[i, ] == 1
      b <- my[j, ] == 1
      int <- sum(a & b)
      sx <- sum(a)
      sy <- sum(b)
      uni <- sum(a | b)
      out[[length(out) + 1]] <- data.frame(
        idx_x = i, idx_y = j,
        intersection = int, size_x = sx, size_y = sy, union = uni,
        jaccard = if (uni == 0) NA_real_ else int / uni,
        overlap = if (min(sx, sy) == 0) NA_real_ else int / min(sx, sy),
        p_val = stats::phyper(int - 1, sx, n - sx, sy, lower.tail = FALSE)
      )
    }
  }
  res <- do.call(rbind, out)
  res[order(res$idx_x, res$idx_y), ]
}

# Sort a results tibble into canonical order for comparisons.
sort_results <- function(res) {
  res <- res[order(res$feature_x, res$feature_y), ]
  rownames(res) <- NULL
  res
}
