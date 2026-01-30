library(testthat)
library(Matrix)
library(dplyr)
library(data.table)

# --- Helpers ---

# Create dummy Coreact Object
make_obj <- function(n_rows, n_cols, prefix = "G") {
  # rand.x = function(n) rep(1, n) forces all non-zero values to be exactly 1
  mat <- Matrix::rsparsematrix(n_rows, n_cols, density = 0.5,
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

  new_coreact_data(mat, meta, name = paste0("Obj_", prefix))
}

# Write temp TSV for pipeline test
write_pipeline_tsv <- function(rows, cols, prefix) {
  # sample(0:1) generates binary integers
  mat <- matrix(sample(0:1, rows * cols, replace = TRUE), nrow = rows)
  df <- as.data.frame(mat)
  colnames(df) <- paste0("S", 1:cols)

  # Add ID column
  df$FeatureID <- paste0(prefix, 1:rows)
  df$Desc <- paste0("Desc_", df$FeatureID)

  # Reorder: ID, Desc, Samples...
  df <- df[, c("FeatureID", "Desc", paste0("S", 1:cols))]

  tmp <- tempfile(fileext = ".tsv")
  write.table(df, file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  return(tmp)
}

# --- Tests for run_coreact (The Engine Wrapper) ---

test_that("run_coreact computes results correctly", {
  # Manually create binary sparse matrix
  m_x <- Matrix::Matrix(c(1,0, 0,1), 2, 2, sparse=TRUE)
  rownames(m_x) <- c("X1", "X2"); colnames(m_x) <- c("S1", "S2")
  o_x <- new_coreact_data(m_x, data.frame(id=c("X1","X2"), row.names=c("X1","X2")))

  m_y <- Matrix::Matrix(c(1,0, 0,1), 2, 2, sparse=TRUE)
  rownames(m_y) <- c("Y1", "Y2"); colnames(m_y) <- c("S1", "S2")
  o_y <- new_coreact_data(m_y, data.frame(id=c("Y1","Y2"), row.names=c("Y1","Y2")))

  res <- run_coreact(o_x, o_y, filter_config = list(min_intersection = 1))

  expect_equal(nrow(res), 2)
  expect_true("X1" %in% res$feature_x)
  expect_true("Y1" %in% res$feature_y)
  expect_equal(res$jaccard, c(1, 1))
})

test_that("run_coreact handles swap optimization correctly", {
  o_x <- make_obj(2, 5, "X")
  o_y <- make_obj(10, 5, "Y")
  colnames(o_x$mat) <- paste0("S", 1:5)
  colnames(o_y$mat) <- paste0("S", 1:5)

  # Expect message about swapping because nrow(X) < nrow(Y)
  expect_message(
    res <- run_coreact(o_x, o_y, filter_config = list(min_intersection = 0)),
    "Swapping X and Y"
  )

  # Check that results are mapped back correctly (Prefixes should match original input)
  expect_true(all(grepl("^X", res$feature_x)))
  expect_true(all(grepl("^Y", res$feature_y)))
  expect_gt(nrow(res), 0)
})

test_that("run_coreact respects chunk_size", {
  o_x <- make_obj(50, 5, "X")
  o_y <- make_obj(50, 5, "Y")
  colnames(o_x$mat) <- colnames(o_y$mat)

  expect_no_error(
    res <- run_coreact(o_x, o_y, chunk_size = 10, filter_config = list(min_intersection=0))
  )
  expect_gt(nrow(res), 0)
})

# --- Tests for coreact_pipeline (The Full Workflow) ---

test_that("coreact_pipeline runs end-to-end", {
  path_x <- write_pipeline_tsv(10, 4, "GeneX")
  path_y <- write_pipeline_tsv(10, 4, "GeneY")
  out_path <- file.path(tempdir(), "final_results.tsv")

  expect_message(
    coreact_pipeline(
      paths = c(path_x, path_y),
      out_path = out_path,
      meta_cols = c("FeatureID", "Desc"),
      feature_ids = "FeatureID",
      filter_config = list(min_intersection = 1),
      fdr_threshold = 1.0 # Keep everything for test
    ),
    "Pipeline completed successfully"
  )

  expect_true(file.exists(out_path))
  res <- data.table::fread(out_path)
  expect_true(nrow(res) > 0)
  expect_true("p_adj" %in% colnames(res))

  side_x <- sub("\\.tsv$", "_metadata_x.tsv", out_path)
  expect_true(file.exists(side_x))
})

test_that("coreact_pipeline validates inputs", {
  p1 <- write_pipeline_tsv(5, 3, "A")
  p2 <- write_pipeline_tsv(5, 4, "B") # Different n columns

  # We must specify BOTH metadata columns (1 and 2).
  # Otherwise column 2 ("Desc") is treated as data, creating a character matrix.
  expect_error(
    coreact_pipeline(
      paths = c(p1, p2),
      out_path = tempfile(),
      meta_cols = 1:2,
      feature_ids = 1
    ),
    "Sample identifiers.*not identical"
  )
})

test_that("coreact_pipeline handles empty results gracefully", {
  rows=5; cols=4
  df <- as.data.frame(matrix(0, rows, cols))
  colnames(df) <- paste0("S", 1:cols)
  df$ID <- paste0("G", 1:rows)
  df <- df[, c("ID", paste0("S", 1:cols))]

  p1 <- tempfile(fileext=".tsv")
  write.table(df, p1, sep="\t", row.names=FALSE, quote=FALSE)

  out_empty <- file.path(tempdir(), "empty_results.tsv")

  expect_warning(
    coreact_pipeline(
      paths = c(p1, p1),
      out_path = out_empty,
      meta_cols = 1, feature_ids = 1,
      min_prevalence = 0,
      filter_config = list(min_intersection = 1)
    ),
    "No significant results found"
  )

  expect_true(file.exists(out_empty))
  res <- data.table::fread(out_empty)
  expect_equal(nrow(res), 0)
  # Check we have columns even though empty
  expect_true("feature_x" %in% colnames(res))
})
