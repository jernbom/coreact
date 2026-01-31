library(testthat)
library(Matrix)
library(dplyr)
library(data.table)

# --- Helpers ---

# Create dummy Coreact Object
# Uses rsparsematrix with a specific rand.x to ensure all non-zero values are exactly 1
# This is required to pass the strict 'is_binary_matrix' check in new_coreact_data
make_obj <- function(n_rows, n_cols, prefix = "G") {
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
# Accepts explicit column names to allow testing of mismatched files
write_pipeline_tsv <- function(rows, cols_arg, prefix) {

  if (is.numeric(cols_arg) && length(cols_arg) == 1) {
    col_names <- paste0("S", 1:cols_arg)
  } else {
    col_names <- cols_arg
  }
  n_cols <- length(col_names)

  # sample(0:1) generates binary integers
  mat <- matrix(sample(0:1, rows * n_cols, replace = TRUE), nrow = rows)
  df <- as.data.frame(mat)
  colnames(df) <- col_names

  # Add ID column
  df$FeatureID <- paste0(prefix, 1:rows)
  df$Desc <- paste0("Desc_", df$FeatureID)

  # Reorder: ID, Desc, Samples...
  df <- df[, c("FeatureID", "Desc", col_names)]

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
      names = c("DataX", "DataY"),
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

  # Check Sidecar Metadata
  side_x <- sub("\\.tsv$", "_metadata_x.tsv", out_path)
  expect_true(file.exists(side_x))
})

test_that("coreact_pipeline supports sample_cols filtering", {
  # Scenario: Two files have different samples, but share S3 and S4.
  # X: S1, S2, S3, S4
  # Y: S3, S4, S5, S6
  path_x <- write_pipeline_tsv(10, c("S1", "S2", "S3", "S4"), "X")
  path_y <- write_pipeline_tsv(10, c("S3", "S4", "S5", "S6"), "Y")
  out_path <- file.path(tempdir(), "filtered_results.tsv")

  # 1. Without sample_cols -> Should FAIL strict consistency check
  expect_error(
    coreact_pipeline(
      paths = c(path_x, path_y),
      out_path = out_path,
      meta_cols = c("FeatureID", "Desc"),
      feature_ids = 1
    ),
    "Sample identifiers.*not identical"
  )

  # 2. With sample_cols -> Should PASS
  # Note: fdr_threshold=1.0 ensures we don't filter out results, validating the pipeline ran
  expect_message(
    coreact_pipeline(
      paths = c(path_x, path_y),
      out_path = out_path,
      meta_cols = c("FeatureID", "Desc"),
      feature_ids = 1,
      sample_cols = c("S3", "S4"),
      fdr_threshold = 1.0
    ),
    "Pipeline completed successfully"
  )

  # Check that output is generated
  expect_true(file.exists(out_path))
})

test_that("coreact_pipeline validates mismatched inputs", {
  p1 <- write_pipeline_tsv(5, 3, "A")
  p2 <- write_pipeline_tsv(5, 4, "B") # Different n columns

  # We must specify BOTH metadata columns (1 and 2).
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
