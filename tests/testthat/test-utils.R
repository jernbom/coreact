library(testthat)
library(Matrix)
library(data.table)

# --- Helper to create a dummy S3 object for testing ---
# We recreate a minimal valid object here to ensure 'utils' tests
# are relatively isolated, though they rely on the S3 subset method `[.coreact_data`
# being present in the package namespace.
make_dummy_obj <- function() {
  # 5 Features, 4 Samples
  # Counts:
  # F1: 1, 1, 1, 1 (Sum 4)
  # F2: 1, 1, 1, 0 (Sum 3)
  # F3: 1, 1, 0, 0 (Sum 2)
  # F4: 1, 0, 0, 0 (Sum 1)
  # F5: 0, 0, 0, 0 (Sum 0)

  vals <- c(rep(1, 4), rep(1, 3), rep(1, 2), 1)
  # i (row indices, 0-based in C, 1-based in R construction)
  i <- c(rep(0, 4), rep(1, 3), rep(2, 2), 3)
  # j (col indices)
  j <- c(0:3, 0:2, 0:1, 0)

  mat <- Matrix::sparseMatrix(i = i, j = j, x = vals, dims = c(5, 4), index1 = FALSE)
  rownames(mat) <- paste0("F", 1:5)
  colnames(mat) <- paste0("S", 1:4)

  meta <- data.frame(
    id = paste0("F", 1:5),
    symbol = paste0("Gene", 1:5),
    row.names = paste0("F", 1:5),
    stringsAsFactors = FALSE
  )

  structure(list(mat = mat, meta = meta, name = "TestObj"), class = "coreact_data")
}


test_that("resolve_xy_args handles input variations correctly", {

  # 1. Single Atomic Value (Recycle)
  res1 <- resolve_xy_args(10, "arg")
  expect_equal(res1$x, 10)
  expect_equal(res1$y, 10)

  # 2. List of length 1 (Recycle)
  res2 <- resolve_xy_args(list("A"), "arg")
  expect_equal(res2$x, "A")
  expect_equal(res2$y, "A")

  # 3. List of length 2 (Specific)
  res3 <- resolve_xy_args(list(1, 2), "arg")
  expect_equal(res3$x, 1)
  expect_equal(res3$y, 2)

  # 4. Error Conditions
  expect_error(resolve_xy_args(list(1, 2, 3), "test_arg"), "Argument 'test_arg'.*length 1 or 2")
  expect_error(resolve_xy_args(list(), "test_arg"), "Argument 'test_arg'.*length 1 or 2")
})


test_that("filter_by_prevalence filters correctly", {
  obj <- make_dummy_obj()

  # 1. Zero threshold (Keep all)
  res0 <- filter_by_prevalence(obj, 0)
  expect_equal(nrow(res0$mat), 5)

  # 2. Percentage Threshold (0.5 = 50% of 4 samples = 2 samples)
  # Should keep F1(4), F2(3), F3(2). Drop F4(1), F5(0).
  res_pct <- filter_by_prevalence(obj, 0.5)
  expect_equal(nrow(res_pct$mat), 3)
  expect_equal(rownames(res_pct$mat), c("F1", "F2", "F3"))

  # 3. Absolute Threshold (3 samples)
  # Should keep F1(4), F2(3). Drop rest.
  res_abs <- filter_by_prevalence(obj, 3)
  expect_equal(nrow(res_abs$mat), 2)
  expect_equal(rownames(res_abs$mat), c("F1", "F2"))

  # 4. Error: Filter removes everything
  # Threshold 5 (Max prevalence is 4)
  expect_error(filter_by_prevalence(obj, 5), "Filter removed all features")
})


test_that("write_metadata_sidecar writes correct files", {
  obj <- make_dummy_obj()

  # Create a temporary directory for file output
  tmp_dir <- tempdir()

  # Case 1: Standard file with extension
  out_path <- file.path(tmp_dir, "results.tsv")
  keep_ids <- c("F1", "F3") # Keep a subset

  # Run function
  expect_message(
    write_metadata_sidecar(obj, keep_ids, out_path, suffix = "x"),
    "Writing x metadata"
  )

  # Verify File Creation
  expected_file <- file.path(tmp_dir, "results_metadata_x.tsv")
  expect_true(file.exists(expected_file))

  # Verify Content
  df <- data.table::fread(expected_file)
  expect_equal(nrow(df), 2)
  expect_equal(df$feature_id, c("F1", "F3"))
  expect_equal(df$symbol, c("Gene1", "Gene3")) # Check metadata column came along

  # Case 2: File without extension
  out_path_noext <- file.path(tmp_dir, "results_base")
  write_metadata_sidecar(obj, keep_ids, out_path_noext, suffix = "y")

  expected_file_noext <- file.path(tmp_dir, "results_base_metadata_y.tsv")
  expect_true(file.exists(expected_file_noext))

  # Case 3: No metadata columns in object (should still write feature_id)
  obj_nometa <- obj
  obj_nometa$meta <- obj_nometa$meta[, character(0)] # empty cols

  write_metadata_sidecar(obj_nometa, keep_ids, out_path, suffix = "z")
  expected_file_z <- file.path(tmp_dir, "results_metadata_z.tsv")

  df_z <- data.table::fread(expected_file_z)
  expect_equal(colnames(df_z), "feature_id")
  expect_equal(nrow(df_z), 2)
})
