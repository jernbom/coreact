library(testthat)
library(Matrix)

# --- Helper to create dummy data for tests ---
create_dummy_data <- function(n_rows = 5, n_cols = 3) {
  mat <- Matrix::rsparsematrix(n_rows, n_cols, density = 0.5,
                               rand.x = function(n) stats::rbinom(n, 1, 1))

  row_ids <- paste0("Gene_", seq_len(n_rows))
  col_ids <- paste0("Sample_", seq_len(n_cols))
  rownames(mat) <- row_ids
  colnames(mat) <- col_ids

  meta <- data.frame(
    id = row_ids,
    desc = paste("Description", seq_len(n_rows)),
    row.names = row_ids,
    stringsAsFactors = FALSE
  )

  return(list(mat = mat, meta = meta))
}

test_that("new_coreact_data creates valid objects", {
  d <- create_dummy_data()
  obj <- new_coreact_data(d$mat, d$meta, name = "TestSet")

  expect_s3_class(obj, "coreact_data")
  expect_equal(obj$name, "TestSet")
  expect_identical(rownames(obj$mat), rownames(obj$meta))
})

test_that("new_coreact_data enforces binary data", {
  d <- create_dummy_data()

  # Inject a non-binary value (e.g., 5)
  bad_mat <- d$mat
  bad_mat[1, 1] <- 5

  expect_error(new_coreact_data(bad_mat, d$meta), "must be a binary matrix")

  # Inject a negative value
  bad_mat2 <- d$mat
  bad_mat2[1, 1] <- -1
  expect_error(new_coreact_data(bad_mat2, d$meta), "must be a binary matrix")
})

test_that("new_coreact_data enforces strict row name presence", {
  d <- create_dummy_data()

  # Case 1: Matrix missing names
  mat_unnamed <- d$mat
  rownames(mat_unnamed) <- NULL
  expect_error(new_coreact_data(mat_unnamed, d$meta), "Input 'mat' must have row names")

  # Case 2: Metadata missing names
  meta_unnamed <- d$meta
  rownames(meta_unnamed) <- NULL
  expect_error(new_coreact_data(d$mat, meta_unnamed), "must be identical")
})

test_that("new_coreact_data enforces identical row names", {
  d <- create_dummy_data()

  # Case 1: Totally different names
  bad_meta <- d$meta
  rownames(bad_meta) <- paste0("Wrong_", 1:nrow(bad_meta))
  expect_error(new_coreact_data(d$mat, bad_meta), "must be identical")

  # Case 2: Same names, different order
  rev_meta <- d$meta[nrow(d$meta):1, ]
  expect_error(new_coreact_data(d$mat, rev_meta), "must be identical")
})

test_that("subsetting works correctly", {
  d <- create_dummy_data(n_rows = 10, n_cols = 5)
  obj <- new_coreact_data(d$mat, d$meta)

  sub_obj <- obj[1:3, ]
  expect_equal(nrow(sub_obj$mat), 3)
  expect_match(sub_obj$name, "_subset")
  # Ensure constructor validated the result
  expect_identical(rownames(sub_obj$mat), rownames(sub_obj$meta))
})

test_that("cbind works correctly", {
  d <- create_dummy_data(n_rows = 5, n_cols = 4)
  obj <- new_coreact_data(d$mat, d$meta, name = "ObjA")

  obj1 <- obj[, 1:2]
  obj2 <- obj[, 3:4]

  res <- cbind.coreact_data(obj1, obj2)

  expect_equal(ncol(res$mat), 4)
  expect_equal(res$name, "ObjA_subset+ObjA_subset")

  # Failure: Mismatched IDs
  d_diff <- create_dummy_data(n_rows = 5, n_cols = 2)
  rownames(d_diff$mat) <- paste0("Other_", 1:5)
  rownames(d_diff$meta) <- paste0("Other_", 1:5)
  obj_diff <- new_coreact_data(d_diff$mat, d_diff$meta)

  expect_error(cbind.coreact_data(obj1, obj_diff), "Feature IDs")
})
