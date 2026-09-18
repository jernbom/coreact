
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
  rownames(bad_meta) <- paste0("Wrong_", seq_len(nrow(bad_meta)))
  expect_error(new_coreact_data(d$mat, bad_meta), "must be identical")

  # Case 2: Same names, different order
  rev_meta <- d$meta[rev(seq_len(nrow(d$meta))), ]
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

test_that("new_coreact_data accepts a base matrix and converts to dgCMatrix", {
  mat <- matrix(c(1, 0, 1, 1, 0, 0), nrow = 2,
                dimnames = list(c("F1", "F2"), c("S1", "S2", "S3")))
  meta <- data.frame(symbol = c("A", "B"), row.names = c("F1", "F2"))

  obj <- new_coreact_data(mat, meta, name = "base")
  expect_s4_class(obj$mat, "dgCMatrix")
  expect_equal(as.matrix(obj$mat), mat)

  # Logical matrices are fine too
  obj_l <- new_coreact_data(mat == 1, meta, name = "logical")
  expect_equal(as.matrix(obj_l$mat), mat)

  # Non-binary / NA / character base matrices are rejected
  expect_error(new_coreact_data(mat * 2, meta), "must be a binary matrix")
  mat_na <- mat
  mat_na[1, 1] <- NA
  expect_error(new_coreact_data(mat_na, meta), "must not contain NA")
  mat_chr <- matrix("1", 2, 3, dimnames = dimnames(mat))
  expect_error(new_coreact_data(mat_chr, meta), "numeric or logical")

  # Other argument validation
  expect_error(new_coreact_data(mat, meta, name = c("a", "b")), "single string")
  expect_error(new_coreact_data(mat, as.list(meta)), "must be a data.frame")
  expect_error(new_coreact_data(list(), meta), "must be a Matrix")
})

test_that("dim(), nrow() and ncol() work on coreact_data", {
  d <- create_dummy_data(n_rows = 7, n_cols = 4)
  obj <- new_coreact_data(d$mat, d$meta)
  expect_equal(dim(obj), c(7, 4))
  expect_equal(nrow(obj), 7)
  expect_equal(ncol(obj), 4)
})

test_that("subsetting to zero features yields a valid empty object", {
  d <- create_dummy_data(n_rows = 5, n_cols = 3)
  obj <- new_coreact_data(d$mat, d$meta)
  empty <- obj[integer(0), ]
  expect_s3_class(empty, "coreact_data")
  expect_equal(nrow(empty), 0)
  expect_equal(ncol(empty), 3)
  expect_output(print(empty), "0 Features x 3 Samples")
})

test_that("subsetting accepts names, logicals and negative indices", {
  d <- create_dummy_data(n_rows = 5, n_cols = 3)
  obj <- new_coreact_data(d$mat, d$meta)

  by_name <- obj[c("Gene_2", "Gene_4"), ]
  expect_equal(rownames(by_name$mat), c("Gene_2", "Gene_4"))
  expect_equal(by_name$meta$desc, c("Description 2", "Description 4"))

  by_lgl <- obj[c(TRUE, FALSE, TRUE, FALSE, TRUE), ]
  expect_equal(rownames(by_lgl$mat), c("Gene_1", "Gene_3", "Gene_5"))

  by_neg <- obj[-1, "Sample_2"]
  expect_equal(dim(by_neg), c(4, 1))
  expect_equal(colnames(by_neg$mat), "Sample_2")
})

test_that("print method reports dimensions, features and metadata columns", {
  d <- create_dummy_data(n_rows = 5, n_cols = 3)
  obj <- new_coreact_data(d$mat, d$meta, name = "Printed")
  out <- capture.output(ret <- print(obj))
  expect_identical(ret, obj)
  expect_match(out[1], "\\[Printed\\]")
  expect_match(out[2], "5 Features x 3 Samples")
  expect_true(any(grepl("Gene_1, Gene_2, Gene_3 \\.\\.\\.", out)))
  expect_true(any(grepl("Metadata columns: id, desc", out)))
})

test_that("cbind validates classes and refuses duplicated sample names", {
  d <- create_dummy_data(n_rows = 5, n_cols = 4)
  obj <- new_coreact_data(d$mat, d$meta, name = "ObjA")

  expect_error(cbind(obj, d$mat), "must be 'coreact_data' objects")

  # Same object twice -> every sample name duplicated
  expect_error(cbind(obj, obj), "sample names would be duplicated")

  # Disjoint sample names combine; S3 dispatch through cbind() works
  obj2 <- obj
  colnames(obj2$mat) <- paste0("Other_", 1:4)
  combined <- cbind(obj, obj2)
  expect_equal(ncol(combined), 8)
  expect_equal(colnames(combined$mat), c(paste0("Sample_", 1:4), paste0("Other_", 1:4)))
  expect_equal(combined$name, "ObjA+ObjA")

  # Single object is returned as is; nothing returns NULL
  expect_identical(cbind.coreact_data(obj), obj)
  expect_null(cbind.coreact_data())
})
