library(testthat)
library(Matrix)

# --- Helper to create temp TSV files ---
write_temp_tsv <- function(df) {
  tmp <- tempfile(fileext = ".tsv")
  write.table(df, file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  return(tmp)
}

# --- Helper Data Frame ---
# 3 Genes, 2 Metadata cols (ID, Symbol), 2 Samples (S1, S2)
# Contains counts (10, 5, 2) which are NOT valid unless binarized.
create_df <- function() {
  data.frame(
    GeneID = c("G1", "G2", "G3"),
    Symbol = c("Alpha", "Beta", "Gamma"),
    S1 = c(10, 0, 5),
    S2 = c(0, 2, 5),
    stringsAsFactors = FALSE
  )
}

test_that("import_coreact_tsv imports valid data correctly", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # 1. Standard Import with Binarization
  obj <- import_coreact_tsv(
    path,
    meta_cols = c("GeneID", "Symbol"),
    feature_id = "GeneID",
    name = "TestImport",
    binarize = TRUE # Required because data contains 10, 5, etc.
  )

  expect_s3_class(obj, "coreact_data")
  expect_equal(obj$name, "TestImport")

  # Check Dimensions: 3 Features x 2 Samples
  expect_equal(dim(obj$mat), c(3, 2))
  expect_equal(dim(obj$meta), c(3, 2))

  # Check Content
  expect_equal(rownames(obj$mat), c("G1", "G2", "G3"))
  expect_equal(colnames(obj$mat), c("S1", "S2"))
  expect_equal(obj$meta$Symbol, c("Alpha", "Beta", "Gamma"))

  # Check Values are now Binary
  # Original S1 was 10, 0, 5 -> Should be 1, 0, 1
  expect_equal(as.vector(obj$mat[, "S1"]), c(1, 0, 1))

  # Check Matrix Type (Sparse)
  expect_true(methods::is(obj$mat, "sparseMatrix"))
})

test_that("import_coreact_tsv handles different column specifications", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # 1. Numeric Indices for Meta
  obj_idx <- import_coreact_tsv(path, meta_cols = 1:2, feature_id = 1, binarize = TRUE)
  expect_equal(colnames(obj_idx$meta), c("GeneID", "Symbol"))
  expect_equal(rownames(obj_idx$mat), c("G1", "G2", "G3"))

  # 2. Composite Feature ID (Concatenation)
  obj_comp <- import_coreact_tsv(
    path,
    meta_cols = 1:2,
    feature_id = c("GeneID", "Symbol"),
    feature_id_sep = "_",
    binarize = TRUE
  )

  expected_ids <- c("G1_Alpha", "G2_Beta", "G3_Gamma")
  expect_equal(rownames(obj_comp$mat), expected_ids)
  expect_equal(rownames(obj_comp$meta), expected_ids)
})

test_that("import_coreact_tsv handles implicit IDs", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # If feature_id is NULL, it should use row numbers
  obj <- import_coreact_tsv(path, meta_cols = 1:2, feature_id = NULL, binarize = TRUE)

  expect_equal(rownames(obj$mat), c("1", "2", "3"))
  expect_equal(rownames(obj$meta), c("1", "2", "3"))
})

test_that("import_coreact_tsv enforce binary requirement", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # Fail Case: Importing counts without binarize = TRUE
  # The strict constructor (new_coreact_data) should reject the matrix.
  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, feature_id = 1, binarize = FALSE),
    "must be a binary matrix"
  )
})

test_that("import_coreact_tsv validates inputs", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # 1. File Not Found
  expect_error(import_coreact_tsv("nonexistent_file.tsv"), "File not found")

  # 2. Missing Metadata Columns
  expect_error(import_coreact_tsv(path, meta_cols = "WrongCol"), "Metadata columns not found")

  # 3. Missing Feature ID Columns
  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, feature_id = "WrongID"),
    "Feature ID columns not found"
  )

  # 4. Feature ID not inside Meta Cols
  expect_error(
    import_coreact_tsv(path, meta_cols = "GeneID", feature_id = "S1"),
    "must be present in 'meta_cols'"
  )

  # 5. Out of bounds indices
  expect_error(import_coreact_tsv(path, meta_cols = 99), "indices out of bounds")
})

test_that("import_coreact_tsv validates uniqueness", {
  # Create data with duplicate IDs
  df_dup <- create_df()
  df_dup <- rbind(df_dup, df_dup[1, ]) # Duplicate first row
  path <- write_temp_tsv(df_dup)

  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, feature_id = "GeneID"),
    "Constructed Feature IDs are not unique"
  )
})
