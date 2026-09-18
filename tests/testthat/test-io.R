library(testthat)
library(Matrix)

# --- Helper to create temp TSV files ---
write_temp_tsv <- function(df) {
  tmp <- tempfile(fileext = ".tsv")
  write.table(df, file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  return(tmp)
}

# --- Helper Data Frame ---
# 3 Genes, 2 Metadata cols (ID, Symbol), 3 Samples (S1, S2, S3)
# Contains counts (10, 5, 2) which are NOT valid unless binarized.
create_df <- function() {
  data.frame(
    GeneID = c("G1", "G2", "G3"),
    Symbol = c("Alpha", "Beta", "Gamma"),
    S1 = c(10, 0, 5),
    S2 = c(0, 2, 5),
    S3 = c(1, 1, 0),
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

  # Check Dimensions: 3 Features x 3 Samples (Default imports all non-meta)
  expect_equal(dim(obj$mat), c(3, 3))
  expect_equal(dim(obj$meta), c(3, 2))

  # Check Content
  expect_equal(rownames(obj$mat), c("G1", "G2", "G3"))
  expect_equal(colnames(obj$mat), c("S1", "S2", "S3"))
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

test_that("import_coreact_tsv handles sample_cols filtering", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # Columns: 1=GeneID, 2=Symbol, 3=S1, 4=S2, 5=S3

  # 1. Filter by Name (Subset to S1 and S3 only)
  obj_name <- import_coreact_tsv(
    path,
    meta_cols = c("GeneID", "Symbol"),
    sample_cols = c("S1", "S3"),
    feature_id = "GeneID",
    binarize = TRUE
  )
  expect_equal(colnames(obj_name$mat), c("S1", "S3"))
  expect_equal(dim(obj_name$mat), c(3, 2)) # 3 rows, 2 columns

  # 2. Filter by Index (Original Indices: S2 is 4, S3 is 5)
  # Indices must refer to the ORIGINAL file structure
  obj_idx <- import_coreact_tsv(
    path,
    meta_cols = 1:2,
    sample_cols = c(4, 5),
    feature_id = 1,
    binarize = TRUE
  )
  expect_equal(colnames(obj_idx$mat), c("S2", "S3"))

  # 3. Reordering via Indices
  # Requesting S3 (5) then S1 (3)
  obj_reorder <- import_coreact_tsv(
    path,
    meta_cols = 1:2,
    sample_cols = c(5, 3),
    feature_id = 1,
    binarize = TRUE
  )
  expect_equal(colnames(obj_reorder$mat), c("S3", "S1"))
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

  # 5. Out of bounds metadata indices
  expect_error(import_coreact_tsv(path, meta_cols = 99), "indices out of bounds")
})

test_that("import_coreact_tsv validates sample_cols", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # 1. Missing sample column name
  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, sample_cols = c("S1", "BadCol")),
    "Sample columns not found"
  )

  # 2. Out of bounds sample indices
  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, sample_cols = c(3, 99)),
    "Sample column indices out of bounds"
  )

  # 3. No sample columns selected (e.g. empty match or explicit empty)
  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, sample_cols = character(0)),
    "No sample columns selected"
  )
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

test_that("import_coreact_tsv rejects non-numeric sample columns with a clear message", {
  df <- create_df()
  path <- write_temp_tsv(df)

  # Only GeneID is metadata, so the text column Symbol falls into the samples
  expect_error(
    import_coreact_tsv(path, meta_cols = "GeneID", feature_id = "GeneID", binarize = TRUE),
    "Non-numeric columns: Symbol"
  )
})

test_that("import_coreact_tsv rejects missing values", {
  df <- create_df()
  df$S2[2] <- NA
  path <- write_temp_tsv(df)

  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, feature_id = 1, binarize = TRUE),
    "1 missing \\(NA\\) values"
  )
})

test_that("import_coreact_tsv rejects duplicate and overlapping column selections", {
  df <- create_df()
  path <- write_temp_tsv(df)

  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, sample_cols = c("S1", "S1", "S2"), binarize = TRUE),
    "Duplicate sample columns requested: S1"
  )
  expect_error(
    import_coreact_tsv(path, meta_cols = 1:2, sample_cols = c("Symbol", "S1"), binarize = TRUE),
    "both metadata and sample columns: Symbol"
  )
})

test_that("import_coreact_tsv always yields a dgCMatrix", {
  # A square, symmetric 0/1 block: Matrix::Matrix() would return dsCMatrix
  df <- data.frame(
    id = c("a", "b", "c"),
    S1 = c(1, 1, 0), S2 = c(1, 1, 1), S3 = c(0, 1, 1),
    stringsAsFactors = FALSE
  )
  path <- write_temp_tsv(df)
  obj <- import_coreact_tsv(path, meta_cols = 1, feature_id = 1)
  expect_s4_class(obj$mat, "dgCMatrix")
  expect_equal(as.matrix(obj$mat), matrix(c(1, 1, 0, 1, 1, 1, 0, 1, 1), 3,
                                          dimnames = list(c("a", "b", "c"), c("S1", "S2", "S3"))))

  # Diagonal block: would be ddiMatrix
  df2 <- data.frame(id = c("a", "b"), S1 = c(1, 0), S2 = c(0, 1), stringsAsFactors = FALSE)
  obj2 <- import_coreact_tsv(write_temp_tsv(df2), meta_cols = 1, feature_id = 1)
  expect_s4_class(obj2$mat, "dgCMatrix")

  # All-zero data is still a valid (empty) dgCMatrix
  df3 <- data.frame(id = c("a", "b"), S1 = c(0, 0), S2 = c(0, 0), stringsAsFactors = FALSE)
  obj3 <- import_coreact_tsv(write_temp_tsv(df3), meta_cols = 1, feature_id = 1)
  expect_s4_class(obj3$mat, "dgCMatrix")
  expect_equal(sum(obj3$mat), 0)
})

test_that("as_binary_sparse round-trips dense 0/1 matrices", {
  withr::local_seed(8)
  m <- matrix(rbinom(60, 1, 0.3), 6, 10,
              dimnames = list(paste0("f", 1:6), paste0("s", 1:10)))
  sp <- as_binary_sparse(m)
  expect_s4_class(sp, "dgCMatrix")
  expect_equal(as.matrix(sp), m)
  expect_equal(dimnames(sp), dimnames(m))
})

test_that("binarize = TRUE thresholds at >= 1 and drops negatives", {
  df <- create_df()
  df$S1 <- c(-3, 0.5, 2)
  path <- write_temp_tsv(df)
  obj <- import_coreact_tsv(path, meta_cols = 1:2, feature_id = 1, binarize = TRUE)
  expect_equal(as.vector(obj$mat[, "S1"]), c(0, 0, 1))
})
