# Helpers make_obj(), write_pipeline_tsv(), brute_force_coreact() and
# sort_results() live in helper-data.R.

# --- Tests for run_coreact (The Engine Wrapper) ---

test_that("run_coreact computes results correctly", {
  # Manually create binary sparse matrix
  m_x <- Matrix::Matrix(c(1, 0, 0, 1), 2, 2, sparse = TRUE)
  rownames(m_x) <- c("X1", "X2")
  colnames(m_x) <- c("S1", "S2")
  o_x <- new_coreact_data(m_x, data.frame(id = c("X1", "X2"), row.names = c("X1", "X2")))

  m_y <- Matrix::Matrix(c(1, 0, 0, 1), 2, 2, sparse = TRUE)
  rownames(m_y) <- c("Y1", "Y2")
  colnames(m_y) <- c("S1", "S2")
  o_y <- new_coreact_data(m_y, data.frame(id = c("Y1", "Y2"), row.names = c("Y1", "Y2")))

  res <- run_coreact(o_x, o_y, filter_config = list(min_intersection = 1))

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 2)
  expect_true("X1" %in% res$feature_x)
  expect_true("Y1" %in% res$feature_y)
  expect_equal(res$jaccard, c(1, 1))
  expect_named(res, names(empty_result_tibble()))
})

test_that("run_coreact matches the brute-force reference (all pairs)", {
  withr::local_seed(20240917)
  o_x <- make_obj(30, 12, "X", density = 0.4)
  o_y <- make_obj(8, 12, "Y", density = 0.6)

  res <- run_coreact(o_x, o_y, filter_config = list())
  ref <- brute_force_coreact(o_x$mat, o_y$mat)

  expect_equal(nrow(res), nrow(o_x) * nrow(o_y))

  res <- as.data.frame(res[order(res$idx_x, res$idx_y), ])
  rownames(res) <- NULL
  rownames(ref) <- NULL
  expect_equal(res$feature_x, rownames(o_x$mat)[ref$idx_x])
  expect_equal(res$feature_y, rownames(o_y$mat)[ref$idx_y])
  for (col in c("intersection", "size_x", "size_y", "union", "jaccard", "overlap", "p_val")) {
    expect_equal(res[[col]], ref[[col]], info = col)
  }
})

test_that("run_coreact handles swap optimization correctly", {
  withr::local_seed(42)
  o_x <- make_obj(2, 5, "X")
  o_y <- make_obj(10, 5, "Y")

  # Expect message about swapping because nrow(X) < nrow(Y)
  expect_message(
    res <- run_coreact(o_x, o_y, filter_config = list(min_intersection = 0)),
    "Swapping X and Y"
  )

  # Check that results are mapped back correctly (Prefixes should match original input)
  expect_true(all(grepl("^X", res$feature_x)))
  expect_true(all(grepl("^Y", res$feature_y)))
  expect_equal(nrow(res), 20)
  expect_named(res, names(empty_result_tibble()))

  # size_x must be the prevalence of feature_x in the *original* X, not of the swapped role
  expect_equal(res$size_x, as.numeric(Matrix::rowSums(o_x$mat))[res$idx_x])
  expect_equal(res$size_y, as.numeric(Matrix::rowSums(o_y$mat))[res$idx_y])
  expect_equal(res$feature_x, rownames(o_x$mat)[res$idx_x])
  expect_equal(res$feature_y, rownames(o_y$mat)[res$idx_y])
})

test_that("run_coreact(x, y) and run_coreact(y, x) agree after relabelling", {
  withr::local_seed(7)
  o_a <- make_obj(25, 10, "A", density = 0.4)
  o_b <- make_obj(6, 10, "B", density = 0.5)

  suppressMessages({
    r_ab <- run_coreact(o_a, o_b, filter_config = list(min_intersection = 1))
    r_ba <- run_coreact(o_b, o_a, filter_config = list(min_intersection = 1))
  })

  # Relabel r_ba so that x refers to A and y to B
  r_ba_rel <- r_ba
  names(r_ba_rel)[match(c("feature_x", "feature_y", "idx_x", "idx_y", "size_x", "size_y"),
                        names(r_ba_rel))] <-
    c("feature_y", "feature_x", "idx_y", "idx_x", "size_y", "size_x")
  r_ba_rel <- r_ba_rel[, names(r_ab)]

  expect_equal(sort_results(r_ab), sort_results(r_ba_rel))
})

test_that("run_coreact respects chunk_size and results are chunking-invariant", {
  withr::local_seed(99)
  o_x <- make_obj(50, 5, "X")
  o_y <- make_obj(50, 5, "Y")

  suppressMessages({
    res_10 <- run_coreact(o_x, o_y, chunk_size = 10, filter_config = list(min_intersection = 0))
    res_1 <- run_coreact(o_x, o_y, chunk_size = 1, filter_config = list(min_intersection = 0))
    res_7 <- run_coreact(o_x, o_y, chunk_size = 7, filter_config = list(min_intersection = 0))
    res_auto <- run_coreact(o_x, o_y, chunk_size = NULL, filter_config = list(min_intersection = 0))
  })

  expect_equal(nrow(res_10), 2500)
  expect_equal(sort_results(res_10), sort_results(res_1))
  expect_equal(sort_results(res_10), sort_results(res_7))
  expect_equal(sort_results(res_10), sort_results(res_auto))
})

test_that("run_coreact parallel execution matches serial", {
  skip_on_os("windows")
  skip_on_cran()
  withr::local_seed(2024)
  o_x <- make_obj(120, 8, "X", density = 0.3)
  o_y <- make_obj(40, 8, "Y", density = 0.5)

  suppressMessages({
    res_serial <- run_coreact(o_x, o_y, n_cores = 1, chunk_size = 25,
                              filter_config = list(min_intersection = 1))
    res_par <- run_coreact(o_x, o_y, n_cores = 2, chunk_size = 25,
                           filter_config = list(min_intersection = 1))
  })

  expect_gt(nrow(res_serial), 0)
  expect_equal(sort_results(res_serial), sort_results(res_par))
})

test_that("run_coreact on Windows falls back to serial with a warning", {
  skip_on_os(c("mac", "linux", "solaris"))
  withr::local_seed(1)
  o_x <- make_obj(10, 5, "X")
  o_y <- make_obj(10, 5, "Y")
  expect_warning(
    suppressMessages(run_coreact(o_x, o_y, n_cores = 2, filter_config = list())),
    "Running serially"
  )
})

test_that("run_coreact returns the empty schema when nothing passes or inputs are empty", {
  withr::local_seed(3)
  o_x <- make_obj(10, 5, "X")
  o_y <- make_obj(10, 5, "Y")

  suppressMessages(res <- run_coreact(o_x, o_y, filter_config = list(min_intersection = 100)))
  expect_equal(nrow(res), 0)
  expect_named(res, names(empty_result_tibble()))
  expect_type(res$feature_x, "character")
  expect_type(res$p_val, "double")

  # Zero-feature input (e.g. everything filtered upstream) -> empty schema, no error
  o_empty <- o_x[integer(0), ]
  expect_equal(nrow(o_empty), 0)
  suppressMessages(res_e <- run_coreact(o_empty, o_y, filter_config = list()))
  expect_equal(nrow(res_e), 0)
  expect_named(res_e, names(empty_result_tibble()))
})

test_that("run_coreact validates its inputs", {
  withr::local_seed(5)
  o_x <- make_obj(10, 5, "X")
  o_y <- make_obj(10, 5, "Y")

  expect_error(run_coreact(o_x$mat, o_y), "coreact_data")
  expect_error(run_coreact(o_x, list()), "coreact_data")

  expect_error(suppressMessages(run_coreact(o_x, o_y, chunk_size = 0)), "chunk_size")
  expect_error(suppressMessages(run_coreact(o_x, o_y, chunk_size = -5)), "chunk_size")
  expect_error(suppressMessages(run_coreact(o_x, o_y, chunk_size = c(1, 2))), "chunk_size")
  expect_error(suppressMessages(run_coreact(o_x, o_y, n_cores = 0)), "n_cores")

  # Misspelt filter names used to be silently ignored
  expect_error(run_coreact(o_x, o_y, filter_config = list(min_intersect = 2)), "Unknown 'filter_config'")
  expect_error(run_coreact(o_x, o_y, filter_config = list(2)), "must be named")
  expect_error(run_coreact(o_x, o_y, filter_config = list(min_jaccard = -1)), "non-negative")
  expect_error(run_coreact(o_x, o_y, filter_config = list(min_jaccard = "a")), "non-negative")
  expect_error(run_coreact(o_x, o_y, filter_config = "min_intersection"), "named list")
  expect_no_error(validate_filter_config(NULL))
  expect_no_error(validate_filter_config(list()))
})

test_that("run_coreact refuses mismatched or reordered samples", {
  withr::local_seed(6)
  o_x <- make_obj(10, 5, "X")
  o_y <- make_obj(10, 5, "Y")

  # Same names, different order
  o_y_rev <- o_y
  colnames(o_y_rev$mat) <- rev(colnames(o_y_rev$mat))
  expect_error(run_coreact(o_x, o_y_rev, filter_config = list()), "different order")

  # Different sample sets
  o_y_other <- o_y
  colnames(o_y_other$mat) <- paste0("T", 1:5)
  expect_error(run_coreact(o_x, o_y_other, filter_config = list()), "not identical")

  # Different number of samples
  o_y_short <- o_y[, 1:3]
  expect_error(run_coreact(o_x, o_y_short, filter_config = list()), "not identical")

  # Unnamed samples: only the count can be checked
  o_x_nn <- o_x
  o_y_nn <- o_y_short
  colnames(o_x_nn$mat) <- NULL
  colnames(o_y_nn$mat) <- NULL
  expect_error(run_coreact(o_x_nn, o_y_nn, filter_config = list()), "Sample count differs")
})

# --- Tests for coreact_pipeline (The Full Workflow) ---

test_that("coreact_pipeline runs end-to-end and returns what it wrote", {
  withr::local_seed(11)
  tmp <- withr::local_tempdir()
  path_x <- write_pipeline_tsv(10, 4, "GeneX", dir = tmp)
  path_y <- write_pipeline_tsv(10, 4, "GeneY", dir = tmp)
  out_path <- file.path(tmp, "final_results.tsv")

  expect_message(
    res <- coreact_pipeline(
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
  on_disk <- data.table::fread(out_path)
  expect_gt(nrow(on_disk), 0)
  expect_true("p_adj" %in% colnames(on_disk))

  # Returned (invisibly) tibble matches the file
  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), nrow(on_disk))
  expect_equal(names(res), names(on_disk))
  expect_equal(res$p_adj, on_disk$p_adj)
  expect_equal(res$p_adj, stats::p.adjust(res$p_val, method = "BH"))

  # Check Sidecar Metadata
  side_x <- file.path(tmp, "final_results_metadata_x.tsv")
  side_y <- file.path(tmp, "final_results_metadata_y.tsv")
  expect_true(file.exists(side_x))
  expect_true(file.exists(side_y))
  meta_x <- data.table::fread(side_x)
  expect_setequal(meta_x$feature_id, unique(res$feature_x))
  expect_named(meta_x, c("feature_id", "FeatureID", "Desc"))
})

test_that("coreact_pipeline applies the FDR threshold and prevalence filters", {
  withr::local_seed(12)
  tmp <- withr::local_tempdir()
  path_x <- write_pipeline_tsv(30, 12, "X", dir = tmp)
  path_y <- write_pipeline_tsv(30, 12, "Y", dir = tmp)

  suppressWarnings(suppressMessages({
    all_res <- coreact_pipeline(
      paths = c(path_x, path_y), out_path = file.path(tmp, "all.tsv"),
      meta_cols = 1:2, feature_ids = 1, fdr_threshold = NULL
    )
    strict <- coreact_pipeline(
      paths = c(path_x, path_y), out_path = file.path(tmp, "strict.tsv"),
      meta_cols = 1:2, feature_ids = 1, fdr_threshold = 0.5
    )
  }))
  expect_true(all(strict$p_adj <= 0.5))
  expect_equal(nrow(strict), sum(all_res$p_adj <= 0.5))

  # Prevalence filter reports counts and reduces the feature space
  expect_message(
    prev <- coreact_pipeline(
      paths = c(path_x, path_y), out_path = file.path(tmp, "prev.tsv"),
      meta_cols = 1:2, feature_ids = 1, min_prevalence = c(0.5, 8), fdr_threshold = NULL
    ),
    "Prevalence filter"
  )
  expect_true(all(prev$size_x >= 6))
  expect_true(all(prev$size_y >= 8))
})

test_that("coreact_pipeline supports sample_cols filtering", {
  withr::local_seed(13)
  tmp <- withr::local_tempdir()
  # Scenario: Two files have different samples, but share S3 and S4.
  # X: S1, S2, S3, S4
  # Y: S3, S4, S5, S6
  path_x <- write_pipeline_tsv(10, c("S1", "S2", "S3", "S4"), "X", dir = tmp)
  path_y <- write_pipeline_tsv(10, c("S3", "S4", "S5", "S6"), "Y", dir = tmp)
  out_path <- file.path(tmp, "filtered_results.tsv")

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
  expect_message(
    res <- coreact_pipeline(
      paths = c(path_x, path_y),
      out_path = out_path,
      meta_cols = c("FeatureID", "Desc"),
      feature_ids = 1,
      sample_cols = c("S3", "S4"),
      fdr_threshold = 1.0
    ),
    "Pipeline completed successfully"
  )

  expect_true(file.exists(out_path))
  expect_true(all(res$size_x <= 2))
})

test_that("coreact_pipeline validates mismatched inputs", {
  withr::local_seed(14)
  tmp <- withr::local_tempdir()
  p1 <- write_pipeline_tsv(5, 3, "A", dir = tmp)
  p2 <- write_pipeline_tsv(5, 4, "B", dir = tmp) # Different n columns

  expect_error(
    coreact_pipeline(
      paths = c(p1, p2),
      out_path = tempfile(tmpdir = tmp),
      meta_cols = 1:2,
      feature_ids = 1
    ),
    "Sample identifiers.*not identical"
  )
})

test_that("coreact_pipeline validates arguments before reading files", {
  tmp <- withr::local_tempdir()
  out <- file.path(tmp, "out.tsv")
  ok <- c("a.tsv", "b.tsv") # never read: validation fails first

  expect_error(coreact_pipeline(paths = 1:2, out_path = out), "'paths'")
  expect_error(coreact_pipeline(paths = c("a", "b", "c"), out_path = out), "'paths'")
  expect_error(coreact_pipeline(paths = ok, out_path = out, names = "only_one"), "'names'")
  expect_error(coreact_pipeline(paths = ok, out_path = c(out, out)), "'out_path'")
  expect_error(coreact_pipeline(paths = ok, out_path = out, binarize = "yes"), "'binarize'")
  expect_error(coreact_pipeline(paths = ok, out_path = out, binarize = c(TRUE, FALSE, TRUE)), "'binarize'")
  expect_error(coreact_pipeline(paths = ok, out_path = out, min_prevalence = NA), "'min_prevalence'")
  expect_error(coreact_pipeline(paths = ok, out_path = out, min_prevalence = 1:3), "'min_prevalence'")
  expect_error(coreact_pipeline(paths = ok, out_path = out, meta_cols = list(1, 2, 3)), "'meta_cols'")
  expect_error(coreact_pipeline(paths = ok, out_path = out,
                                filter_config = list(min_intersect = 1)), "Unknown 'filter_config'")
  # Only after all of the above does the file check run
  expect_error(coreact_pipeline(paths = ok, out_path = out), "File not found")
})

test_that("coreact_pipeline handles empty results gracefully", {
  tmp <- withr::local_tempdir()
  rows <- 5
  cols <- 4
  df <- as.data.frame(matrix(0, rows, cols))
  colnames(df) <- paste0("S", seq_len(cols))
  df$ID <- paste0("G", seq_len(rows))
  df <- df[, c("ID", paste0("S", seq_len(cols)))]

  p1 <- file.path(tmp, "zeros.tsv")
  write.table(df, p1, sep = "\t", row.names = FALSE, quote = FALSE)

  out_empty <- file.path(tmp, "empty_results.tsv")

  expect_warning(
    res <- coreact_pipeline(
      paths = c(p1, p1),
      out_path = out_empty,
      meta_cols = 1, feature_ids = 1,
      min_prevalence = 0,
      filter_config = list(min_intersection = 1)
    ),
    "No significant results found"
  )

  expect_true(file.exists(out_empty))
  on_disk <- data.table::fread(out_empty)
  expect_equal(nrow(on_disk), 0)
  # Same schema as a non-empty result, including p_adj
  expect_named(on_disk, c(names(empty_result_tibble()), "p_adj"))
  expect_equal(nrow(res), 0)
  expect_true("p_adj" %in% names(res))

  # Sidecars exist and are empty apart from the header
  meta_x <- data.table::fread(file.path(tmp, "empty_results_metadata_x.tsv"))
  expect_equal(nrow(meta_x), 0)
  expect_named(meta_x, c("feature_id", "ID"))
})

test_that("coreact_pipeline creates missing output directories", {
  withr::local_seed(15)
  tmp <- withr::local_tempdir()
  path_x <- write_pipeline_tsv(6, 4, "X", dir = tmp)
  out_path <- file.path(tmp, "nested", "deeper", "res.tsv")

  suppressWarnings(suppressMessages(
    coreact_pipeline(paths = c(path_x, path_x), out_path = out_path,
                     meta_cols = 1:2, feature_ids = 1, fdr_threshold = NULL)
  ))
  expect_true(file.exists(out_path))
  expect_true(file.exists(file.path(tmp, "nested", "deeper", "res_metadata_y.tsv")))
})
