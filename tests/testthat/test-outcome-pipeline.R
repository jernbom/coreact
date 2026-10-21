# coreact_outcome_pipeline(): from the files coreact_pipeline() wrote to the
# outcome table, by way of the provenance sidecar.

# Write a one-row-per-sample outcome TSV, with an unrelated column in between
write_outcome_tsv <- function(ids, values, dir, name = "outcome.tsv") {
  path <- file.path(dir, name)
  utils::write.table(data.frame(sample = ids, age = seq_along(ids), status = values),
                     path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  path
}

# A two-input coreact_pipeline() run on disk, with per-side settings that the
# rebuild has to reproduce: different prevalence bounds, one of them a fraction
pipeline_fixture <- function(dir, seed = 1, n_samples = 80) {
  withr::with_seed(seed, {
    m <- make_signal_matrices(10, 6, n_samples, density = 0.3, flip = 0.2)
    path_x <- write_matrix_tsv(m$x, dir)
    path_y <- write_matrix_tsv(m$y, dir)
    status <- stats::rbinom(n_samples, 1, stats::plogis(-0.5 + m$x[1, ] + 0.5 * m$y[2, ]))
  })
  results <- file.path(dir, "cooc.tsv")
  res <- suppressMessages(coreact_pipeline(c(path_x, path_y), results, meta_cols = 1:2,
                                           feature_ids = 1, min_prevalence = c(0.1, 12),
                                           fdr_threshold = 1))
  ids <- colnames(m$x)
  list(path_x = path_x, path_y = path_y, results = results, res = res, ids = ids,
       outcome = stats::setNames(status, ids))
}

# The objects coreact_pipeline() analysed, rebuilt by hand
in_memory_inputs <- function(path, min_prevalence) {
  obj <- suppressMessages(import_coreact_tsv(path, meta_cols = 1:2, feature_id = 1,
                                             name = if (min_prevalence < 1) "Input_X" else "Input_Y"))
  suppressMessages(filter_by_prevalence(obj, min_prevalence))
}

test_that("the disk route reproduces the in-memory analysis exactly", {
  dir <- withr::local_tempdir()
  f <- pipeline_fixture(dir)
  outcome_file <- write_outcome_tsv(f$ids, f$outcome, dir)
  out_path <- file.path(dir, "assoc.tsv")

  from_disk <- suppressMessages(coreact_outcome_pipeline(f$results, outcome_file, out_path,
                                                         outcome_col = "status",
                                                         sample_id_col = "sample"))

  x <- in_memory_inputs(f$path_x, 0.1)
  y <- in_memory_inputs(f$path_y, 12)
  in_memory <- suppressMessages(run_coreact_outcome(f$res, x, y, outcome = f$outcome))
  expect_identical(from_disk, in_memory)
})

test_that("the rebuilt inputs are the objects coreact analysed, indices and all", {
  dir <- withr::local_tempdir()
  f <- pipeline_fixture(dir, seed = 2)
  res <- read_coreact_results(f$results)
  rebuilt <- suppressMessages(rebuild_inputs(attr(res, "provenance")))

  expect_identical(rebuilt$x, in_memory_inputs(f$path_x, 0.1))
  expect_identical(rebuilt$y, in_memory_inputs(f$path_y, 12))
  expect_identical(rownames(rebuilt$x$mat)[res$idx_x], res$feature_x)
  expect_identical(rownames(rebuilt$y$mat)[res$idx_y], res$feature_y)
  # And the results read back can be adjusted again
  expect_equal(adjust_fdr(res)$p_adj, f$res$p_adj)
})

test_that("a within-library result is rebuilt from its single input", {
  dir <- withr::local_tempdir()
  withr::local_seed(3)
  path <- write_pipeline_tsv(12, 70, prefix = "F", dir = dir)
  results <- file.path(dir, "within.tsv")
  res <- suppressMessages(coreact_pipeline(path, results, meta_cols = 1:2, fdr_threshold = 1))
  obj <- suppressMessages(import_coreact_tsv(path, meta_cols = 1:2, name = "Input"))
  o <- stats::setNames(stats::rbinom(70, 1, 0.5), colnames(obj$mat))
  outcome_file <- write_outcome_tsv(names(o), o, dir)

  from_disk <- suppressMessages(coreact_outcome_pipeline(results, outcome_file,
                                                         file.path(dir, "assoc.tsv"),
                                                         outcome_col = "status"))
  in_memory <- suppressMessages(run_coreact_outcome(res, obj, outcome = o))
  expect_identical(from_disk, in_memory)
  expect_true(attr(from_disk, "within_library"))
})

test_that("both files are written, and the settings read back", {
  dir <- withr::local_tempdir()
  f <- pipeline_fixture(dir, seed = 4)
  outcome_file <- write_outcome_tsv(f$ids, f$outcome, dir)
  out_path <- file.path(dir, "assoc.tsv")
  out <- suppressMessages(coreact_outcome_pipeline(f$results, outcome_file, out_path,
                                                   outcome_col = 3, ci = "none",
                                                   near_nested_max = function(n) n / 8))

  on_disk <- data.table::fread(out_path)
  expect_equal(nrow(on_disk), nrow(out))
  expect_identical(names(on_disk), names(out))

  prov <- read_provenance_sidecar(sidecar_path(out_path, "provenance"))
  expect_equal(prov$near_nested_max, 80 / 8)          # the resolved value, not the function
  expect_identical(prov$ci, "none")
  expect_identical(prov$contrasts, contrast_catalogue)
  expect_identical(prov$outcome_col, 3)
  expect_identical(prov$sample_id_col, 1)
  expect_identical(prov$results_path, normalizePath(f$results, winslash = "/"))
  expect_identical(prov$n_used, 80L)
})

test_that("columns chosen by index give the same result as by name", {
  dir <- withr::local_tempdir()
  f <- pipeline_fixture(dir, seed = 5)
  outcome_file <- write_outcome_tsv(f$ids, f$outcome, dir)
  run <- function(...) {
    suppressMessages(coreact_outcome_pipeline(f$results, outcome_file, file.path(dir, "a.tsv"),
                                              ci = "none", ...))
  }
  expect_identical(run(outcome_col = 3, sample_id_col = 1),
                   run(outcome_col = "status", sample_id_col = "sample"))
})

test_that("sample IDs are read as text, so leading zeros survive", {
  dir <- withr::local_tempdir()
  withr::local_seed(6)
  ids <- sprintf("%03d", 1:40)                        # "001", "002", ...
  m <- matrix(stats::rbinom(6 * 40, 1, 0.4), 6, dimnames = list(paste0("F", 1:6), ids))
  path <- write_matrix_tsv(m, dir)
  results <- file.path(dir, "r.tsv")
  suppressMessages(coreact_pipeline(path, results, meta_cols = 1:2, feature_ids = 1, fdr_threshold = 1))
  outcome_file <- write_outcome_tsv(ids, stats::rbinom(40, 1, 0.5), dir)

  out <- suppressMessages(coreact_outcome_pipeline(results, outcome_file, file.path(dir, "a.tsv"),
                                                   outcome_col = "status", ci = "none"))
  expect_equal(out$n_used[1], 40L)
})

test_that("a continuous outcome and missing values come through the file", {
  dir <- withr::local_tempdir()
  f <- pipeline_fixture(dir, seed = 7)
  values <- withr::with_seed(8, stats::rnorm(length(f$ids)))
  values[c(2, 9)] <- NA                                # empty cells in the TSV
  outcome_file <- write_outcome_tsv(f$ids, values, dir)

  from_disk <- suppressMessages(coreact_outcome_pipeline(f$results, outcome_file,
                                                         file.path(dir, "a.tsv"),
                                                         outcome_col = "status",
                                                         outcome_type = "continuous"))
  expect_equal(from_disk$n_used[1], length(f$ids) - 2L)
  x <- in_memory_inputs(f$path_x, 0.1)
  y <- in_memory_inputs(f$path_y, 12)
  # write.table() keeps 15 significant digits; compare against what the file holds
  stored <- as.numeric(as.character(values))
  in_memory <- suppressMessages(run_coreact_outcome(f$res, x, y, outcome = stats::setNames(stored, f$ids),
                                                    outcome_type = "continuous"))
  expect_identical(from_disk, in_memory)
})

test_that("problems with the files are reported before anything is fitted", {
  dir <- withr::local_tempdir()
  f <- pipeline_fixture(dir, seed = 9)
  outcome_file <- write_outcome_tsv(f$ids, f$outcome, dir)
  out_path <- file.path(dir, "a.tsv")
  run <- function(results = f$results, outcome = outcome_file, out = out_path, ...) {
    coreact_outcome_pipeline(results, outcome, out, ...)
  }

  expect_error(run(outcome_col = "nope"), "column 'nope' not found")
  expect_error(run(outcome_col = 9), "has 3 columns")
  expect_error(run(outcome_col = "sample"), "the same column")
  expect_error(run(outcome_col = c("a", "b")), "single column name")
  expect_error(run(outcome = file.path(dir, "absent.tsv"), outcome_col = "status"), "File not found")
  expect_error(run(out = f$results, outcome_col = "status"), "results file itself")

  text <- write_outcome_tsv(f$ids, ifelse(f$outcome == 1, "case", "control"), dir, "text.tsv")
  expect_error(run(outcome = text, outcome_col = "status"), "holds text")

  short <- write_outcome_tsv(f$ids[-1], f$outcome[-1], dir, "short.tsv")
  expect_error(suppressMessages(run(outcome = short, outcome_col = "status")), "no value for")

  # Results without their sidecar cannot be rebuilt
  bare <- file.path(dir, "bare.tsv")
  file.copy(f$results, bare)
  expect_error(run(results = bare, outcome_col = "status"), "No provenance sidecar")
})

test_that("an input that has moved is named in the error", {
  dir <- withr::local_tempdir()
  f <- pipeline_fixture(dir, seed = 10)
  outcome_file <- write_outcome_tsv(f$ids, f$outcome, dir)
  file.rename(f$path_y, file.path(dir, "moved.tsv"))
  expect_error(suppressMessages(coreact_outcome_pipeline(f$results, outcome_file,
                                                         file.path(dir, "a.tsv"),
                                                         outcome_col = "status")),
               "Y input .* is no longer at")
})
