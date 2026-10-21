# Interface of run_coreact_outcome(): inputs, checks, output shape, testing.

# --- Output shape ---

test_that("the output is rectangular, with the documented columns and types", {
  d <- outcome_fixture(seed = 1)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary))

  expect_identical(names(out), names(empty_outcome_tibble()))
  expect_identical(vapply(out, typeof, ""), vapply(empty_outcome_tibble(), typeof, ""))
  expect_identical(levels(out$occupancy_class), occupancy_levels)
  expect_identical(levels(out$contrast), contrast_catalogue)
  expect_identical(levels(out$beta_estimability), estimability_levels)
  expect_equal(ncol(out), 24)
  # Every contrast for every pair, none dropped
  expect_equal(nrow(out), 7 * nrow(d$res))
  # By default only tested rows get intervals
  expect_identical(attr(out, "ci"), "tested")
  expect_true(all(is.na(out$ci_lo[out$role == "reported"])))
})

test_that("only tested rows carry tests, and p_adj adjusts over them alone", {
  d <- outcome_fixture(seed = 2)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary,
                                              fdr_method = "BH"))
  reported <- out$role == "reported"
  expect_true(all(is.na(out$p_val[reported])))
  expect_true(all(is.na(out$p_adj[reported])))
  expect_true(all(is.na(out$statistic[reported])))
  expect_true(all(is.na(out$df[reported])))

  tested <- !reported
  expect_true(all(out$df[tested] == 1L))
  expect_equal(out$p_adj[tested], stats::p.adjust(out$p_val[tested], "BH"))
})

test_that("any_significant is per pair, and NA for a pair with nothing tested", {
  d <- outcome_fixture(seed = 3)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary))
  pair <- pair_index(out)

  untested <- out$occupancy_class %in% c("disjoint", "degenerate")
  expect_true(any(untested))
  expect_true(all(is.na(out$any_significant[untested])))
  expect_false(anyNA(out$any_significant[!untested]))

  expected <- as.vector(tapply(out$role == "tested" & !is.na(out$p_adj) & out$p_adj <= 0.05, pair, any))
  expect_identical(out$any_significant[!untested], expected[pair][!untested])
})

test_that("the family-size guard stops a run whose tested set has drifted", {
  d <- outcome_fixture(seed = 4)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary))
  pair <- pair_index(out)
  expect_no_error(adjust_outcome_fdr(out, pair))

  dropped <- out
  dropped$role[which(dropped$role == "tested")[1]] <- "reported"
  expect_error(adjust_outcome_fdr(dropped, pair), "Internal error")

  lost <- out
  lost$p_val[which(lost$role == "tested")[1]] <- NA
  expect_error(adjust_outcome_fdr(lost, pair), "Internal error")
})

test_that("the settings and coreact's provenance are recorded", {
  d <- outcome_fixture(seed = 5)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$continuous,
                                              outcome_type = "continuous", near_nested_max = 8,
                                              contrasts = c("interaction", "11_vs_00"),
                                              ci = "none", fdr_method = "holm"))
  expect_equal(attr(out, "near_nested_max"), 8)
  expect_identical(attr(out, "contrasts"), c("interaction", "11_vs_00"))
  expect_identical(attr(out, "ci"), "none")
  expect_identical(attr(out, "fdr_method"), "holm")
  expect_identical(attr(out, "outcome_type"), "continuous")
  expect_false(attr(out, "within_library"))
  expect_identical(attr(out, "coreact_provenance"), attr(d$res, "provenance"))
  expect_true(is.numeric(attr(out, "frac_unbounded_tested")))
})

test_that("a result with no pairs gives the empty table", {
  d <- outcome_fixture(seed = 6)
  none <- suppressMessages(run_coreact(d$x, d$y, effect_filters = list(min_intersection = 1e6)))
  out <- suppressMessages(run_coreact_outcome(none, d$x, d$y, outcome = d$binary))
  expect_equal(nrow(out), 0)
  expect_identical(vapply(out, typeof, ""), vapply(empty_outcome_tibble(), typeof, ""))
})

test_that("many separated tested rows are reported in a message", {
  # Twelve samples: nearly every cell is all cases or all controls
  f <- pair_fixture(cells4(3, 3, 3, 3), cases = cells4(3, 0, 0, 3))
  expect_message(run_coreact_outcome(f$res, f$x, f$y, outcome = f$outcome), "separated")
})

# --- Within a library ---

test_that("a within-library result is analysed with data_y omitted", {
  d <- outcome_fixture(seed = 7)
  within <- suppressMessages(run_coreact(d$x, effect_filters = list()))
  out <- suppressMessages(run_coreact_outcome(within, d$x, outcome = d$binary))
  expect_equal(nrow(out), 7 * nrow(within))
  expect_true(attr(out, "within_library"))
})

test_that("the declared mode must match the recorded one, in both directions", {
  d <- outcome_fixture(seed = 8)
  within <- suppressMessages(run_coreact(d$x, effect_filters = list()))
  expect_error(run_coreact_outcome(within, d$x, d$y, outcome = d$binary),
               "within-library coreact result: omit 'data_y'")
  expect_error(run_coreact_outcome(d$res, d$x, outcome = d$binary),
               "two-input coreact result: supply 'data_y'")
})

test_that("an edited within-library result with both orientations is refused", {
  d <- outcome_fixture(seed = 9)
  within <- suppressMessages(run_coreact(d$x, effect_filters = list()))
  flipped <- within
  flipped[1, c("idx_x", "idx_y", "feature_x", "feature_y")] <-
    within[1, c("idx_y", "idx_x", "feature_y", "feature_x")]
  expect_error(run_coreact_outcome(flipped, d$x, outcome = d$binary), "one orientation per pair")
})

test_that("a pair listed twice is refused", {
  d <- outcome_fixture(seed = 10)
  twice <- d$res[c(1, 1, 2), ]
  attr(twice, "provenance") <- attr(d$res, "provenance")
  expect_error(run_coreact_outcome(twice, d$x, d$y, outcome = d$binary), "more than once")
})

test_that("between two datasets with row-number IDs, (1, 2) and (2, 1) are different pairs", {
  withr::local_seed(11)
  make <- function(name) {
    m <- matrix(stats::rbinom(4 * 60, 1, 0.4), 4, dimnames = list(as.character(1:4), paste0("S", 1:60)))
    new_coreact_data(m, data.frame(row.names = rownames(m)), name = name)
  }
  a <- make("A")
  b <- make("B")
  res <- suppressMessages(run_coreact(a, b, effect_filters = list()))
  o <- stats::setNames(stats::rbinom(60, 1, 0.5), paste0("S", 1:60))
  out <- suppressMessages(run_coreact_outcome(res, a, b, outcome = o))
  pairs <- unique(paste(out$feature_x, out$feature_y))
  expect_true(all(c("1 2", "2 1") %in% pairs))
})

# --- The objects must be those the result was computed on ---

test_that("data_x and data_y swapped are caught even under row-number IDs", {
  # The index check alone passes here: both sides have rows "1" to "4"
  withr::local_seed(12)
  make <- function(name, p) {
    m <- matrix(stats::rbinom(4 * 60, 1, p), 4, dimnames = list(as.character(1:4), paste0("S", 1:60)))
    new_coreact_data(m, data.frame(row.names = rownames(m)), name = name)
  }
  a <- make("A", 0.3)
  b <- make("B", 0.5)
  res <- suppressMessages(run_coreact(a, b, effect_filters = list()))
  o <- stats::setNames(stats::rbinom(60, 1, 0.5), paste0("S", 1:60))
  expect_error(run_coreact_outcome(res, b, a, outcome = o), "is not the data this result was computed on")
})

test_that("an object with the recorded name and counts but other contents is refused", {
  # One positive moved to another sample: name, size and number of positives
  # all match the record, so only the content hash can tell
  d <- outcome_fixture(seed = 22)
  mat <- as.matrix(d$x$mat)
  row <- which(rowSums(mat) > 0 & rowSums(mat) < ncol(mat))[1]
  from <- which(mat[row, ] == 1)[1]
  to <- which(mat[row, ] == 0)[1]
  mat[row, c(from, to)] <- c(0, 1)
  edited <- new_coreact_data(mat, d$x$meta, name = d$x$name)
  expect_equal(sum(edited$mat), sum(d$x$mat))
  expect_error(run_coreact_outcome(d$res, edited, d$y, outcome = d$binary),
               "same name and counts as recorded, but different contents")
})

test_that("objects filtered differently from those of the result are refused", {
  d <- outcome_fixture(seed = 13)
  expect_error(run_coreact_outcome(d$res, d$x[-1, ], d$y, outcome = d$binary), "data_x")

  # Same name, size and positive cells, rows in another order: the content
  # hash tells
  perm <- rev(seq_len(nrow(d$x$mat)))
  reordered <- new_coreact_data(d$x$mat[perm, ], d$x$meta[perm, , drop = FALSE], name = d$x$name)
  expect_error(run_coreact_outcome(d$res, reordered, d$y, outcome = d$binary),
               "different contents")
})

test_that("a result whose indices were edited is refused, naming the pair", {
  # The objects are right; the result is not: two rows' idx_x exchanged
  d <- outcome_fixture(seed = 23)
  edited <- d$res
  rows <- which(edited$idx_x != edited$idx_x[1])[1]
  edited$idx_x[c(1, rows)] <- edited$idx_x[c(rows, 1)]
  expect_error(run_coreact_outcome(edited, d$x, d$y, outcome = d$binary), "indices point at")
})

test_that("the pairs analysed are the ones named in the result", {
  d <- outcome_fixture(seed = 14)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary))
  first <- out[!duplicated(pair_index(out)), ]
  for (i in seq_len(nrow(first))) {
    x <- as.numeric(d$x$mat[first$feature_x[i], ])
    y <- as.numeric(d$y$mat[first$feature_y[i], ])
    expect_identical(c(first$n00[i], first$n10[i], first$n01[i], first$n11[i]),
                     as.integer(c(sum(!x & !y), sum(x & !y), sum(!x & y), sum(x & y))))
  }
})

test_that("a result without its provenance is refused, with the way to restore it", {
  d <- outcome_fixture(seed = 15)
  bare <- d$res
  attr(bare, "provenance") <- NULL
  expect_error(run_coreact_outcome(bare, d$x, d$y, outcome = d$binary), "read_coreact_results")
})

test_that("a result read back from disk gives the same table as in memory", {
  withr::local_seed(16)
  dir <- withr::local_tempdir()
  path <- write_pipeline_tsv(12, 50, prefix = "F", dir = dir)
  out_path <- file.path(dir, "res.tsv")
  in_memory <- suppressMessages(coreact_pipeline(path, out_path, meta_cols = 1:2, fdr_threshold = 1))
  from_disk <- read_coreact_results(out_path)
  obj <- suppressMessages(import_coreact_tsv(path, meta_cols = 1:2, name = "Input"))
  o <- stats::setNames(stats::rbinom(50, 1, 0.5), colnames(obj$mat))

  a <- suppressMessages(run_coreact_outcome(in_memory, obj, outcome = o))
  b <- suppressMessages(run_coreact_outcome(from_disk, obj, outcome = o))
  expect_identical(a, b)
})

# --- Arguments and the outcome ---

test_that("positional and repeated data objects are refused with the remedy", {
  d <- outcome_fixture(seed = 17)
  within <- suppressMessages(run_coreact(d$x, effect_filters = list()))
  # run_coreact_outcome(res, x, outcome) binds the outcome to data_y
  expect_error(run_coreact_outcome(within, d$x, d$binary), "pass 'outcome' by name")
  expect_error(run_coreact_outcome(d$res, d$x, d$x, outcome = d$binary), "hold the same data")
})

test_that("the outcome is matched by sample ID and validated", {
  d <- outcome_fixture(seed = 18)
  run <- function(o, type = "binary") {
    suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = o, outcome_type = type))
  }
  o <- d$binary
  expect_error(run(unname(o)), "named by sample ID")
  expect_error(run(o[-1]), "no value for S1")
  expect_error(run(c(o, S999 = 1)), "unknown samples S999")
  expect_error(run(c(o, o[1])), "more than once")
  expect_error(run(replace(o, 1, 2)), "coded 0/1")
  expect_error(run(stats::setNames(factor(c("a", "b", "c"))[(seq_along(o) %% 3) + 1], names(o))),
               "two levels")

  # Matching is by name: reversing the vector changes nothing
  expect_identical(run(rev(o)), run(o))
  # Logical and factor codings agree with 0/1; the factor's second level is 1
  expect_identical(run(o == 1), run(o))
  fct <- stats::setNames(factor(ifelse(o == 1, "case", "control"), levels = c("control", "case")), names(o))
  expect_identical(run(fct), run(o))
})

test_that("the outcome is checked on the analysis sample", {
  d <- outcome_fixture(seed = 19)
  run <- function(o, type) {
    run_coreact_outcome(d$res, d$x, d$y, outcome = o, outcome_type = type)
  }
  b <- d$binary
  # Only one value left once missing outcomes are dropped
  one_value <- b
  one_value[b == 0] <- NA
  expect_error(run(one_value, "binary"), "single value")
  expect_error(run(stats::setNames(rep(NA_integer_, length(b)), names(b)), "binary"), "non-missing")

  cont <- d$continuous
  expect_error(run(replace(cont, 1, Inf), "continuous"), "finite")
  expect_error(run(stats::setNames(rep(2, length(cont)), names(cont)), "continuous"), "constant")
  few <- cont
  few[-(1:4)] <- NA
  expect_error(run(few, "continuous"), "at least 5")
  expect_error(run(stats::setNames(as.character(cont), names(cont)), "continuous"), "numeric")
  expect_warning(suppressMessages(run(b + 0, "continuous")), "only two values")
})

test_that("invalid arguments are refused before any fitting", {
  d <- outcome_fixture(seed = 20)
  run <- function(...) run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary, ...)
  expect_error(run(ci = "bogus"), "should be one of")
  expect_error(run(covariates = data.frame(a = 1)), "covariates")
  expect_error(run(contrasts = "12_vs_00"), "Unknown contrast")
  expect_error(run(contrasts = character(0)), "non-empty")
  expect_error(run(fdr_method = "bogus"), "fdr_method")
  expect_error(run(near_nested_max = -1), "near_nested_max")
  expect_error(run(near_nested_max = function(n) NA), "near_nested_max")
  expect_error(run(n_cores = 0), "n_cores")
  expect_error(run(chunk_size = 0), "chunk_size")
  expect_error(run(outcome_type = "count"), "should be one of")
  expect_error(run_coreact_outcome(as.list(d$res), d$x, d$y, outcome = d$binary), "coreact result")
})

test_that("pairs that fall below their selection thresholds are counted", {
  d <- outcome_fixture(seed = 21)
  filtered <- suppressMessages(filter_by_prevalence(d$x, 40))
  res <- suppressMessages(run_coreact(filtered, d$y, effect_filters = list()))
  o <- d$binary
  o[1:60] <- NA                           # drops enough samples to matter

  expect_message(out <- run_coreact_outcome(res, filtered, d$y, outcome = o),
                 "fall below the prevalence or intersection threshold")
  first <- out[!duplicated(pair_index(out)), ]
  expect_identical(attr(out, "n_pairs_below_selection"),
                   as.integer(sum(first$n10 + first$n11 < 40)))
  expect_gt(attr(out, "n_pairs_below_selection"), 0)
})

test_that("the in-memory flow filters, adjusts and passes only the co-reactive pairs on", {
  # filter_by_prevalence() is exported, takes the pipeline's argument names,
  # and its bounds reach the provenance
  expect_true("filter_by_prevalence" %in% getNamespaceExports("coreact"))
  d <- outcome_fixture(n_samples = 200, seed = 24)
  x <- suppressMessages(filter_by_prevalence(d$x, min_prevalence = 0.1))
  y <- suppressMessages(filter_by_prevalence(d$y, min_prevalence = 20, max_prevalence = 0.9))
  res <- adjust_fdr(suppressMessages(run_coreact(x, y, effect_filters = list())))
  prov <- attr(res, "provenance")
  expect_equal(prov$x_min_prevalence_count, 20)
  expect_equal(c(prov$y_min_prevalence_count, prov$y_max_prevalence_count), c(20, 180))

  # Row-filtering on the co-occurrence FDR keeps the provenance, and the
  # outcome step analyses exactly the pairs kept
  sig <- res[res$p_adj <= 0.05, ]
  expect_gt(nrow(sig), 0)
  expect_lt(nrow(sig), nrow(res))
  out <- suppressMessages(run_coreact_outcome(sig, x, y, outcome = d$binary))
  expect_setequal(unique(paste(out$feature_x, out$feature_y)), paste(sig$feature_x, sig$feature_y))
})
