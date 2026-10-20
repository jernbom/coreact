# Comparing a library with itself: one orientation per pair, no self-pairs,
# and a family of n * (n - 1) / 2.

test_that("run_coreact(obj) keeps one orientation per pair and no self-pairs", {
  withr::with_seed(42, {
    obj <- make_obj(12, 40, prefix = "F", density = 0.4)
    res <- suppressMessages(run_coreact(obj, effect_filters = list()))
  })

  n <- nrow(obj$mat)
  expect_equal(nrow(res), n * (n - 1) / 2)
  expect_true(all(res$idx_x < res$idx_y))
  expect_false(any(res$feature_x == res$feature_y))
  expect_equal(attr(res, "n_pairs"), n * (n - 1) / 2)

  # Every unordered pair present exactly once
  key <- paste(pmin(res$idx_x, res$idx_y), pmax(res$idx_x, res$idx_y))
  expect_equal(anyDuplicated(key), 0L)
})

test_that("within-library statistics equal the upper triangle of all pairs", {
  withr::with_seed(7, {
    obj <- make_obj(9, 30, prefix = "F", density = 0.5)
    res <- suppressMessages(run_coreact(obj, effect_filters = list()))
  })

  ref <- brute_force_coreact(obj$mat, obj$mat)
  ref <- ref[ref$idx_x < ref$idx_y, ]
  ref <- ref[order(ref$idx_x, ref$idx_y), ]
  got <- res[order(res$idx_x, res$idx_y), ]

  expect_equal(got$idx_x, ref$idx_x)
  expect_equal(got$idx_y, ref$idx_y)
  expect_equal(got$intersection, ref$intersection)
  expect_equal(got$union, ref$union)
  expect_equal(got$jaccard, ref$jaccard)
  expect_equal(got$overlap, ref$overlap)
  expect_equal(got$p_val, ref$p_val)
})

test_that("n_pairs counts the features retained, not those supplied", {
  withr::with_seed(3, obj <- make_obj(20, 30, prefix = "F", density = 0.1))

  filtered <- suppressMessages(filter_by_prevalence(obj, 5))
  n_kept <- nrow(filtered$mat)
  expect_lt(n_kept, nrow(obj$mat))        # the filter must actually bite

  res <- suppressMessages(run_coreact(filtered, effect_filters = list()))
  expect_equal(attr(res, "n_pairs"), n_kept * (n_kept - 1) / 2)
  expect_equal(nrow(res), n_kept * (n_kept - 1) / 2)
})

test_that("the kept triangle does not depend on chunk size or worker count", {
  withr::with_seed(11, obj <- make_obj(37, 40, prefix = "F", density = 0.4))

  runs <- lapply(c(1, 5, 13, 37, 100), function(cs) {
    sort_results(suppressMessages(
      run_coreact(obj, chunk_size = cs, effect_filters = list())
    ))
  })
  for (r in runs[-1]) expect_equal(r, runs[[1]])

  skip_on_os("windows")
  par <- sort_results(suppressMessages(
    run_coreact(obj, n_cores = 2, chunk_size = 7, effect_filters = list())
  ))
  expect_equal(par, runs[[1]])
})

test_that("a library with fewer than two features yields the empty schema", {
  withr::with_seed(5, obj <- make_obj(4, 20, prefix = "F", density = 0.5))

  for (n_feat in c(0, 1)) {
    small <- obj[seq_len(n_feat), ]
    res <- suppressMessages(run_coreact(small, effect_filters = list()))
    expect_equal(nrow(res), 0)
    expect_equal(attr(res, "n_pairs"), 0)
    expect_named(res, names(suppressMessages(
      run_coreact(obj, effect_filters = list(min_intersection = 1e6))
    )))
    # The family is empty, so the adjustment has nothing to spread over
    expect_equal(nrow(adjust_fdr(res)), 0)
  }
})

test_that("supplying one library as both arguments is refused", {
  withr::with_seed(9, obj <- make_obj(6, 20, prefix = "F", density = 0.5))

  expect_error(suppressMessages(run_coreact(obj, obj, effect_filters = list())),
               "same data")
  expect_error(suppressMessages(run_coreact(obj, obj, effect_filters = list())),
               "'data_x' only")

  # A separately built but equal object is the same data too
  same <- new_coreact_data(obj$mat, obj$meta, name = "other_name")
  expect_error(suppressMessages(run_coreact(obj, same, effect_filters = list())),
               "same data")

  # data_y must still be a coreact_data when it is supplied
  expect_error(run_coreact(obj, list()), "coreact_data")
})

test_that("coreact_pipeline refuses the same file twice and accepts one path", {
  withr::with_seed(13, {
    dir <- withr::local_tempdir()
    path <- write_pipeline_tsv(14, 25, prefix = "F", dir = dir)
    out <- file.path(dir, "within.tsv")

    # The same path twice, and a byte-identical copy under another name
    copy <- file.path(dir, "copy.tsv")
    file.copy(path, copy)

    expect_error(
      suppressMessages(coreact_pipeline(paths = c(path, path), out_path = out,
                                        meta_cols = 1:2, fdr_threshold = 1)),
      "same file"
    )
    expect_error(
      suppressMessages(coreact_pipeline(paths = c(path, copy), out_path = out,
                                        meta_cols = 1:2, fdr_threshold = 1)),
      "same data"
    )

    res <- suppressMessages(coreact_pipeline(paths = path, out_path = out,
                                             meta_cols = 1:2, fdr_threshold = 1))
  })

  expect_true(all(res$idx_x < res$idx_y))
  expect_false(any(res$feature_x == res$feature_y))
  expect_equal(nrow(res), 14 * 13 / 2)
  expect_true(file.exists(out))
  # Both sidecars are written, from the one library
  expect_true(file.exists(file.path(dir, "within_metadata_x.tsv")))
  expect_true(file.exists(file.path(dir, "within_metadata_y.tsv")))
})

test_that("one path requires single-valued per-dataset arguments", {
  withr::with_seed(17, {
    dir <- withr::local_tempdir()
    path <- write_pipeline_tsv(6, 20, prefix = "F", dir = dir)
    out <- file.path(dir, "within.tsv")
  })

  expect_error(
    coreact_pipeline(paths = path, out_path = out, meta_cols = 1:2,
                     min_prevalence = c(1, 2)),
    "'min_prevalence'"
  )
  expect_error(
    coreact_pipeline(paths = path, out_path = out,
                     meta_cols = list(1:2, 1:2)),
    "'meta_cols'"
  )
  expect_error(
    coreact_pipeline(paths = path, out_path = out, meta_cols = 1:2,
                     names = c("A", "B")),
    "length 1"
  )
  expect_error(
    coreact_pipeline(paths = character(0), out_path = out),
    "length 1 or 2"
  )
})


# One library supplied as two objects cannot be refused, because the objects
# differ. It is warned about instead.

test_that("sharing identical features across two inputs warns", {
  withr::with_seed(21, obj <- make_obj(16, 40, prefix = "F", density = 0.4))

  lo <- suppressMessages(filter_by_prevalence(obj, 2))
  hi <- suppressMessages(filter_by_prevalence(obj, 20))
  shared <- intersect(rownames(lo$mat), rownames(hi$mat))
  expect_gt(length(shared), 0)

  expect_warning(suppressMessages(run_coreact(lo, hi, effect_filters = list())),
                 sprintf("share %d feature", length(shared)))
  expect_warning(suppressMessages(run_coreact(lo, hi, effect_filters = list())),
                 "'data_x' alone")
})

test_that("the warning ignores uninformative and non-identical matches", {
  n_s <- 20
  ids <- as.character(1:4)

  # Row 1: identical, a single positive -> uninformative under row-number IDs
  # Row 2: identical, all zero          -> uninformative
  # Row 3: identical, prevalence 3      -> the only one that should count
  # Row 4: differs in exactly one sample
  mk <- function(row4_flip) {
    m <- matrix(0, 4, n_s, dimnames = list(ids, paste0("S", seq_len(n_s))))
    m[1, 1] <- 1
    m[3, 1:3] <- 1
    m[4, 5:8] <- 1
    if (row4_flip) m[4, 8] <- 0
    new_coreact_data(m, data.frame(row.names = ids), name = "lib")
  }
  a <- mk(FALSE)
  b <- mk(TRUE)

  expect_warning(suppressMessages(run_coreact(a, b, effect_filters = list())),
                 "share 1 feature")

  # Dropping the one informative shared row leaves nothing to warn about
  expect_silent(suppressMessages(
    run_coreact(a[c(1, 2, 4), ], b[c(1, 2, 4), ], effect_filters = list())
  ))
})

test_that("two libraries with no shared feature IDs never warn", {
  withr::with_seed(23, {
    o_x <- make_obj(8, 30, prefix = "X", density = 0.4)
    o_y <- make_obj(8, 30, prefix = "Y", density = 0.4)
  })
  expect_silent(suppressMessages(run_coreact(o_x, o_y, effect_filters = list())))
})
