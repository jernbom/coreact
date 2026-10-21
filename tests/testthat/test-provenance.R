# Provenance: how a result was computed, recorded in memory and in a sidecar,
# so that a result read back from disk can be adjusted and reproduced.

# --- Cutoffs travel on the filtered object ---

test_that("filter_by_prevalence records the bounds it applied, as counts", {
  withr::with_seed(1, obj <- make_obj(30, 40, density = 0.4))

  frac <- suppressMessages(filter_by_prevalence(obj, 0.25, 0.75))
  # 0.25 * 40 = 10 exactly; 0.75 * 40 = 30 exactly
  expect_equal(attr(frac, "prevalence_cutoffs"), c(min = 10, max = 30))

  odd <- suppressMessages(filter_by_prevalence(obj, 0.26, 0.74))
  # ceiling(10.4) = 11 for the lower bound, floor(29.6) = 29 for the upper
  expect_equal(attr(odd, "prevalence_cutoffs"), c(min = 11, max = 29))

  counts <- suppressMessages(filter_by_prevalence(obj, 3, NULL))
  expect_equal(attr(counts, "prevalence_cutoffs"), c(min = 3, max = Inf))
})

test_that("an unfiltered pass records no filter, not unknown", {
  withr::with_seed(2, obj <- make_obj(5, 10))
  same <- filter_by_prevalence(obj, 0, NULL)
  expect_equal(attr(same, "prevalence_cutoffs"), c(min = 0, max = Inf))
})

test_that("cutoffs are unknown once the object no longer matches them", {
  withr::with_seed(3, obj <- make_obj(10, 30, prefix = "F", density = 0.4))
  filtered <- suppressMessages(filter_by_prevalence(obj, 2))

  # Subsetting rebuilds the object, so the cutoffs no longer describe it
  expect_null(attr(filtered[1:5, ], "prevalence_cutoffs"))

  res <- suppressMessages(run_coreact(filtered[1:5, ], effect_filters = list()))
  prov <- attr(res, "provenance")
  expect_true(is.na(prov$x_min_prevalence_count))
  expect_true(is.na(prov$x_max_prevalence_count))

  # An object never filtered at all is unknown too, never "no filter"
  raw <- attr(suppressMessages(run_coreact(obj, effect_filters = list())), "provenance")
  expect_true(is.na(raw$x_min_prevalence_count))
})

# --- run_coreact() records provenance ---

test_that("the fingerprint describes each side as supplied, before the swap", {
  withr::with_seed(4, {
    o_x <- make_obj(4, 30, prefix = "X", density = 0.5, name = "small")
    o_y <- make_obj(9, 30, prefix = "Y", density = 0.5, name = "large")
  })
  # nrow(x) < nrow(y), so the engine swaps the two internally
  prov <- attr(suppressMessages(run_coreact(o_x, o_y, effect_filters = list())), "provenance")

  expect_identical(prov$x_name, "small")
  expect_identical(prov$y_name, "large")
  expect_equal(prov$x_n_features, 4)
  expect_equal(prov$y_n_features, 9)
  expect_equal(prov$x_n_positive, sum(o_x$mat))
  expect_equal(prov$y_n_positive, sum(o_y$mat))
  expect_false(prov$within_library)
})

test_that("provenance records the family, samples, filters and version", {
  withr::with_seed(5, {
    o_x <- make_obj(6, 25, prefix = "X")
    o_y <- make_obj(5, 25, prefix = "Y")
  })
  res <- suppressMessages(run_coreact(o_x, o_y, max_p_val = 0.5,
                                      effect_filters = list(min_jaccard = 0.1)))
  prov <- attr(res, "provenance")

  expect_equal(prov$n_pairs, 30)
  expect_equal(prov$n_samples, 25)
  expect_equal(prov$max_p_val, 0.5)
  expect_identical(prov$coreact_version, as.character(utils::packageVersion("coreact")))
  # All three effect filters are recorded; a missing one is 0, i.e. off
  expect_equal(prov[c("effect_min_intersection", "effect_min_jaccard", "effect_min_overlap")],
               list(effect_min_intersection = 0, effect_min_jaccard = 0.1,
                    effect_min_overlap = 0))

  # No pruning is recorded as NULL, not as a number
  open <- attr(suppressMessages(run_coreact(o_x, o_y, effect_filters = list())), "provenance")
  expect_true("max_p_val" %in% names(open))
  expect_null(open$max_p_val)
})

test_that("the n_pairs attribute is kept exactly as before", {
  withr::with_seed(6, {
    o_x <- make_obj(7, 20, prefix = "X")
    o_y <- make_obj(3, 20, prefix = "Y")
  })
  res <- suppressMessages(run_coreact(o_x, o_y, effect_filters = list(min_intersection = 1e6)))
  # Even an empty result carries both
  expect_equal(nrow(res), 0)
  expect_equal(attr(res, "n_pairs"), 21)
  expect_equal(attr(res, "provenance")$n_pairs, 21)
  expect_equal(nrow(adjust_fdr(res)), 0)
})

test_that("a within-library result records one side only", {
  withr::with_seed(7, obj <- make_obj(8, 20, prefix = "F"))
  prov <- attr(suppressMessages(run_coreact(obj, effect_filters = list())), "provenance")

  expect_true(prov$within_library)
  expect_equal(prov$n_pairs, 8 * 7 / 2)
  expect_false(any(startsWith(names(prov), "y_")))
})

# --- The sidecar ---

test_that("sidecars are named after the results file", {
  expect_identical(sidecar_path("out/res.tsv", "provenance"), "out/res_provenance.tsv")
  expect_identical(sidecar_path("out/res", "provenance"), "out/res_provenance.tsv")
  # The metadata sidecars keep their names
  expect_identical(sidecar_path("out/res.tsv", "metadata_x"), "out/res_metadata_x.tsv")
})

test_that("every value reads back with its type and exact value", {
  dir <- withr::local_tempdir()
  prov <- list(
    third = 1 / 3,                      # needs 17 significant digits
    tenth = 0.1,                        # 15 suffice, and it stays legible
    million_dbl = 1e6,
    million_int = 1000000L,             # must stay distinct from the double
    inf = Inf, neg_inf = -Inf, na_dbl = NA_real_, na_lgl = NA,
    flag = TRUE, none = NULL,
    cols_by_index = 1:4,
    cols_by_name = c("2", "gene name", "a\tb", "NA", ""),
    empty_chr = character(0), empty_int = integer(0)
  )
  path <- file.path(dir, "r.tsv")
  suppressMessages(write_provenance_sidecar(prov, path))
  back <- read_provenance_sidecar(sidecar_path(path, "provenance"))

  expect_identical(back, prov)
  expect_identical(typeof(back$million_int), "integer")
  expect_identical(typeof(back$million_dbl), "double")
  # Column 2 by index and a column named "2" stay different requests
  expect_identical(back$cols_by_index[2], 2L)
  expect_identical(back$cols_by_name[1], "2")

  sidecar <- readLines(sidecar_path(path, "provenance"))
  expect_true(any(grepl("^tenth\tdouble\t1\t0\\.1$", sidecar)))
})

test_that("the sidecar refuses values it could not read back faithfully", {
  dir <- withr::local_tempdir()
  path <- file.path(dir, "r.tsv")
  expect_error(write_provenance_sidecar(list(a = NA_character_), path), "missing string")
  expect_error(write_provenance_sidecar(list(a = list(1)), path), "unsupported type")
})

# --- Round trip through coreact_pipeline() and read_coreact_results() ---

test_that("a pipeline result reads back with its types and provenance", {
  dir <- withr::local_tempdir()
  out <- file.path(dir, "res.tsv")
  path_x <- system.file("extdata", "input_x.tsv", package = "coreact")
  path_y <- system.file("extdata", "input_y.tsv", package = "coreact")

  res <- suppressMessages(coreact_pipeline(
    paths = c(path_x, path_y), out_path = out,
    meta_cols = list(c("u_pep_id", "pep_nr", "protein_id", "gene_name"),
                     c("protein_id", "gene_name")),
    feature_ids = list("u_pep_id", c("protein_id", "gene_name")),
    min_prevalence = c(0.1, 2), fdr_threshold = 1
  ))
  back <- read_coreact_results(out)

  # Values, column types and both attributes. Doubles in the results file
  # carry 15 significant digits, so they agree to tolerance only.
  expect_equal(back, res)
  expect_identical(vapply(back, typeof, ""), vapply(res, typeof, ""))
  expect_identical(attr(back, "provenance"), attr(res, "provenance"))
  expect_identical(attr(back, "n_pairs"), attr(res, "n_pairs"))

  # Per-side arguments are stored per side, as resolved
  prov <- attr(back, "provenance")
  expect_identical(prov$x_meta_cols, c("u_pep_id", "pep_nr", "protein_id", "gene_name"))
  expect_identical(prov$y_feature_ids, c("protein_id", "gene_name"))
  expect_equal(prov$x_min_prevalence, 0.1)
  expect_equal(prov$y_min_prevalence, 2)
  expect_null(prov$x_max_prevalence)
  expect_equal(prov$fdr_threshold, 1)
  expect_identical(prov$x_path, normalizePath(path_x, winslash = "/"))
})

test_that("adjust_fdr() works on a result read back from disk", {
  withr::with_seed(8, {
    dir <- withr::local_tempdir()
    path_x <- write_pipeline_tsv(25, 40, prefix = "X", dir = dir)
    path_y <- write_pipeline_tsv(15, 40, prefix = "Y", dir = dir)
  })
  out <- file.path(dir, "res.tsv")
  res <- suppressMessages(coreact_pipeline(c(path_x, path_y), out, meta_cols = 1:2,
                                           fdr_threshold = 1))
  back <- read_coreact_results(out)

  expect_equal(adjust_fdr(back)$p_adj, res$p_adj)
  expect_equal(adjust_fdr(back, method = "holm")$p_adj,
               stats::p.adjust(res$p_val, "holm", n = attr(res, "n_pairs")))
})

test_that("row-number feature IDs come back as text that still indexes the inputs", {
  withr::with_seed(9, {
    dir <- withr::local_tempdir()
    path <- write_pipeline_tsv(10, 30, prefix = "F", dir = dir)
  })
  out <- file.path(dir, "within.tsv")
  # feature_ids = NULL: IDs are the row numbers "1", "2", ...
  res <- suppressMessages(coreact_pipeline(path, out, meta_cols = 1:2, fdr_threshold = 1))
  back <- read_coreact_results(out)

  expect_type(back$feature_x, "character")
  expect_identical(back$feature_x, res$feature_x)
  obj <- suppressMessages(import_coreact_tsv(path, meta_cols = 1:2))
  expect_identical(rownames(obj$mat)[back$idx_x], back$feature_x)
  expect_true(attr(back, "provenance")$within_library)
})

test_that("an empty results file reads back with the right column types", {
  dir <- withr::local_tempdir()
  df <- data.frame(ID = paste0("G", 1:4), matrix(0, 4, 5, dimnames = list(NULL, paste0("S", 1:5))))
  path <- file.path(dir, "zeros.tsv")
  utils::write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  out <- file.path(dir, "empty.tsv")

  res <- suppressWarnings(suppressMessages(
    coreact_pipeline(path, out, meta_cols = 1, effect_filters = list(min_intersection = 1))
  ))
  back <- read_coreact_results(out)

  expect_equal(nrow(back), 0)
  expect_identical(vapply(back, typeof, ""), vapply(res, typeof, ""))
  expect_equal(attr(back, "n_pairs"), 4 * 3 / 2)
})

test_that("a results file without a sidecar reads back with a warning", {
  withr::with_seed(10, {
    dir <- withr::local_tempdir()
    path <- write_pipeline_tsv(6, 20, prefix = "F", dir = dir)
  })
  out <- file.path(dir, "res.tsv")
  suppressMessages(coreact_pipeline(path, out, meta_cols = 1:2, fdr_threshold = 1))
  file.remove(sidecar_path(out, "provenance"))

  expect_warning(back <- read_coreact_results(out), "No provenance sidecar")
  expect_null(attr(back, "n_pairs"))
  expect_error(adjust_fdr(back), "n_pairs")
})

test_that("read_coreact_results() refuses files that are not coreact results", {
  dir <- withr::local_tempdir()
  other <- file.path(dir, "other.tsv")
  utils::write.table(data.frame(a = 1, b = 2), other, sep = "\t", row.names = FALSE)

  expect_error(read_coreact_results(other), "not a coreact results file")
  expect_error(read_coreact_results(file.path(dir, "absent.tsv")), "File not found")
  expect_error(read_coreact_results(c("a", "b")), "single file path")
})

# --- The content hash ---

test_that("the content hash is fixed across platforms and R versions", {
  # Pinned: CI computes it on Linux, macOS and Windows. A non-ASCII name checks
  # that names are hashed as UTF-8 everywhere.
  m <- Matrix::sparseMatrix(i = c(1, 2, 3, 3), j = c(1, 3, 2, 4), x = 1, dims = c(3, 4),
                            dimnames = list(c("A", "caf\u00e9", "PROT_1|GENE_A"), paste0("S", 1:4)))
  expect_identical(matrix_md5(m), "38f0017dcbaf4351c05c86211c2b1db8")
})

test_that("the content hash sees the data, not how it is stored", {
  m <- Matrix::sparseMatrix(i = c(1, 2, 3, 3), j = c(1, 3, 2, 4), x = 1, dims = c(3, 4),
                            dimnames = list(c("A", "caf\u00e9", "C"), paste0("S", 1:4)))
  h <- matrix_md5(m)

  # The same matrix stored differently: an explicit zero, a pattern matrix, a
  # name in another encoding
  stored_zero <- Matrix::sparseMatrix(i = c(1, 1, 2, 3, 3), j = c(2, 1, 3, 2, 4), x = c(0, 1, 1, 1, 1),
                                      dims = c(3, 4), dimnames = dimnames(m))
  latin1 <- m
  rownames(latin1)[2] <- iconv("caf\u00e9", "UTF-8", "latin1")
  expect_identical(matrix_md5(stored_zero), h)
  expect_identical(matrix_md5(methods::as(m, "nMatrix")), h)
  expect_identical(matrix_md5(latin1), h)

  # Any change to the data: a positive moved, a name changed, columns in
  # another order (nothing is reordered before hashing)
  moved <- m
  moved[1, 1] <- 0
  moved[1, 2] <- 1
  renamed <- m
  rownames(renamed)[1] <- "B"
  expect_false(identical(matrix_md5(moved), h))
  expect_false(identical(matrix_md5(renamed), h))
  expect_false(identical(matrix_md5(m[, c(2, 1, 3, 4)]), h))
  # Hashing leaves the matrix as it was
  expect_identical(m, Matrix::sparseMatrix(i = c(1, 2, 3, 3), j = c(1, 3, 2, 4), x = 1, dims = c(3, 4),
                                           dimnames = list(c("A", "caf\u00e9", "C"), paste0("S", 1:4))))
})

test_that("provenance records each side's content hash", {
  withr::with_seed(31, {
    o_x <- make_obj(5, 20, prefix = "X")
    o_y <- make_obj(4, 20, prefix = "Y")
  })
  prov <- attr(suppressMessages(run_coreact(o_x, o_y, effect_filters = list())), "provenance")
  expect_identical(prov$x_md5, matrix_md5(o_x$mat))
  expect_identical(prov$y_md5, matrix_md5(o_y$mat))
})
