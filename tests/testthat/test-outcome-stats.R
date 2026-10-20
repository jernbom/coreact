# Statistical behaviour of run_coreact_outcome().

# --- Classification never looks at the outcome ---

test_that("permuting the outcome leaves every classification unchanged", {
  # The single property whose violation would invalidate every adjusted
  # p-value: the choice of what to test must not depend on the outcome.
  d <- outcome_fixture(seed = 11)
  keep <- c("feature_x", "feature_y", "n00", "n10", "n01", "n11", "m_off", "n_used",
            "occupancy_class", "contrast", "role")
  for (type in c("binary", "continuous")) {
    o <- d[[type]]
    o[c(3, 17, 40, 41)] <- NA
    # Permute among the samples that have an outcome: moving a missing value
    # would change the analysis sample, which is not what is being tested
    shuffled <- o
    has <- !is.na(o)
    shuffled[has] <- withr::with_seed(12, sample(unname(o[has])))
    a <- suppressWarnings(suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = o,
                                                               outcome_type = type)))
    b <- suppressWarnings(suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = shuffled,
                                                               outcome_type = type)))
    expect_identical(a[keep], b[keep], info = type)
  }
})

# --- Occupancy classes ---

test_that("each class is assigned and tested as specified", {
  fixtures <- list(
    disjoint = list(n = cells4(20, 10, 10, 0), tested = character(0)),
    degenerate = list(n = cells4(0, 0, 0, 20), tested = character(0)),
    identical = list(n = cells4(20, 0, 0, 15), tested = "11_vs_00"),
    nested_y_in_x = list(n = cells4(20, 10, 0, 15), tested = "11_vs_10"),
    nested_x_in_y = list(n = cells4(20, 0, 10, 15), tested = "11_vs_01"),
    no_double_negative = list(n = cells4(0, 10, 12, 15), tested = c("11_vs_10", "11_vs_01")),
    near_nested = list(n = cells4(40, 15, 5, 12), tested = c("interaction", "11_vs_10")),
    full = list(n = cells4(40, 15, 12, 12), tested = "interaction")
  )
  for (cls in names(fixtures)) {
    n <- fixtures[[cls]]$n
    out <- run_fixture(pair_fixture(n, cases = floor(n / 2)))

    expect_true(all(out$occupancy_class == cls), info = cls)
    expect_identical(as.character(out$contrast[out$role == "tested"]), fixtures[[cls]]$tested,
                     info = cls)
    # Every contrast is reported; estimable exactly when its cells are occupied
    expect_identical(as.character(out$contrast), contrast_catalogue, info = cls)
    occupied <- vapply(contrast_weights, function(w) all(n[w != 0] > 0), logical(1))
    expect_identical(out$beta_estimability != "not_identified", unname(occupied), info = cls)
  }
})

test_that("near_nested tests the conditional contrast that avoids the thinner cell", {
  thin <- function(n) {
    out <- run_fixture(pair_fixture(n, cases = floor(n / 2)))
    as.character(out$contrast[out$role == "tested"])
  }
  expect_identical(thin(cells4(40, 15, 5, 12)), c("interaction", "11_vs_10"))
  expect_identical(thin(cells4(40, 5, 15, 12)), c("interaction", "11_vs_01"))
  # The one documented asymmetry: a tie goes to 11_vs_10, in either orientation
  expect_identical(thin(cells4(40, 6, 6, 12)), c("interaction", "11_vs_10"))
})

test_that("near_nested_max is a count that moves pairs between near_nested and full", {
  f <- pair_fixture(cells4(40, 15, 6, 12), cases = cells4(20, 7, 3, 6))
  expect_true(all(run_fixture(f, near_nested_max = 10)$occupancy_class == "near_nested"))
  expect_true(all(run_fixture(f, near_nested_max = 6)$occupancy_class == "full"))
  # A function is applied to n_used, the samples with an outcome, and the value
  # it returns is recorded
  f$outcome[1:3] <- NA
  out <- run_fixture(f, near_nested_max = function(n) n / 10)
  expect_equal(out$n_used[1], 70L)
  expect_equal(attr(out, "near_nested_max"), 70 / 10)
  expect_true(all(out$occupancy_class == "near_nested"))
})

test_that("narrowing contrasts never drops a tested one", {
  nested <- pair_fixture(cells4(20, 10, 0, 15), cases = cells4(10, 5, 0, 7))
  out <- run_fixture(nested, contrasts = "interaction")
  expect_identical(as.character(out$contrast), c("interaction", "11_vs_10"))
  expect_identical(as.character(out$role), c("reported", "tested"))

  full <- pair_fixture(cells4(40, 15, 12, 12), cases = cells4(20, 7, 6, 6))
  out <- run_fixture(full, contrasts = c("10_vs_00", "10_vs_00"))
  expect_identical(as.character(out$contrast), c("interaction", "10_vs_00"))
})

test_that("all sixteen occupancy patterns are tested only where estimable", {
  # Enumerating patterns, not classes, is what would have caught a class whose
  # tested contrast references an empty cell
  patterns <- expand.grid(rep(list(c(0L, 6L)), 4))
  for (i in seq_len(nrow(patterns))) {
    n <- stats::setNames(unlist(patterns[i, ]), cell_order)
    cls <- classify_occupancy(n, 10)
    for (cn in cls$tested) {
      w <- contrast_weights[[cn]]
      expect_true(all(n[w != 0] > 0), info = paste(c(names(n)[n > 0], cn), collapse = " "))
    }
    if (length(cls$tested) == 0) expect_true(cls$class %in% c("disjoint", "degenerate"))
    if (sum(n) == 0) next

    out <- run_fixture(pair_fixture(n, cases = n %/% 2L))
    expected <- if (cls$class %in% one_test_classes) 1 else if (cls$class %in% two_test_classes) 2 else 0
    expect_equal(sum(out$role == "tested" & !is.na(out$p_val)), expected,
                 info = paste(names(n)[n > 0], collapse = " "))
  }
})

# --- Estimates against closed forms ---

test_that("binary estimates and SEs match the closed forms of the saturated model", {
  n <- cells4(200, 60, 60, 20)
  cases <- cells4(40, 18, 25, 12)
  out <- run_fixture(pair_fixture(n, cases = cases))
  logit <- log(cases / (n - cases))
  woolf <- 1 / cases + 1 / (n - cases)
  for (cn in contrast_catalogue) {
    w <- contrast_weights[[cn]]
    row <- out[out$contrast == cn, ]
    expect_equal(row$beta, sum(w * logit), tolerance = 1e-8, info = cn)
    expect_equal(row$se, sqrt(sum(w^2 * woolf)), tolerance = 1e-6, info = cn)
  }
  # The interaction is the log ratio of the eight outcome-stratified counts
  i <- out[out$contrast == "interaction", ]
  expect_equal(i$beta, log((12 * 42 * 35 * 40) / (8 * 18 * 25 * 160)), tolerance = 1e-8)
})

test_that("continuous estimates and SEs match the closed forms", {
  n <- cells4(120, 40, 30, 15)
  o <- withr::with_seed(5, stats::rnorm(sum(n), rep(c(0, 0.4, 0.3, 1.2), n)))
  out <- run_fixture(pair_fixture(n, outcome = o), outcome_type = "continuous")
  cell <- rep(1:4, n)
  means <- tapply(o, cell, mean)
  sigma <- sqrt(sum((o - means[cell])^2) / (sum(n) - 4))
  for (cn in contrast_catalogue) {
    w <- contrast_weights[[cn]]
    row <- out[out$contrast == cn, ]
    expect_equal(row$beta, sum(w * means), tolerance = 1e-10, info = cn)
    expect_equal(row$se, sigma * sqrt(sum(w^2 / n)), tolerance = 1e-10, info = cn)
  }
})

test_that("a continuous outcome constant within every cell returns the limit p-values", {
  # Decimal values that are not exactly representable, so lm() leaves rounding
  # residue instead of an exact zero. Cells 10 and 11 hold the same value.
  n <- cells4(40, 15, 5, 12)
  o <- rep(c(0.1, 0.3, 0.7, 0.3), n)
  out <- run_fixture(pair_fixture(n, outcome = o), outcome_type = "continuous")
  tested <- out[out$role == "tested", ]
  expect_identical(as.character(tested$contrast), c("interaction", "11_vs_10"))

  # Merging changes the fit for the interaction (0.3 - 0.3 - 0.7 + 0.1 != 0)...
  expect_identical(tested$p_val[tested$contrast == "interaction"], 0)
  # ...and not for 11_vs_10, whose cells are equal. Read without a tolerance,
  # the rounding residue would make this highly significant.
  expect_identical(tested$p_val[tested$contrast == "11_vs_10"], 1)
  # The run completed, so the family-size guard held
  expect_equal(sum(!is.na(out$p_adj)), 2)
  # With no spread the estimates are exact differences of cell values and the
  # standard errors are zero, not rounding residue
  expect_equal(out$beta[out$contrast == "interaction"], -0.6, tolerance = 1e-12)
  expect_true(all(out$se[out$beta_estimability == "estimable"] == 0))
})

test_that("the continuous F test equals the squared t statistic of its contrast", {
  n <- cells4(120, 40, 30, 15)
  o <- withr::with_seed(6, stats::rnorm(sum(n), rep(c(0, 0.4, 0.3, 1.2), n)))
  out <- run_fixture(pair_fixture(n, outcome = o), outcome_type = "continuous")
  i <- out[out$contrast == "interaction", ]
  expect_equal(i$statistic, (i$beta / i$se)^2, tolerance = 1e-10)
  expect_equal(i$p_val, 2 * stats::pt(-abs(i$beta / i$se), sum(n) - 4), tolerance = 1e-10)
})

# --- Separation ---

test_that("the four separation scenarios are reported as specified", {
  n <- cells4(200, 60, 60, 20)
  run <- function(cases) {
    out <- run_fixture(pair_fixture(n, cases = cases))
    out[out$contrast == "interaction", ]
  }

  # Cell 11 all cases: zero only on the denominator side
  a <- run(cells4(40, 20, 20, 20))
  expect_identical(as.character(a$beta_estimability), "unbounded_positive")
  expect_identical(a$beta, Inf)
  expect_identical(a$se, Inf)
  expect_lt(a$p_val, 1e-3)                 # the likelihood ratio survives

  # Cells 10 and 11 all cases: zeros on both sides, so no limit exists
  b <- run(cells4(40, 60, 20, 20))
  expect_identical(as.character(b$beta_estimability), "indeterminate")
  expect_true(is.na(b$beta))
  expect_false(is.na(b$p_val))

  # Cell 10 all controls and cell 11 all cases: both zeros on the denominator side
  both_den <- run(cells4(40, 0, 20, 20))
  expect_identical(as.character(both_den$beta_estimability), "unbounded_positive")
  expect_identical(both_den$beta, Inf)
  expect_lt(both_den$p_val, 1e-6)

  # No zero cell
  none <- run(cells4(40, 20, 20, 12))
  expect_identical(as.character(none$beta_estimability), "estimable")
  expect_true(is.finite(none$beta))
})

test_that("a separated fit leaves contrasts away from the separated cell estimable", {
  out <- run_fixture(pair_fixture(cells4(200, 60, 60, 20), cases = cells4(40, 20, 20, 20)))
  away <- out[out$contrast %in% c("10_vs_00", "01_vs_00", "10_vs_01"), ]
  expect_true(all(away$beta_estimability == "estimable"))
  expect_true(all(is.finite(away$beta)))
  # The negative direction: cell 11 all controls
  neg <- run_fixture(pair_fixture(cells4(200, 60, 60, 20), cases = cells4(40, 20, 20, 0)))
  expect_identical(as.character(neg$beta_estimability[neg$contrast == "interaction"]),
                   "unbounded_negative")
  expect_identical(neg$beta[neg$contrast == "interaction"], -Inf)
})

test_that("separation of a conditional contrast takes its direction from its own cells", {
  n <- cells4(200, 60, 0, 20)            # nested_y_in_x, tested on 11_vs_10
  up <- run_fixture(pair_fixture(n, cases = cells4(40, 20, 0, 20)))
  down <- run_fixture(pair_fixture(n, cases = cells4(40, 20, 0, 0)))
  t_up <- up[up$role == "tested", ]
  t_down <- down[down$role == "tested", ]
  expect_identical(as.character(t_up$contrast), "11_vs_10")
  expect_identical(as.character(t_up$beta_estimability), "unbounded_positive")
  expect_identical(as.character(t_down$beta_estimability), "unbounded_negative")
  expect_lt(t_up$p_val, 1e-3)
  expect_lt(t_down$p_val, 1e-3)
})

test_that("an empty predictor cell is aliasing, never separation", {
  # A nested pair (cell 01 empty) whose occupied cell 11 is all cases
  out <- run_fixture(pair_fixture(cells4(200, 60, 0, 20), cases = cells4(40, 20, 0, 20)))
  # The empty cell has no cases and no controls, which must not read as zeros
  # on both sides of the interaction
  expect_identical(as.character(out$beta_estimability[out$contrast == "interaction"]), "not_identified")
  # Only the occupied cell's zero is counted: 11 has no controls
  expect_true(all(out$zero_cells_8 == 1L))

  clean <- run_fixture(pair_fixture(cells4(200, 60, 60, 20), cases = cells4(40, 20, 20, 12)))
  expect_true(all(clean$zero_cells_8 == 0L))
})

test_that("separation is decided from the counts, not by a failing detector", {
  build <- function(n, cases) {
    cell <- rep(1:4, n)
    list(cell = cell,
         o = unlist(lapply(1:4, function(k) rep(c(1L, 0L), c(cases[k], n[k] - cases[k])))))
  }
  # Cells 10 and 11 all cases: glm reports convergence and an innocent-looking
  # interaction, although the fit is separated
  ind <- build(c(200, 60, 60, 20), c(40, 60, 20, 20))
  fit <- fit_cell_means(ind$cell, ind$o, "binary", cell_counts(ind$cell))
  naive <- contrast_estimate(fit$model, contrast_weights$interaction, cell_counts(ind$cell))
  expect_true(fit$model$converged)            # "converged" does not rule separation out
  expect_lt(abs(naive$beta), 10)              # nor does a bound on |beta|

  # One separated cell of 3 samples at N = 2998: glm stops well short of the
  # bound, so even a fitted-probability threshold of 1e-6 misses it
  one <- build(c(2000, 500, 495, 3), c(600, 150, 148, 3))
  fit_one <- fit_cell_means(one$cell, one$o, "binary", cell_counts(one$cell))
  expect_true(fit_one$model$converged)
  expect_gt(1 - max(stats::fitted(fit_one$model)), 1e-6)

  # The count rule gets both right
  x_of <- function(cell) as.integer(cell %in% c(2, 4))
  y_of <- function(cell) as.integer(cell %in% c(3, 4))
  r_ind <- analyse_pair(x_of(ind$cell), y_of(ind$cell), ind$o, "binary", 10, "interaction")
  r_one <- analyse_pair(x_of(one$cell), y_of(one$cell), one$o, "binary", 10, "interaction")
  expect_identical(r_ind$beta_estimability[r_ind$contrast == "interaction"], "indeterminate")
  expect_identical(r_one$beta_estimability[r_one$contrast == "interaction"], "unbounded_positive")
  expect_identical(r_one$beta[r_one$contrast == "interaction"], Inf)
})

# --- Size under the null ---

test_that("tests hold their size under an additive truth", {
  n <- c(330, 100, 50, 20)
  cell <- rep(1:4, n)
  x <- cell %in% c(2, 4)
  y <- cell %in% c(3, 4)
  interaction_p <- function(o, type) {
    r <- analyse_pair(x, y, o, type, 10, "interaction")
    r$p_val[r$contrast == "interaction"]
  }
  withr::local_seed(2024)
  p_bin <- replicate(400, interaction_p(stats::rbinom(500, 1, stats::plogis(-1 + 0.5 * x + 0.4 * y)),
                                        "binary"))
  p_con <- replicate(400, interaction_p(stats::rnorm(500, 0.5 * x + 0.4 * y), "continuous"))
  # 400 replicates: the standard error of a 5% rate is about 1.1%
  for (p in list(p_bin, p_con)) {
    expect_gt(mean(p < 0.05), 0.02)
    expect_lt(mean(p < 0.05), 0.085)
  }
})

# --- The analysis sample ---

test_that("a cell whose outcomes are all missing is treated as empty", {
  n <- cells4(40, 15, 6, 12)
  f <- pair_fixture(n, cases = floor(n / 2))
  expect_true(all(run_fixture(f)$occupancy_class == "near_nested"))

  in_01 <- as.numeric(f$x$mat[1, ]) == 0 & as.numeric(f$y$mat[1, ]) == 1
  f$outcome[in_01] <- NA
  out <- run_fixture(f)                    # without complete-case filtering, this errors

  expect_true(all(out$n01 == 0))
  expect_true(all(out$n_used == sum(n) - 6))
  expect_true(all(out$occupancy_class == "nested_y_in_x"))
  expect_setequal(as.character(out$contrast[out$beta_estimability == "not_identified"]),
                  c("interaction", "11_vs_01", "01_vs_00", "10_vs_01"))
})

# --- Swapping x and y ---

test_that("swapping x and y exchanges the directional contrasts and nothing else", {
  d <- outcome_fixture(seed = 21)
  res_ba <- suppressMessages(run_coreact(d$y, d$x, effect_filters = list()))
  for (type in c("binary", "continuous")) {
    o <- d[[type]]
    ab <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = o, outcome_type = type))
    ba <- suppressMessages(run_coreact_outcome(res_ba, d$y, d$x, outcome = o, outcome_type = type))

    # Ties in the thin cell are the one documented asymmetry; leave them out
    ties <- unique(paste(ab$feature_x, ab$feature_y)[ab$occupancy_class == "near_nested" & ab$n01 == ab$n10])
    ab <- ab[!paste(ab$feature_x, ab$feature_y) %in% ties, ]
    ba <- ba[!paste(ba$feature_y, ba$feature_x) %in% ties, ]

    swap_contrast <- c("interaction" = "interaction", "11_vs_10" = "11_vs_01", "11_vs_01" = "11_vs_10",
                       "11_vs_00" = "11_vs_00", "10_vs_00" = "01_vs_00", "01_vs_00" = "10_vs_00",
                       "10_vs_01" = "10_vs_01")
    swap_class <- c(nested_y_in_x = "nested_x_in_y", nested_x_in_y = "nested_y_in_x")
    key_ab <- paste(ab$feature_x, ab$feature_y, ab$contrast)
    key_ba <- paste(ba$feature_y, ba$feature_x, swap_contrast[as.character(ba$contrast)])
    ba <- ba[match(key_ab, key_ba), ]

    cls_ba <- as.character(ba$occupancy_class)
    cls_ba[cls_ba %in% names(swap_class)] <- swap_class[cls_ba[cls_ba %in% names(swap_class)]]
    expect_identical(as.character(ab$occupancy_class), cls_ba, info = type)
    expect_identical(as.character(ab$role), as.character(ba$role), info = type)
    expect_identical(as.character(ab$beta_estimability), as.character(ba$beta_estimability), info = type)

    sign <- ifelse(ab$contrast == "10_vs_01", -1, 1)
    expect_equal(ab$beta, sign * ba$beta, tolerance = 1e-7, info = type)
    expect_equal(ab$se, ba$se, tolerance = 1e-7, info = type)
    expect_equal(ab$p_val, ba$p_val, tolerance = 1e-7, info = type)
    expect_equal(ab$p_adj, ba$p_adj, tolerance = 1e-7, info = type)
  }
})

# --- Against a brute-force reference ---

test_that("the engine matches a treatment-coded reference fit for every pair", {
  d <- outcome_fixture(seed = 31)
  for (type in c("binary", "continuous")) {
    o <- d[[type]]
    out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = o, outcome_type = type))
    pair <- pair_index(out)
    for (p in unique(pair)) {
      rows <- out[pair == p, ]
      x <- as.numeric(d$x$mat[rows$feature_x[1], ])
      y <- as.numeric(d$y$mat[rows$feature_y[1], ])
      cls <- as.character(rows$occupancy_class[1])
      info <- paste(type, rows$feature_x[1], rows$feature_y[1], cls)

      ref <- reference_tests(x, y, o, type, cls)
      tested <- rows[rows$role == "tested", ]
      expect_setequal(as.character(tested$contrast), as.character(names(ref)))
      for (cn in names(ref)) {
        expect_equal(tested$p_val[tested$contrast == cn], ref[[cn]], tolerance = 1e-5, info = info)
      }

      for (i in seq_len(nrow(rows))) {
        want <- closed_form_beta(contrast_weights[[as.character(rows$contrast[i])]], x, y, o, type)
        if (is.na(want)) {
          # An empty cell or an all-case / all-control cell: no finite estimate
          expect_false(rows$beta_estimability[i] == "estimable", info = paste(info, rows$contrast[i]))
        } else {
          expect_equal(rows$beta[i], want, tolerance = 1e-6, info = paste(info, rows$contrast[i]))
        }
      }
    }
  }
})

test_that("results do not depend on chunk size or on running in parallel", {
  d <- outcome_fixture(seed = 41)
  run <- function(...) {
    suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary, ...))
  }
  ref <- run(chunk_size = 1)
  expect_identical(run(chunk_size = 7), ref)
  expect_identical(run(), ref)

  skip_on_os("windows")
  expect_identical(run(n_cores = 2, chunk_size = 5), ref)
})

test_that("on Windows, more than one core falls back to serial with a warning", {
  skip_on_os(c("mac", "linux", "solaris"))
  f <- pair_fixture(cells4(40, 15, 12, 12), cases = cells4(20, 7, 6, 6))
  expect_warning(run_fixture(f, n_cores = 2), "Running serially")
})
