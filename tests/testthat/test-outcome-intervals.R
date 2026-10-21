# Confidence intervals of run_coreact_outcome().

# Independent references: for a binary outcome, the profile interval of a
# treatment-coded glm fitted to the contrast's own cells (the likelihood
# factorises over cells, so the other cells do not matter); for the interaction,
# the x:y term of glm(o ~ x * y). For a continuous outcome, confint() on an lm
# with the contrast as a coefficient, fitted to all cells (the pooled residual
# standard deviation uses them all).
reference_interval <- function(contrast, x, y, o, type) {
  cell <- c("00", "10", "01", "11")[1 + x + 2 * y]
  quiet <- function(expr) suppressWarnings(suppressMessages(expr))
  if (contrast == "interaction") {
    d <- data.frame(o = o, x = x, y = y)
    m <- if (type == "binary") quiet(stats::glm(o ~ x * y, family = stats::binomial(), data = d))
         else stats::lm(o ~ x * y, data = d)
    # confint() gives a vector for a glm coefficient and a 1 x 2 matrix for an lm one
    return(as.numeric(quiet(stats::confint(m, "x:y"))))
  }
  ends <- strsplit(contrast, "_vs_", fixed = TRUE)[[1]]
  if (type == "binary") {
    keep <- cell %in% ends
    d <- data.frame(o = o[keep], a = as.integer(cell[keep] == ends[1]))
    m <- quiet(stats::glm(o ~ a, family = stats::binomial(), data = d))
    return(as.numeric(quiet(stats::confint(m, "a"))))
  }
  d <- data.frame(o = o, g = stats::relevel(factor(cell), ref = ends[2]))
  as.numeric(stats::confint(stats::lm(o ~ g, data = d))[paste0("g", ends[1]), ])
}

test_that("ci chooses which estimable rows get an interval", {
  d <- outcome_fixture(seed = 51)
  run <- function(ci) {
    suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary, ci = ci))
  }
  estimable <- function(out) out$beta_estimability == "estimable"

  tested <- run("tested")
  has <- !is.na(tested$ci_lo)
  expect_identical(has, estimable(tested) & tested$role == "tested")

  all_rows <- run("all")
  expect_identical(!is.na(all_rows$ci_lo), estimable(all_rows))
  # Within the estimable rows, intervals bracket the estimate
  e <- estimable(all_rows)
  expect_true(all(all_rows$ci_lo[e] < all_rows$beta[e] & all_rows$beta[e] < all_rows$ci_hi[e]))

  none <- run("none")
  expect_true(all(is.na(none$ci_lo)) && all(is.na(none$ci_hi)))
  # Intervals change nothing else in the table; only the recorded setting differs
  same <- setdiff(names(none), c("ci_lo", "ci_hi"))
  without_ci <- function(t) {
    attr(t, "ci") <- NULL
    t[same]
  }
  expect_identical(without_ci(none), without_ci(all_rows))
})

test_that("binary intervals are profile intervals, matching an independent fit", {
  d <- outcome_fixture(seed = 52)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$binary, ci = "all"))
  rows <- which(out$beta_estimability == "estimable")
  for (i in rows) {
    x <- as.numeric(d$x$mat[out$feature_x[i], ])
    y <- as.numeric(d$y$mat[out$feature_y[i], ])
    want <- reference_interval(as.character(out$contrast[i]), x, y, d$binary, "binary")
    expect_equal(c(out$ci_lo[i], out$ci_hi[i]), want, tolerance = 1e-5,
                 info = paste(out$feature_x[i], out$feature_y[i], out$contrast[i]))
  }
  # Not Wald: a profile interval is not symmetric about the estimate
  e <- out[rows, ]
  expect_false(isTRUE(all.equal(e$beta - e$ci_lo, e$ci_hi - e$beta)))
})

test_that("continuous intervals are exact t intervals on the pooled fit", {
  d <- outcome_fixture(seed = 53)
  out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d$continuous,
                                              outcome_type = "continuous", ci = "all"))
  rows <- which(out$beta_estimability == "estimable")
  for (i in rows) {
    x <- as.numeric(d$x$mat[out$feature_x[i], ])
    y <- as.numeric(d$y$mat[out$feature_y[i], ])
    want <- reference_interval(as.character(out$contrast[i]), x, y, d$continuous, "continuous")
    expect_equal(c(out$ci_lo[i], out$ci_hi[i]), want, tolerance = 1e-10,
                 info = paste(out$feature_x[i], out$feature_y[i], out$contrast[i]))
  }
})

test_that("an interval excludes zero exactly when its test is significant", {
  # The interval inverts the test: profile likelihood against the likelihood
  # ratio test, the t interval against the F test
  for (type in c("binary", "continuous")) {
    found <- 0
    for (seed in 61:63) {
      d <- outcome_fixture(n_samples = 220, seed = seed)
      out <- suppressMessages(run_coreact_outcome(d$res, d$x, d$y, outcome = d[[type]],
                                                  outcome_type = type))
      t <- out[!is.na(out$ci_lo), ]
      # The profile interval is interpolated, so leave out p-values at the edge
      t <- t[abs(t$p_val - 0.05) > 0.002, ]
      excludes <- t$ci_lo > 0 | t$ci_hi < 0
      expect_identical(excludes, t$p_val < 0.05, info = paste(type, seed))
      found <- found + sum(excludes)
    }
    expect_gt(found, 0)                    # both kinds of row occurred
  }
})

test_that("separated rows have no interval, and their neighbours keep theirs", {
  f <- pair_fixture(cells4(200, 60, 60, 20), cases = cells4(40, 20, 20, 20))
  out <- run_fixture(f, ci = "all")
  separated <- out$beta_estimability != "estimable"
  expect_true(all(c("interaction", "11_vs_10") %in% as.character(out$contrast[separated])))
  expect_true(all(is.na(out$ci_lo[separated]) & is.na(out$ci_hi[separated])))

  # Contrasts away from cell 11 profile normally, although the refit includes
  # the diverging cell
  x <- as.numeric(f$x$mat[1, ])
  y <- as.numeric(f$y$mat[1, ])
  for (cn in c("10_vs_00", "01_vs_00", "10_vs_01")) {
    row <- out[out$contrast == cn, ]
    expect_equal(c(row$ci_lo, row$ci_hi), reference_interval(cn, x, y, f$outcome, "binary"),
                 tolerance = 1e-5, info = cn)
  }
})

test_that("with no spread in any cell, the interval is the point estimate", {
  n <- cells4(40, 15, 5, 12)
  out <- run_fixture(pair_fixture(n, outcome = rep(c(0.1, 0.3, 0.7, 0.3), n)),
                     outcome_type = "continuous", ci = "all")
  e <- out[out$beta_estimability == "estimable", ]
  expect_identical(e$ci_lo, e$beta)
  expect_identical(e$ci_hi, e$beta)
})
