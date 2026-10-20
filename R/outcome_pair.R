# Per-pair statistics for run_coreact_outcome().
#
# Each pair's two binary features cross into a 2x2 predictor table with cells
# 00, 10, 01, 11 (x then y). The outcome is modelled by one mean per OCCUPIED
# cell, and every quantity reported is a linear contrast of those cell means.
# Which contrasts are estimable, and which are tested, follows from which cells
# are occupied -- never from the outcome.

# --- Constants ---------------------------------------------------------------

# Cell order used throughout. With x and y coded 0/1, 1 + x + 2 * y indexes it.
cell_order <- c("00", "10", "01", "11")

# Weights on (mu00, mu10, mu01, mu11). The names are deliberately not syntactic,
# so that a wide pivot yields beta_11_vs_10 rather than beta_beta_11_vs_10.
contrast_weights <- list(
  "interaction" = c(1, -1, -1, 1),
  "11_vs_10" = c(0, -1, 0, 1),
  "11_vs_01" = c(0, 0, -1, 1),
  "11_vs_00" = c(-1, 0, 0, 1),
  "10_vs_00" = c(-1, 1, 0, 0),
  "01_vs_00" = c(-1, 0, 1, 0),
  "10_vs_01" = c(0, 1, -1, 0)
)
contrast_catalogue <- names(contrast_weights)

# Every contrast but the interaction is zero exactly when its two cells share a
# mean, so its reduced model merges those two levels. The interaction's reduced
# model is the additive one.
contrast_merge <- list(
  "11_vs_10" = c("10", "11"),
  "11_vs_01" = c("01", "11"),
  "11_vs_00" = c("00", "11"),
  "10_vs_00" = c("00", "10"),
  "01_vs_00" = c("00", "01"),
  "10_vs_01" = c("10", "01")
)

estimability_levels <- c("estimable", "unbounded_positive", "unbounded_negative",
                         "indeterminate", "not_identified")

occupancy_levels <- c("disjoint", "degenerate", "identical", "nested_y_in_x",
                      "nested_x_in_y", "no_double_negative", "near_nested", "full")

# Classes contributing one and two rows to the FDR family; the rest contribute
# none. Used by the family-size guard.
one_test_classes <- c("identical", "nested_y_in_x", "nested_x_in_y", "full")
two_test_classes <- c("near_nested", "no_double_negative")

# A residual sum of squares at or below this fraction of the analysis sample's
# total sum of squares is zero: lm() leaves rounding residue rather than an
# exact zero when every cell is constant.
zero_rss_rel_tol <- .Machine$double.eps


# --- Counting ----------------------------------------------------------------

#' Predictor Cell Counts of a Pair (Internal)
#'
#' @param cell Integer vector, `1 + x + 2 * y`, indexing `cell_order`.
#' @return Named integer vector `c("00", "10", "01", "11")`.
#' @keywords internal
cell_counts <- function(cell) {
  stats::setNames(tabulate(cell, nbins = 4L), cell_order)
}

#' Outcome-Stratified Cell Counts of a Pair (Internal)
#'
#' @param cell As for [cell_counts()].
#' @param outcome Integer vector of 0 / 1.
#' @return 4 x 2 integer matrix: rows `cell_order`, columns outcome `"0"`, `"1"`.
#' @keywords internal
cell_counts_8 <- function(cell, outcome) {
  matrix(tabulate(cell + 4L * outcome, nbins = 8L), nrow = 4L,
         dimnames = list(cell_order, c("0", "1")))
}


# --- Classification ----------------------------------------------------------

#' Classify a Pair by the Occupancy of Its 2x2 Table (Internal)
#'
#' @description
#' Decides which contrasts a pair is tested on, from its four predictor counts
#' alone. Never from the outcome: that is what makes choosing the estimand a
#' form of independent filtering, which costs no multiple-testing budget.
#'
#' @details
#' The order of the checks matters. `degenerate` (fewer than two occupied
#' cells) must come before `identical`, whose condition does not require
#' cell 00 to be occupied: otherwise a pair whose every sample is
#' double-positive would be tested on `11_vs_00`, a contrast it cannot
#' estimate, and the family-size guard would stop the whole run.
#'
#' `near_nested` tests, as well as the interaction, the conditional contrast
#' that avoids the thinner off-diagonal cell; a tie goes to `11_vs_10`.
#'
#' @param n Named cell counts from [cell_counts()].
#' @param near_nested_max Count below which the thinner off-diagonal cell makes
#'   a fully occupied pair `near_nested` rather than `full`.
#' @return List with `class`, `tested` (character) and `m_off`.
#' @keywords internal
classify_occupancy <- function(n, near_nested_max) {
  m_off <- min(n[["01"]], n[["10"]])
  out <- function(class, tested) list(class = class, tested = tested, m_off = m_off)

  if (n[["11"]] == 0) return(out("disjoint", character(0)))
  if (sum(n > 0) < 2) return(out("degenerate", character(0)))
  if (n[["01"]] == 0 && n[["10"]] == 0) return(out("identical", "11_vs_00"))
  if (n[["01"]] == 0) return(out("nested_y_in_x", "11_vs_10"))
  if (n[["10"]] == 0) return(out("nested_x_in_y", "11_vs_01"))
  # Three occupied cells: no interaction, and neither conditional contrast is
  # privileged, so both are tested
  if (n[["00"]] == 0) return(out("no_double_negative", c("11_vs_10", "11_vs_01")))

  avoids_thin <- if (n[["01"]] <= n[["10"]]) "11_vs_10" else "11_vs_01"
  if (m_off < near_nested_max) return(out("near_nested", c("interaction", avoids_thin)))
  out("full", "interaction")
}


# --- Fitting -----------------------------------------------------------------

#' Fit One Mean per Occupied Cell (Internal)
#'
#' @param cell As for [cell_counts()].
#' @param outcome Outcome vector on the analysis sample.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @param n Cell counts.
#' @return List with `model` and `data`, a self-contained model frame.
#' @keywords internal
fit_cell_means <- function(cell, outcome, outcome_type, n) {
  data <- data.frame(outcome = outcome, g = factor(cell_order[cell], levels = cell_order[n > 0]))
  list(model = fit_level_means(data, outcome_type), data = data)
}

#' Fit One Mean per Level of `g` (Internal)
#'
#' A factor with a single level cannot enter a model formula, so a single
#' group gets an intercept-only model, which is the same single mean. That
#' happens for a pair with one occupied cell, and for a reduced model that
#' merges the only two occupied cells.
#'
#' @param data Model frame with columns `outcome` and `g`.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @return The fitted model.
#' @keywords internal
fit_level_means <- function(data, outcome_type) {
  if (nlevels(data$g) == 1) return(fit_model(outcome ~ 1, data, outcome_type))
  fit_model(outcome ~ 0 + g, data, outcome_type)
}

#' Fit a Binary or Continuous Model (Internal)
#'
#' Logistic `glm()` for a binary outcome, `lm()` for a continuous one. `lm()`
#' rather than a Gaussian `glm()` because its intervals are the exact t
#' intervals that match the F test.
#'
#' @param formula,data Passed on.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @return The fitted model.
#' @keywords internal
fit_model <- function(formula, data, outcome_type) {
  if (outcome_type == "binary") {
    # Separation makes glm() warn about fitted probabilities of 0 or 1; it is
    # decided from the counts instead, see contrast_determinacy()
    suppressWarnings(stats::glm(formula, family = stats::binomial(), data = data))
  } else {
    stats::lm(formula, data = data)
  }
}


# --- Estimability, estimates and tests ---------------------------------------

#' Direction of a Binary Contrast, From the Counts (Internal)
#'
#' @description
#' Decides whether a contrast has a finite estimate, and if not which way it
#' diverges, from the outcome-by-cell counts. This is exact in the unadjusted
#' model: the cell-means fit is saturated, so each fitted probability is its
#' cell's observed proportion, and a cell is separated exactly when it has no
#' cases or no controls.
#'
#' @details
#' On the logit scale a contrast is `sum_k w_k * (log n_k1 - log n_k0)`. A cell
#' with weight +1 puts its cases in the numerator and its controls in the
#' denominator, a cell with weight -1 the reverse. Zeros only in the
#' denominator send the contrast to `+Inf`, only in the numerator to `-Inf`,
#' and on both sides leave it without a limit.
#'
#' Separation is not detected from the fitted probabilities. `glm()` stops at
#' a distance from 0 or 1 that grows with the sample size, so any fixed
#' threshold misses single-cell separation; `fit$converged` is `TRUE` on
#' separated fits; and a bound on `|beta|` misses fits whose diverging cells
#' cancel. Covariate adjustment will need a different check.
#'
#' Called only for identified contrasts, whose every cell is occupied. That is
#' load-bearing: an empty predictor cell has no cases and no controls, and
#' would turn aliasing into a spurious `indeterminate`.
#'
#' @param w Contrast weights.
#' @param counts8 Matrix from [cell_counts_8()].
#' @return One of `estimable`, `unbounded_positive`, `unbounded_negative`,
#'   `indeterminate`.
#' @keywords internal
contrast_determinacy <- function(w, counts8) {
  pos <- which(w > 0)
  neg <- which(w < 0)
  num <- c(counts8[pos, "1"], counts8[neg, "0"])
  den <- c(counts8[pos, "0"], counts8[neg, "1"])
  zero_num <- any(num == 0)
  zero_den <- any(den == 0)
  if (zero_num && zero_den) "indeterminate"
  else if (zero_den) "unbounded_positive"
  else if (zero_num) "unbounded_negative"
  else "estimable"
}

#' Estimate and Standard Error of a Contrast (Internal)
#'
#' @param model The cell-means fit.
#' @param w Contrast weights.
#' @param n Cell counts; only occupied cells have coefficients.
#' @param zero_se Logical. For a continuous fit with no spread in any cell
#'   (see [has_zero_spread()]), the standard error is zero; it is set so rather
#'   than computed from rounding residue.
#' @return List with `beta` and `se`.
#' @keywords internal
contrast_estimate <- function(model, w, n, zero_se = FALSE) {
  occupied <- n > 0
  coefs <- paste0("g", cell_order[occupied])
  w <- w[occupied]
  beta <- sum(w * stats::coef(model)[coefs])
  if (zero_se) return(list(beta = beta, se = 0))
  v <- stats::vcov(model)[coefs, coefs, drop = FALSE]
  list(beta = beta, se = sqrt(drop(crossprod(w, v %*% w))))
}

#' Whether a Continuous Fit Has No Spread in Any Cell (Internal)
#'
#' True when the residual sum of squares is at most `.Machine$double.eps` times
#' the total sum of squares of the analysis sample. `lm()` leaves rounding
#' residue rather than an exact zero when every cell is constant, and that
#' residue must not be read as variance.
#'
#' @param fit From [fit_cell_means()], continuous outcome.
#' @return Logical.
#' @keywords internal
has_zero_spread <- function(fit) {
  y <- fit$data$outcome
  stats::deviance(fit$model) <= zero_rss_rel_tol * sum((y - mean(y))^2)
}

#' Test One Contrast Against Its Reduced Model (Internal)
#'
#' A 1-df likelihood ratio test for a binary outcome, the exact F test for a
#' continuous one, never a Wald test. Both are two-sided: the statistics are
#' squared quantities, blind to the contrast's sign. The reduced model is the
#' cell-means fit with the contrast constrained to zero: two cells merged, or
#' for the interaction the additive model. Both are built from the same model
#' frame, with x and y reconstructed from the cell labels, so they cannot drift
#' apart. The statistics are computed from `deviance()` (the residual deviance
#' of a `glm`, the residual sum of squares of an `lm`) rather than through
#' `anova()`.
#'
#' @details
#' On the continuous path the F test assumes the four cells are equally
#' variable. If every occupied cell is constant, there is no residual variance
#' to scale the difference by: the test then returns `p = 0` when merging the
#' cells changes the fit and `p = 1` when it does not, the limits of the F test
#' as the residual variance goes to zero. "Constant" allows for rounding: a
#' residual sum of squares within `.Machine$double.eps` of the total sum of
#' squares counts as zero. Without that, `lm()`'s residue on decimal values is
#' read as signal, and two cells holding the same value come out significant.
#'
#' @param fit From [fit_cell_means()].
#' @param contrast Contrast name.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @return List with `statistic` and `p`.
#' @keywords internal
contrast_test <- function(fit, contrast, outcome_type) {
  reduced <- fit$data
  if (contrast == "interaction") {
    label <- as.character(reduced$g)
    reduced$x_ <- as.integer(substr(label, 1, 1))
    reduced$y_ <- as.integer(substr(label, 2, 2))
    model0 <- fit_model(outcome ~ x_ + y_, reduced, outcome_type)
  } else {
    merge <- contrast_merge[[contrast]]
    label <- as.character(reduced$g)
    label[label %in% merge] <- paste(merge, collapse = "+")
    reduced$g <- factor(label)
    model0 <- fit_level_means(reduced, outcome_type)
  }

  if (outcome_type == "binary") {
    # Nested maximum-likelihood fits: the difference is never negative except
    # by rounding
    stat <- max(stats::deviance(model0) - stats::deviance(fit$model), 0)
    return(list(statistic = stat, p = stats::pchisq(stat, 1, lower.tail = FALSE)))
  }

  rss1 <- stats::deviance(fit$model)
  rss0 <- stats::deviance(model0)
  df1 <- stats::df.residual(fit$model)
  if (has_zero_spread(fit)) {
    y <- fit$data$outcome
    changes <- rss0 - rss1 > zero_rss_rel_tol * sum((y - mean(y))^2)
    return(if (changes) list(statistic = Inf, p = 0) else list(statistic = 0, p = 1))
  }
  stat <- max(rss0 - rss1, 0) / (rss1 / df1)
  list(statistic = stat, p = stats::pf(stat, 1, df1, lower.tail = FALSE))
}

#' Effective Sample Size of a Contrast (Internal)
#'
#' `1 / sum(w_k^2 / n_k)`, the size of a single group whose mean would have the
#' same standard error as the contrast. Zero when a cell the contrast uses is
#' empty.
#'
#' @param w Contrast weights.
#' @param n Cell counts.
#' @return Numeric.
#' @keywords internal
n_eff <- function(w, n) {
  k <- w != 0
  if (any(n[k] == 0)) return(0)
  1 / sum(w[k]^2 / n[k])
}

#' Outcome-Stratified Effective Sample Size (Internal)
#'
#' The binary analogue over the eight outcome-stratified counts (Woolf). A
#' diagnostic only: it depends on the outcome, so it must never influence
#' classification or which contrast is tested.
#'
#' @param w Contrast weights.
#' @param counts8 Matrix from [cell_counts_8()].
#' @return Numeric.
#' @keywords internal
n_eff_outcome <- function(w, counts8) {
  v <- as.numeric(counts8[w != 0, ])
  if (any(v == 0)) return(0)
  1 / sum(1 / v)
}


# --- One pair ----------------------------------------------------------------

#' Analyse One Pair (Internal)
#'
#' @param x,y 0/1 feature vectors on the analysis sample.
#' @param outcome Outcome on the analysis sample; 0/1 integers if binary.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @param near_nested_max See [classify_occupancy()].
#' @param contrasts Contrasts to report, before the tested ones are added.
#' @return A list of equal-length columns, one element per emitted contrast,
#'   in catalogue order. Pair identity and `n_used` are added by the caller.
#' @keywords internal
analyse_pair <- function(x, y, outcome, outcome_type, near_nested_max, contrasts) {
  cell <- 1L + as.integer(x) + 2L * as.integer(y)
  n <- cell_counts(cell)
  cls <- classify_occupancy(n, near_nested_max)
  binary <- outcome_type == "binary"
  counts8 <- if (binary) cell_counts_8(cell, outcome) else NULL

  # At least one cell is always occupied, so the fit always exists; disjoint
  # and degenerate pairs still report the contrasts they can estimate
  fit <- fit_cell_means(cell, outcome, outcome_type, n)
  flat <- !binary && has_zero_spread(fit)

  # The user's contrasts are unioned with the tested ones, so a hypothesis
  # under test can never be dropped from the table
  active <- contrast_catalogue[contrast_catalogue %in% union(cls$tested, contrasts)]
  k <- length(active)

  beta <- se <- n_eff_col <- n_eff_o <- statistic <- p_val <- rep(NA_real_, k)
  est <- character(k)
  df <- rep(NA_integer_, k)
  role <- ifelse(active %in% cls$tested, "tested", "reported")

  for (i in seq_len(k)) {
    w <- contrast_weights[[active[i]]]
    # Aliasing first: past this point every cell the contrast uses is
    # occupied, which contrast_determinacy() relies on
    est[i] <- if (any(w != 0 & n == 0)) {
      "not_identified"
    } else if (binary) {
      contrast_determinacy(w, counts8)
    } else {
      "estimable"
    }

    if (est[i] == "estimable") {
      e <- contrast_estimate(fit$model, w, n, zero_se = flat)
      beta[i] <- e$beta
      se[i] <- e$se
    } else if (est[i] == "unbounded_positive") {
      beta[i] <- Inf
      se[i] <- Inf
    } else if (est[i] == "unbounded_negative") {
      beta[i] <- -Inf
      se[i] <- Inf
    }

    # Only tested rows get a p-value. A separated row keeps it: separation
    # damages the estimate, not the likelihood ratio
    if (role[i] == "tested" && est[i] != "not_identified") {
      tst <- contrast_test(fit, active[i], outcome_type)
      statistic[i] <- tst$statistic
      p_val[i] <- tst$p
      df[i] <- 1L
    }

    n_eff_col[i] <- n_eff(w, n)
    if (binary) n_eff_o[i] <- n_eff_outcome(w, counts8)
  }

  list(
    n00 = rep(n[["00"]], k), n10 = rep(n[["10"]], k),
    n01 = rep(n[["01"]], k), n11 = rep(n[["11"]], k),
    m_off = rep(as.integer(cls$m_off), k),
    occupancy_class = rep(cls$class, k),
    contrast = active, role = role,
    beta = beta, se = se,
    beta_estimability = est,
    n_eff = n_eff_col, n_eff_outcome = n_eff_o,
    statistic = statistic, df = df, p_val = p_val,
    # Zeros in occupied predictor cells only: an empty predictor cell is
    # aliasing, not separation
    zero_cells_8 = rep(if (binary) sum(counts8[n > 0, , drop = FALSE] == 0L) else NA_integer_, k)
  )
}
