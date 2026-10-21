#' Associate Co-reactive Pairs With a Sample-Level Outcome
#'
#' For each pair in a coreact result, estimates how the two features together
#' associate with an outcome measured on the same samples. The question is the
#' interaction — do x and y *together* go with the outcome — or, where the
#' pair's 2x2 table cannot support an interaction, the contrast of cell means
#' that carries the equivalent question.
#'
#' @param res A result of [run_coreact()] or [coreact_pipeline()], or one read
#'   back with [read_coreact_results()]. It must still carry its `provenance`
#'   attribute (see Provenance below).
#' @param data_x,data_y The `coreact_data` objects `res` was computed on, in the
#'   same order. Omit `data_y` for a within-library result, as for
#'   [run_coreact()]; supplying the same data twice is an error.
#' @param outcome A vector named by sample ID (`colnames(data_x$mat)`); samples
#'   are matched by name, never by position. Binary outcomes are coded `0` / `1`,
#'   logical, or a two-level factor whose second level is the event. `NA` marks
#'   samples to leave out.
#' @param outcome_type `"binary"` (logistic model) or `"continuous"` (linear
#'   model).
#' @param covariates Reserved for covariate adjustment; must be `NULL`.
#' @param near_nested_max Count. A fully occupied pair whose thinner
#'   off-diagonal cell (`min(n01, n10)`) holds fewer samples than this is
#'   `near_nested`, and is tested on a conditional contrast as well as the
#'   interaction. Always read as a count, never as a fraction; it may also be a
#'   function of `n_used`, the number of samples with an outcome, that returns
#'   one. Fix it before looking
#'   at results: it sets the size of the tested family, so adjusted p-values
#'   from runs with different values are not comparable.
#' @param contrasts Which contrasts to report, from `"interaction"`,
#'   `"11_vs_10"`, `"11_vs_01"`, `"11_vs_00"`, `"10_vs_00"`, `"01_vs_00"`,
#'   `"10_vs_01"`. Each pair's tested contrasts are always added, so narrowing
#'   this never removes a hypothesis under test.
#' @param ci Which rows get a 95% confidence interval: `"tested"` (the default)
#'   for tested rows only, `"all"` for every estimable row, `"none"` for none.
#'   On a binary outcome each interval is a profile-likelihood interval that
#'   costs several times the pair's fit and test, so `"all"` is markedly slower;
#'   see Intervals below.
#' @param fdr_method Multiple-testing adjustment, any [stats::p.adjust()]
#'   method.
#' @param n_cores,chunk_size As in [run_coreact()]: a fork cluster above one
#'   core, serial with a warning on Windows. Chunks are rows of `res`.
#'
#' @details
#' ## Estimands
#'
#' The outcome is modelled by one mean per occupied cell of the pair's 2x2
#' table (cells `00`, `10`, `01`, `11`, x then y), on the log-odds scale for a
#' binary outcome and the outcome's own scale for a continuous one. Every
#' reported quantity is a contrast of those means: `interaction` is
#' `mu11 - mu10 - mu01 + mu00`, and `a_vs_b` is `mu_a - mu_b`. A contrast is
#' estimable when every cell it uses is occupied, so which contrasts exist for
#' a pair depends on its occupancy, not on a fitted model.
#'
#' Each pair is classified from its four cell counts alone, and the class
#' decides what is tested: the interaction for a fully occupied pair; for a
#' nested pair (one off-diagonal cell empty) the conditional contrast on the
#' occupied side; `11_vs_00` for `identical` features; nothing for `disjoint`
#' (no double positives) or `degenerate` pairs, which are still reported.
#' `near_nested` pairs are tested on the interaction and on the conditional
#' contrast that avoids the thin cell; `no_double_negative` pairs on both
#' conditional contrasts. Because the choice never looks at the outcome, it
#' costs no multiple-testing budget.
#'
#' ## Testing
#'
#' Only `tested` rows have p-values: a 1-df likelihood ratio test for a binary
#' outcome, the exact F test for a continuous one, each against the cell-means
#' model with the contrast set to zero. Both are two-sided, so a negative
#' interaction is found as readily as a positive one. Wald tests are not used:
#' they lose power badly as a logistic fit approaches separation. Adjusted
#' p-values are computed once over all tested rows of all pairs.
#'
#' The tested family mixes hypotheses, so `p_adj <= 0.05` says that a *pair* is
#' associated with the outcome, not that it has an interaction:
#' `occupancy_class` and `contrast` say which claim a row supports. A
#' `near_nested` pair's two tested rows are strongly correlated and may both be
#' significant; `any_significant` counts the pair once.
#'
#' ## Separation
#'
#' When an occupied cell's samples are all cases or all controls, a logistic
#' estimate diverges. This is decided from the counts, which is exact here,
#' not from how close `glm()` gets to a fitted probability of 0 or 1. A
#' contrast that depends on such a cell is reported as `beta = Inf` or `-Inf`
#' with `se = Inf` when its direction is determined (`unbounded_positive`,
#' `unbounded_negative`), and as `NA` when separated cells pull it both ways
#' (`indeterminate`). An empty predictor cell is a different matter: contrasts
#' using it are `not_identified`. A separated row keeps its p-value and stays
#' in the family: separation damages the estimate, not the likelihood ratio,
#' so significant rows with an infinite `beta` are expected. `n_eff` stays
#' finite on such rows; it measures the predictor cells, which are intact.
#'
#' ## Intervals
#'
#' `ci_lo` and `ci_hi` bound a 95% confidence interval that agrees with the
#' test: it excludes zero exactly when the row's p-value is below 0.05. For a
#' binary outcome it is a profile-likelihood interval, from a refit in which
#' the contrast is a single coefficient; for a continuous outcome the exact t
#' interval. Wald intervals (`beta +/- 1.96 * se`) are not used. Rows without
#' a finite estimate have no interval, and nor, for now, do separated rows:
#' profiling a diverging coefficient gives meaningless bounds.
#'
#' A profile interval costs several times a pair's fit and test, which is why
#' only tested rows get one by default. `ci = "all"` computes one for every
#' estimable row, up to seven per pair; `"none"` skips them.
#'
#' ## Continuous outcomes
#'
#' The linear model pools one residual standard deviation over the four
#' cells, so the F test and every standard error assume the cells are equally
#' variable. Coreact's cells are typically unbalanced, with few
#' double-positive samples: if a small cell is more variable than the large
#' ones, contrasts using it get standard errors that are too small and
#' p-values that are too small. A continuous outcome taking only two values
#' draws a warning, since a binary outcome passed as continuous fits a linear
#' probability model without complaint. If every occupied cell of a pair is
#' constant, its F test returns its limit: `p = 0` when the tested contrast is
#' non-zero, `p = 1` when it is zero.
#'
#' ## Within a library
#'
#' For a result of `run_coreact(data_x)`, which feature of a pair is `x` and
#' which is `y` follows index order and means nothing else. The directional
#' labels — `nested_y_in_x` against `nested_x_in_y`, `11_vs_10` against
#' `11_vs_01`, the sign of `10_vs_01` — then describe orientation only.
#'
#' ## Provenance
#'
#' `res` must carry the `provenance` attribute that [run_coreact()] attaches,
#' which the package uses to check that `data_x` and `data_y` are the objects
#' `res` was computed on and are given in the right order: by name, size,
#' number of positive cells and an MD5 hash of their contents. You need not use it
#' yourself. It is restored by [read_coreact_results()], but lost when results
#' are read with `read_tsv()`, summarised or passed through base `merge()`.
#'
#' Samples with a missing outcome are dropped before anything is counted, so
#' every count refers to the analysis sample (`n_used`). Dropping samples can
#' leave a pair below the prevalence or intersection threshold it was selected
#' under; a message says how many.
#'
#' @return A tibble with one row per pair and reported contrast, never
#'   truncated or filtered: `feature_x`, `feature_y`; the cell counts `n00`,
#'   `n10`, `n01`, `n11` and `m_off = min(n01, n10)` on the analysis sample of
#'   `n_used` samples; `occupancy_class`; `contrast`; `role` (`tested` or
#'   `reported`); `beta` and `se`; `ci_lo`, `ci_hi`; `beta_estimability`
#'   (`estimable`, `unbounded_positive`, `unbounded_negative`, `indeterminate`,
#'   `not_identified`); `n_eff`, the effective sample size of the contrast from
#'   the predictor cells, and for a binary outcome `n_eff_outcome`, from the
#'   outcome-stratified cells (a diagnostic only); `statistic` (deviance
#'   difference or F), `df`, `p_val` and `p_adj`, set on tested rows only;
#'   `any_significant`, whether any of the pair's tested rows has
#'   `p_adj <= 0.05`, `NA` when the pair has none; and for a binary outcome
#'   `zero_cells_8`, the number of zero outcome counts (no cases, or no
#'   controls) in the pair's occupied predictor cells.
#'
#'   Attributes record `near_nested_max` (resolved), `contrasts`, `ci`,
#'   `fdr_method`, `outcome_type`, `within_library`, `frac_unbounded_tested`
#'   (share of tested rows with an infinite estimate),
#'   `n_pairs_below_selection`, and `coreact_provenance`, the provenance of
#'   `res`.
#' @export
#' @examples
#' path_x <- system.file("extdata", "input_x.tsv", package = "coreact")
#' x <- import_coreact_tsv(path_x, meta_cols = 1:4, feature_id = "u_pep_id", name = "X")
#'
#' res <- run_coreact(x)                       # pairs within one library
#'
#' set.seed(1)
#' outcome <- stats::setNames(stats::rbinom(ncol(x$mat), 1, 0.5), colnames(x$mat))
#'
#' out <- run_coreact_outcome(res, x, outcome = outcome)
#' tested <- out[out$role == "tested", ]
#' head(tested[, c("feature_x", "feature_y", "occupancy_class", "contrast",
#'                 "beta", "beta_estimability", "p_val", "p_adj")])
run_coreact_outcome <- function(res,
                                data_x,
                                data_y = NULL,
                                outcome,
                                outcome_type = c("binary", "continuous"),
                                covariates = NULL,
                                near_nested_max = 10,
                                contrasts = c("interaction", "11_vs_10", "11_vs_01", "11_vs_00",
                                              "10_vs_00", "01_vs_00", "10_vs_01"),
                                ci = c("tested", "all", "none"),
                                fdr_method = "BH",
                                n_cores = 1,
                                chunk_size = NULL) {

  # --- 0. Arguments that need no data ---
  outcome_type <- match.arg(outcome_type)
  ci <- match.arg(ci)
  if (!is.null(covariates)) {
    stop("'covariates' is reserved for covariate adjustment, which is not yet available; it must be NULL.")
  }
  contrasts <- validate_contrasts(contrasts)
  validate_fdr_method(fdr_method)
  validate_parallel_args(n_cores, chunk_size)

  # --- 1. The result, its provenance and the objects it was computed on ---
  needed <- c("feature_x", "feature_y", "idx_x", "idx_y")
  if (!is.data.frame(res) || !all(needed %in% names(res))) {
    stop("'res' must be a coreact result with columns feature_x, feature_y, idx_x and idx_y.")
  }
  prov <- attr(res, "provenance")
  if (is.null(prov)) {
    stop("'res' carries no provenance attribute, which is needed to check that data_x and ",
         "data_y are the objects it was computed on. run_coreact() and coreact_pipeline() ",
         "attach it and read_coreact_results() restores it; read_tsv(), summarise() and ",
         "merge() drop it. Read the results with read_coreact_results().")
  }
  if (!inherits(data_x, "coreact_data")) {
    stop("'data_x' must be a 'coreact_data' object (see new_coreact_data()).")
  }
  # data_y sits next to data_x, as in run_coreact(), so a positional call
  # run_coreact_outcome(res, x, outcome) binds the outcome to data_y
  if (!is.null(data_y)) {
    if (!inherits(data_y, "coreact_data")) {
      stop("'data_y' is not a 'coreact_data' object. For a within-library result, ",
           "omit 'data_y' and pass 'outcome' by name.")
    }
    if (identical(data_x$mat, data_y$mat)) {
      stop("'data_x' and 'data_y' hold the same data. For a within-library result, ",
           "supply 'data_x' only.")
    }
  }

  # The mode is declared by omitting data_y and must match the one recorded.
  # Without this check a two-input result with data_y forgotten would be looked
  # up against data_x on both sides, and under row-number feature IDs the index
  # check below can pass.
  within_library <- is.null(data_y)
  if (!identical(within_library, isTRUE(prov$within_library))) {
    stop(if (within_library) {
      "'res' is a two-input coreact result: supply 'data_y' as well."
    } else {
      "'res' is a within-library coreact result: omit 'data_y'."
    })
  }
  # The fingerprints catch data_x and data_y passed in the wrong order, which
  # the index check misses when both carry row-number IDs
  check_fingerprint(data_x, prov, "x")
  if (within_library) {
    data_y <- data_x
  } else {
    check_fingerprint(data_y, prov, "y")
    if (!identical(colnames(data_x$mat), colnames(data_y$mat))) {
      stop("The sample columns of 'data_x' and 'data_y' differ.")
    }
  }
  check_lookup(res, data_x, data_y)
  # Coreact reports each within-library pair once, so this should never fire;
  # it guards a result built or edited by hand
  if (within_library && any(res$idx_x >= res$idx_y)) {
    stop("Within-library input is not reduced to one orientation per pair: expected ",
         "idx_x < idx_y throughout. Was this result produced by run_coreact() on a single input?")
  }
  if (anyDuplicated(data.frame(res$idx_x, res$idx_y))) {
    stop("'res' contains the same pair more than once.")
  }

  # --- 2. Outcome, aligned by sample ID ---
  samples <- colnames(data_x$mat)
  outcome <- align_outcome(outcome, samples, outcome_type)
  if (!is.function(near_nested_max)) check_near_nested_max(near_nested_max)

  # --- 3. Fix the analysis sample before anything is counted ---
  # glm() and lm() drop incomplete cases silently, so a cell counted on all
  # samples can be empty in the fit, its coefficient absent rather than NA.
  keep <- !is.na(outcome)
  n_used <- sum(keep)
  if (n_used == 0) stop("No sample has a non-missing outcome.")
  outcome <- unname(outcome[keep])
  check_analysis_outcome(outcome, outcome_type)
  # Everything that shapes the family is computed on the analysis sample, so a
  # function-valued cutoff is applied to its size
  if (is.function(near_nested_max)) {
    near_nested_max <- near_nested_max(n_used)
    check_near_nested_max(near_nested_max)
  }
  mat_x <- data_x$mat[, keep, drop = FALSE]
  mat_y <- if (within_library) mat_x else data_y$mat[, keep, drop = FALSE]

  # --- 4. Every pair, in chunks of rows of res ---
  idx_x <- as.integer(res$idx_x)
  idx_y <- as.integer(res$idx_y)
  run_chunk <- function(rows) {
    xs <- as.matrix(mat_x[idx_x[rows], , drop = FALSE])
    ys <- as.matrix(mat_y[idx_y[rows], , drop = FALSE])
    lapply(seq_along(rows), function(j) {
      analyse_pair(xs[j, ], ys[j, ], outcome, outcome_type, near_nested_max, contrasts, ci)
    })
  }
  pieces <- map_pair_chunks(seq_len(nrow(res)), run_chunk, n_cores, chunk_size)
  pair <- rep(seq_along(pieces), vapply(pieces, function(p) length(p$contrast), integer(1)))
  out <- assemble_outcome(pieces, pair, as.character(res$feature_x), as.character(res$feature_y),
                          n_used)

  # --- 5. Multiple testing over the tested rows only ---
  out <- adjust_outcome_fdr(out, pair, method = fdr_method)
  out$any_significant <- pair_any_significant(out, pair)

  # --- 6. Diagnostics, reported rather than guessed ---
  tested <- out$role == "tested"
  frac_unbounded <- if (any(tested)) {
    mean(grepl("^unbounded", as.character(out$beta_estimability[tested])))
  } else {
    NA_real_
  }
  if (isTRUE(frac_unbounded > 0.2)) {
    message(sprintf("%.0f%% of tested rows are separated; their effect sizes are unbounded.",
                    100 * frac_unbounded))
  }
  n_below <- count_below_selection(out, pair, prov, within_library)
  if (n_used < length(samples) && n_below > 0) {
    message(sprintf(paste0("After leaving out %d sample(s) with a missing outcome, %d of %d pairs ",
                           "fall below the prevalence or intersection threshold they were selected under."),
                    length(samples) - n_used, n_below, nrow(res)))
  }

  # --- 7. What the table was computed with ---
  # near_nested_max sets the family size and so every adjusted p-value
  attr(out, "near_nested_max") <- near_nested_max
  attr(out, "contrasts") <- contrasts
  attr(out, "ci") <- ci
  attr(out, "fdr_method") <- fdr_method
  attr(out, "outcome_type") <- outcome_type
  attr(out, "within_library") <- within_library
  attr(out, "frac_unbounded_tested") <- frac_unbounded
  attr(out, "n_pairs_below_selection") <- n_below
  attr(out, "coreact_provenance") <- prov
  out
}


# --- Validation --------------------------------------------------------------

#' Validate the Requested Contrasts (Internal)
#'
#' @param contrasts Character vector of contrast names.
#' @return The contrasts, deduplicated and in catalogue order.
#' @keywords internal
validate_contrasts <- function(contrasts) {
  if (!is.character(contrasts) || length(contrasts) == 0 || anyNA(contrasts)) {
    stop("'contrasts' must be a non-empty character vector of contrast names.")
  }
  unknown <- setdiff(contrasts, contrast_catalogue)
  if (length(unknown) > 0) {
    stop(sprintf("Unknown contrast(s) %s. Choose from %s.",
                 paste(sprintf("\"%s\"", unknown), collapse = ", "),
                 paste(sprintf("\"%s\"", contrast_catalogue), collapse = ", ")))
  }
  contrast_catalogue[contrast_catalogue %in% contrasts]
}

#' Validate the Multiple-Testing Method (Internal)
#'
#' @param fdr_method A [stats::p.adjust()] method name.
#' @return `NULL`, invisibly.
#' @keywords internal
validate_fdr_method <- function(fdr_method) {
  if (!is.character(fdr_method) || length(fdr_method) != 1 ||
        !fdr_method %in% stats::p.adjust.methods) {
    stop(sprintf("'fdr_method' must be one of %s.",
                 paste(sprintf("\"%s\"", stats::p.adjust.methods), collapse = ", ")))
  }
  invisible(NULL)
}

#' Validate a Resolved `near_nested_max` (Internal)
#'
#' @param value The count, given directly or returned by the function given.
#' @return `NULL`, invisibly.
#' @keywords internal
check_near_nested_max <- function(value) {
  if (!is.numeric(value) || length(value) != 1 || is.na(value) || value < 0) {
    stop("'near_nested_max' must be a single non-negative count, or a function of the ",
         "number of samples used that returns one.")
  }
  invisible(NULL)
}

#' Check an Object Against Its Recorded Fingerprint (Internal)
#'
#' Compares the object's name, number of features and number of positive cells
#' with those [run_coreact()] recorded for that side, then the hash of its
#' contents ([matrix_md5()]). The counts come first because a mismatch there
#' says how the objects differ; the hash then catches what the counts miss, such
#' as a positive moved from one sample to another. Numbers are compared as
#' numbers: read back from a sidecar they may be double where they were
#' integer.
#'
#' @param obj A `coreact_data` object.
#' @param prov The provenance list.
#' @param side `"x"` or `"y"`.
#' @param what How to refer to the object in the error.
#' @param remedy What to do about a mismatch, for the error.
#' @return `NULL`, invisibly; errors on a mismatch.
#' @keywords internal
check_fingerprint <- function(obj, prov, side, what = sprintf("'data_%s'", side),
                              remedy = "Pass the inputs given to run_coreact(), in the same order.") {
  rec <- function(field) prov[[paste0(side, "_", field)]]
  same_number <- function(a, b) length(b) == 1 && isTRUE(as.numeric(a) == as.numeric(b))
  problems <- c(
    if (!identical(obj$name, as.character(rec("name")))) {
      sprintf("name '%s' where '%s' was recorded", obj$name, paste(rec("name"), collapse = ""))
    },
    if (!same_number(nrow(obj$mat), rec("n_features"))) {
      sprintf("%d features where %s were recorded", nrow(obj$mat), format(rec("n_features")))
    },
    if (!same_number(sum(obj$mat), rec("n_positive"))) {
      sprintf("%s positive cells where %s were recorded", format(sum(obj$mat)), format(rec("n_positive")))
    }
  )
  if (length(problems) == 0 && !identical(matrix_md5(obj$mat), as.character(rec("md5")))) {
    problems <- "the same name and counts as recorded, but different contents"
  }
  if (length(problems) > 0) {
    stop(sprintf("%s is not the data this result was computed on: %s. %s",
                 what, paste(problems, collapse = "; "), remedy))
  }
  invisible(NULL)
}

#' Validate the Feature Lookup (Internal)
#'
#' `idx_x` / `idx_y` index the objects coreact ran on. They are trusted only
#' after the IDs they point at are checked against `feature_x` / `feature_y`;
#' a mismatch is an error naming the first offending pair, never a fall-back
#' to matching by name.
#'
#' @param res A coreact result.
#' @param data_x,data_y `coreact_data` objects.
#' @return `NULL`, invisibly.
#' @keywords internal
check_lookup <- function(res, data_x, data_y) {
  in_range <- function(idx, n) {
    is.numeric(idx) && !anyNA(idx) && all(idx >= 1) && all(idx <= n) && all(idx == round(idx))
  }
  if (!in_range(res$idx_x, nrow(data_x$mat)) || !in_range(res$idx_y, nrow(data_y$mat))) {
    stop("'res' has idx_x / idx_y values that are not row numbers of data_x / data_y: these ",
         "are not the objects this result was computed on.")
  }
  found_x <- rownames(data_x$mat)[res$idx_x]
  found_y <- rownames(data_y$mat)[res$idx_y]
  bad <- which(found_x != as.character(res$feature_x) | found_y != as.character(res$feature_y))
  if (length(bad) > 0) {
    i <- bad[1]
    stop(sprintf(paste0("Row %d of 'res' is the pair (%s, %s), but its indices point at (%s, %s) ",
                        "in data_x / data_y: these are not the objects this result was computed on."),
                 i, res$feature_x[i], res$feature_y[i], found_x[i], found_y[i]))
  }
  invisible(NULL)
}

#' Align and Code the Outcome (Internal)
#'
#' Matches the outcome to the samples by name and codes a binary outcome as
#' 0 / 1 integers, with a two-level factor's second level as 1.
#'
#' @param outcome Named vector.
#' @param samples Sample IDs, in matrix column order.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @return Named vector in `samples` order, `NA` where the outcome is missing.
#' @keywords internal
align_outcome <- function(outcome, samples, outcome_type) {
  ids <- names(outcome)
  if (is.null(ids) || anyNA(ids)) {
    stop("'outcome' must be a vector named by sample ID (the column names of data_x$mat).")
  }
  few <- function(v) paste(utils::head(v, 5), collapse = ", ")
  dup <- unique(ids[duplicated(ids)])
  if (length(dup) > 0) {
    stop(sprintf("'outcome' names these samples more than once: %s.", few(dup)))
  }
  absent <- setdiff(samples, ids)
  extra <- setdiff(ids, samples)
  if (length(absent) > 0 || length(extra) > 0) {
    stop(sprintf("'outcome' does not match the samples of data_x:%s%s",
                 if (length(absent)) sprintf(" no value for %s.", few(absent)) else "",
                 if (length(extra)) sprintf(" unknown samples %s.", few(extra)) else ""))
  }
  outcome <- outcome[samples]

  if (outcome_type == "continuous") {
    if (!is.numeric(outcome)) stop("A continuous outcome must be numeric.")
    return(stats::setNames(as.numeric(outcome), samples))
  }

  coded <- if (is.factor(outcome)) {
    if (nlevels(outcome) != 2) {
      stop(sprintf("A binary outcome given as a factor must have two levels, not %d.", nlevels(outcome)))
    }
    as.integer(outcome) - 1L
  } else if (is.logical(outcome)) {
    as.integer(outcome)
  } else if (is.numeric(outcome)) {
    other <- unique(outcome[!is.na(outcome) & !outcome %in% c(0, 1)])
    if (length(other) > 0) {
      stop(sprintf("A binary outcome must be coded 0/1; found %s.", few(other)))
    }
    as.integer(outcome)
  } else {
    stop("A binary outcome must be coded 0/1, logical, or a two-level factor.")
  }
  stats::setNames(coded, samples)
}

#' Check the Outcome on the Analysis Sample (Internal)
#'
#' Run after samples with a missing outcome are dropped, which can change the
#' answer.
#'
#' @param outcome Outcome on the analysis sample.
#' @param outcome_type `"binary"` or `"continuous"`.
#' @return `NULL`, invisibly.
#' @keywords internal
check_analysis_outcome <- function(outcome, outcome_type) {
  values <- unique(outcome)
  if (outcome_type == "binary") {
    if (length(values) < 2) {
      stop("The binary outcome takes a single value in the analysis sample, so no ",
           "association can be estimated.")
    }
    return(invisible(NULL))
  }
  if (!all(is.finite(outcome))) stop("A continuous outcome must be finite; found Inf in the analysis sample.")
  if (length(values) < 2) stop("The continuous outcome is constant in the analysis sample.")
  # Four cell means leave a residual degree of freedom for every pair only
  # with at least five samples
  if (length(outcome) < 5) {
    stop("A continuous outcome needs at least 5 samples with a value, so that every ",
         "pair's F test has a residual degree of freedom.")
  }
  if (length(values) == 2) {
    warning(sprintf(paste0("outcome_type = \"continuous\" but the outcome takes only two values ",
                           "(%s). A binary outcome passed as continuous fits a linear probability ",
                           "model without complaint; use outcome_type = \"binary\" if that is what it is."),
                    paste(format(sort(values)), collapse = ", ")),
            call. = FALSE)
  }
  invisible(NULL)
}


# --- Execution and assembly --------------------------------------------------

#' Map a Function Over Chunks of Pairs (Internal)
#'
#' The parallel pattern of [run_coreact()], over rows of the result rather
#' than features: a fork cluster above one core, serial with a warning on
#' Windows.
#'
#' @param rows Row numbers of `res`.
#' @param fun Function of a vector of rows, returning a list with one element
#'   per row.
#' @param n_cores,chunk_size As for [run_coreact()].
#' @return A list with one element per row, in row order.
#' @keywords internal
map_pair_chunks <- function(rows, fun, n_cores, chunk_size) {
  if (length(rows) == 0) return(list())
  if (is.null(chunk_size)) {
    chunk_size <- max(50, ceiling(length(rows) / max(n_cores * 20, 10)))
    message(sprintf("Auto-chunking: Size = %d pairs (%d total chunks).",
                    as.integer(chunk_size), as.integer(ceiling(length(rows) / chunk_size))))
  }
  chunks <- split(rows, ceiling(seq_along(rows) / as.integer(chunk_size)))

  done <- if (n_cores > 1 && .Platform$OS.type == "windows") {
    warning("Fork clusters not supported on Windows. Running serially.")
    lapply(chunks, fun)
  } else if (n_cores > 1) {
    cl <- parallel::makeForkCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::parLapply(cl, chunks, fun)
  } else {
    lapply(chunks, fun)
  }
  unlist(done, recursive = FALSE, use.names = FALSE)
}

#' Assemble the Per-Pair Results Into the Output Table (Internal)
#'
#' @param pieces List of per-pair column lists from [analyse_pair()].
#' @param pair Integer pair index for each output row.
#' @param feature_x,feature_y Feature IDs, one per pair.
#' @param n_used Size of the analysis sample.
#' @return The output tibble, with `p_adj` and `any_significant` still `NA`.
#' @keywords internal
assemble_outcome <- function(pieces, pair, feature_x, feature_y, n_used) {
  if (length(pieces) == 0) return(empty_outcome_tibble())
  cols <- data.table::rbindlist(pieces)
  data.table::setDF(cols)
  k <- nrow(cols)
  tibble::tibble(
    feature_x = feature_x[pair], feature_y = feature_y[pair],
    n00 = cols$n00, n10 = cols$n10, n01 = cols$n01, n11 = cols$n11,
    m_off = cols$m_off, n_used = rep(as.integer(n_used), k),
    occupancy_class = factor(cols$occupancy_class, levels = occupancy_levels),
    contrast = factor(cols$contrast, levels = contrast_catalogue),
    role = factor(cols$role, levels = c("tested", "reported")),
    beta = cols$beta, se = cols$se,
    ci_lo = cols$ci_lo, ci_hi = cols$ci_hi,
    beta_estimability = factor(cols$beta_estimability, levels = estimability_levels),
    n_eff = cols$n_eff, n_eff_outcome = cols$n_eff_outcome,
    statistic = cols$statistic, df = cols$df, p_val = cols$p_val,
    p_adj = rep(NA_real_, k), any_significant = rep(NA, k),
    zero_cells_8 = cols$zero_cells_8
  )
}

#' Empty Outcome Table (Internal)
#'
#' The zero-row table for a result with no pairs, with the same columns and
#' types as a non-empty one.
#'
#' @return A zero-row tibble.
#' @keywords internal
empty_outcome_tibble <- function() {
  tibble::tibble(
    feature_x = character(), feature_y = character(),
    n00 = integer(), n10 = integer(), n01 = integer(), n11 = integer(),
    m_off = integer(), n_used = integer(),
    occupancy_class = factor(character(), levels = occupancy_levels),
    contrast = factor(character(), levels = contrast_catalogue),
    role = factor(character(), levels = c("tested", "reported")),
    beta = numeric(), se = numeric(), ci_lo = numeric(), ci_hi = numeric(),
    beta_estimability = factor(character(), levels = estimability_levels),
    n_eff = numeric(), n_eff_outcome = numeric(), statistic = numeric(),
    df = integer(), p_val = numeric(), p_adj = numeric(),
    any_significant = logical(), zero_cells_8 = integer()
  )
}


# --- Multiple testing --------------------------------------------------------

#' Adjust P-values Over the Tested Rows Only (Internal)
#'
#' @description
#' The family is the set of `tested` rows, delimited by `role` and never by
#' which p-values happen to be missing. Not [adjust_fdr()]: that adjusts the
#' co-occurrence family, over `n_pairs`.
#'
#' @details
#' A guard computes the expected family size independently, from the occupancy
#' classes: one test for `identical`, nested and `full` pairs, two for
#' `near_nested` and `no_double_negative`, none for `disjoint` and
#' `degenerate`. A mismatch means a class's tested set changed without the
#' family following, which would shift every adjusted p-value invisibly, so it
#' stops the run.
#'
#' @param out The output table.
#' @param pair Integer pair index per row.
#' @param method A [stats::p.adjust()] method.
#' @return `out` with `p_adj` set on tested rows.
#' @keywords internal
adjust_outcome_fdr <- function(out, pair, method = "BH") {
  family <- out$role == "tested" & !is.na(out$p_val)
  out$p_adj <- NA_real_
  out$p_adj[family] <- stats::p.adjust(out$p_val[family], method = method)

  pair_class <- as.character(out$occupancy_class)[!duplicated(pair)]
  n_expected <- sum(pair_class %in% one_test_classes) + 2 * sum(pair_class %in% two_test_classes)
  if (sum(family) != n_expected) {
    stop(sprintf(paste0("Internal error: %d tested rows have p-values, but the occupancy classes ",
                        "imply a family of %d. Please report this."),
                 sum(family), n_expected))
  }
  out
}

#' Whether Any of a Pair's Tested Rows Is Significant (Internal)
#'
#' `NA`, not `FALSE`, for a pair with no tested rows: "never tested" is not
#' "tested and null".
#'
#' @param out The output table, with `p_adj` set.
#' @param pair Integer pair index per row.
#' @param alpha Threshold on `p_adj`.
#' @return Logical vector, one per row.
#' @keywords internal
pair_any_significant <- function(out, pair, alpha = 0.05) {
  tested <- out$role == "tested"
  # <=, the comparison coreact's own FDR filter uses
  hit <- tested & !is.na(out$p_adj) & out$p_adj <= alpha
  any_tested <- as.vector(tapply(tested, pair, any))
  any_hit <- as.vector(tapply(hit, pair, any))
  ifelse(any_tested, any_hit, NA)[pair]
}

#' Count Pairs That No Longer Meet Their Selection Thresholds (Internal)
#'
#' Coreact selected the pairs on all samples; on the analysis sample a
#' feature's prevalence or a pair's intersection can fall below the threshold
#' it passed. Thresholds recorded as unknown are not checked.
#'
#' @param out The output table.
#' @param pair Integer pair index per row.
#' @param prov The provenance of the coreact result.
#' @param within_library Logical.
#' @return Integer count of pairs.
#' @keywords internal
count_below_selection <- function(out, pair, prov, within_library) {
  first <- !duplicated(pair)
  if (!any(first)) return(0L)
  n10 <- out$n10[first]
  n01 <- out$n01[first]
  n11 <- out$n11[first]
  threshold <- function(v) if (length(v) == 1 && is.finite(v)) v else 0
  min_x <- threshold(prov$x_min_prevalence_count)
  min_y <- threshold(if (within_library) prov$x_min_prevalence_count else prov$y_min_prevalence_count)
  min_int <- threshold(prov$effect_min_intersection)
  as.integer(sum(n10 + n11 < min_x | n01 + n11 < min_y | n11 < min_int))
}
