# Helpers for the outcome-association tests. testthat sources helper-*.R first.

# A one-pair coreact result whose 2x2 table has exactly the cell counts `n`
# (named "00", "10", "01", "11"). For a binary outcome give `cases`, the number
# of outcome-1 samples per cell; otherwise give `outcome` in cell order.
pair_fixture <- function(n, cases = NULL, outcome = NULL) {
  n <- n[c("00", "10", "01", "11")]
  cells <- rep(names(n), n)
  ids <- paste0("S", seq_along(cells))
  if (!is.null(cases)) {
    cases <- cases[names(n)]
    outcome <- unlist(lapply(names(n), function(k) rep(c(1L, 0L), c(cases[[k]], n[[k]] - cases[[k]]))))
  }
  one_row <- function(v, id, name) {
    new_coreact_data(matrix(v, 1, dimnames = list(id, ids)), data.frame(row.names = id), name = name)
  }
  ox <- one_row(as.integer(substr(cells, 1, 1)), "A", "X")
  oy <- one_row(as.integer(substr(cells, 2, 2)), "B", "Y")
  list(res = suppressMessages(run_coreact(ox, oy, effect_filters = list())),
       x = ox, y = oy, outcome = stats::setNames(outcome, ids))
}

# Named cell counts, in the order "00", "10", "01", "11"
cells4 <- function(n00, n10, n01, n11) c("00" = n00, "10" = n10, "01" = n01, "11" = n11)

# Run the outcome step on a fixture, quietly
run_fixture <- function(f, ...) {
  suppressMessages(run_coreact_outcome(f$res, f$x, f$y, outcome = f$outcome, ...))
}

# A two-input result with planted structure, so that several occupancy classes
# occur: nested, identical and disjoint features as well as unstructured ones.
# Returns the objects, the result, and a binary and a continuous outcome that
# depend on some of the features.
outcome_fixture <- function(n_x = 8, n_y = 6, n_samples = 160, seed = 1) {
  withr::with_seed(seed, {
    mx <- matrix(stats::rbinom(n_x * n_samples, 1, 0.3), n_x)
    my <- matrix(stats::rbinom(n_y * n_samples, 1, 0.3), n_y)
    my[1, ] <- mx[1, ] * stats::rbinom(n_samples, 1, 0.5)          # Y1 nested in X1
    my[2, ] <- mx[2, ]                                               # Y2 identical to X2
    my[3, ] <- (1 - mx[3, ]) * stats::rbinom(n_samples, 1, 0.4)    # Y3 disjoint from X3
    my[4, ] <- pmax(mx[4, ], stats::rbinom(n_samples, 1, 0.04))    # X4 nearly nested in Y4
    ids <- paste0("S", seq_len(n_samples))
    dimnames(mx) <- list(paste0("X", seq_len(n_x)), ids)
    dimnames(my) <- list(paste0("Y", seq_len(n_y)), ids)
    lp <- -0.6 + 0.8 * mx[1, ] + 0.5 * my[5, ] + 0.9 * mx[5, ] * my[5, ]
    binary <- stats::setNames(stats::rbinom(n_samples, 1, stats::plogis(lp)), ids)
    continuous <- stats::setNames(stats::rnorm(n_samples, lp), ids)
  })
  ox <- new_coreact_data(mx, data.frame(row.names = rownames(mx)), name = "X")
  oy <- new_coreact_data(my, data.frame(row.names = rownames(my)), name = "Y")
  list(x = ox, y = oy, res = suppressMessages(run_coreact(ox, oy, effect_filters = list())),
       binary = binary, continuous = continuous)
}

# Pair index of each output row, in output order
pair_index <- function(out) {
  key <- paste(out$feature_x, out$feature_y, sep = "\r")
  match(key, unique(key))
}

# Brute-force reference for one pair's tested contrasts: likelihood ratio or F
# tests between treatment-coded models, a different parameterisation from the
# engine's cell-means fit and merged levels.
reference_tests <- function(x, y, o, type, occupancy_class) {
  d <- data.frame(o = o, x = x, y = y, x_only = x * (1 - y), y_only = (1 - x) * y)
  fit <- function(f) {
    if (type == "binary") suppressWarnings(stats::glm(f, family = stats::binomial(), data = d))
    else stats::lm(f, data = d)
  }
  test <- function(full, reduced) {
    m1 <- fit(full)
    m0 <- fit(reduced)
    if (type == "binary") {
      stats::pchisq(max(stats::deviance(m0) - stats::deviance(m1), 0), 1, lower.tail = FALSE)
    } else {
      stats::anova(m0, m1)[["Pr(>F)"]][2]
    }
  }
  switch(occupancy_class,
    full = list(interaction = test(o ~ x * y, o ~ x + y)),
    near_nested = list(interaction = test(o ~ x * y, o ~ x + y),
                       "11_vs_10" = test(o ~ x * y, o ~ x + y_only),
                       "11_vs_01" = test(o ~ x * y, o ~ y + x_only)),
    nested_y_in_x = list("11_vs_10" = test(o ~ x + y, o ~ x)),
    nested_x_in_y = list("11_vs_01" = test(o ~ x + y, o ~ y)),
    identical = list("11_vs_00" = test(o ~ x, o ~ 1)),
    no_double_negative = list("11_vs_10" = test(o ~ x + y, o ~ x),
                              "11_vs_01" = test(o ~ x + y, o ~ y)),
    list()
  )
}

# Closed-form estimate of a contrast of cell means: log-odds for a binary
# outcome, plain means for a continuous one. NA when a cell it uses is empty or,
# for a binary outcome, has only cases or only controls.
closed_form_beta <- function(w, x, y, o, type) {
  cell <- 1 + x + 2 * y
  used <- which(w != 0)
  vals <- vapply(used, function(k) {
    ok <- o[cell == k]
    if (length(ok) == 0) return(NA_real_)
    if (type == "continuous") return(mean(ok))
    if (all(ok == 1) || all(ok == 0)) return(NA_real_)
    log(mean(ok) / (1 - mean(ok)))
  }, numeric(1))
  sum(w[used] * vals)
}
