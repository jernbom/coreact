#' Internal Worker: Compute Co-occurrence Stats for a Chunk
#'
#' Computes intersection, union, Jaccard index, overlap coefficient and a
#' one-sided hypergeometric p-value for every (row of X chunk) x (row of Y)
#' pair, then keeps the pairs that pass the effect-size thresholds in
#' `config` and, if given, the raw p-value ceiling `max_p_val`.
#'
#' @param idx_x_chunk Vector of absolute indices for the rows of X in this chunk.
#' @param mat_x_chunk The actual data sub-matrix for X (dense or sparse).
#' @param mat_y The full matrix Y (dense or sparse).
#' @param prev_x_chunk Vector of prevalence (rowSums) for rows in X chunk.
#' @param prev_y Vector of prevalence (rowSums) for all rows in Y.
#' @param n_samples Integer, total number of samples (ncol).
#' @param config List of name-value pairs of effect-size thresholds. Current
#'   valid options are: `min_intersection`, `min_jaccard`, `min_overlap`. A
#'   threshold of 0 (or a missing entry) disables that filter.
#' @param max_p_val Numeric or `NULL`. If given, pairs whose raw p-value
#'   exceeds it are dropped. Because a BH-adjusted p-value is never smaller
#'   than the raw p-value, dropping pairs with `p_val > fdr_threshold` cannot
#'   change which pairs are significant at that threshold; it only saves
#'   memory.
#' @param upper_triangle Logical. For a within-library comparison, where X and
#'   Y are the same matrix, keep only pairs with `idx_x < idx_y`. That drops
#'   the self-pairs (a feature against itself) and one of each mirrored pair,
#'   which carry the same hypothesis with x and y exchanged.
#'
#' @details
#' The Jaccard index is undefined when both features are absent from every
#' sample (union of size zero), and the overlap coefficient is undefined when
#' either feature is absent from every sample. Such ratios are reported as
#' `NA` and never pass a positive `min_jaccard` / `min_overlap` threshold.
#'
#' @return A tibble of passing pairs with stats, or `NULL` if none pass.
#' @keywords internal
worker_chunk_calc <- function(idx_x_chunk,
                              mat_x_chunk,
                              mat_y,
                              prev_x_chunk,
                              prev_y,
                              n_samples,
                              config,
                              max_p_val = NULL,
                              upper_triangle = FALSE) {

  # --- 1. Intersection (Always Matrix) ---
  # BLAS optimized dense matrix mult (Cross-product)
  mat_int <- as.matrix(Matrix::tcrossprod(mat_x_chunk, mat_y))

  # Initialize Mask (Logical matrix of pairs to keep)
  keep_mask <- matrix(TRUE, nrow = nrow(mat_int), ncol = ncol(mat_int))

  # Filter 0: within-library comparison keeps one orientation per pair.
  # Applied first so that a chunk lying entirely on or below the diagonal
  # exits before any p-value is computed.
  if (upper_triangle) {
    keep_mask <- keep_mask & outer(idx_x_chunk, seq_len(ncol(mat_int)), "<")
    if (!any(keep_mask)) return(NULL)
  }

  # A threshold of 0 keeps everything, so skip the work
  thr <- function(name) {
    v <- config[[name]]
    if (!is.null(v) && v > 0) v else NULL
  }
  min_int <- thr("min_intersection")
  min_jacc <- thr("min_jaccard")
  min_over <- thr("min_overlap")

  # Filter A: Intersection
  if (!is.null(min_int)) {
    keep_mask <- keep_mask & (mat_int >= min_int)
    # Early exit if nothing passes
    if (!any(keep_mask)) return(NULL)
  }

  # --- 2. P-values (Matrix) ---
  # One-sided hypergeometric test for enrichment: P(intersection >= observed).
  # Computed for the whole block so that pairs can be pruned on p before the
  # (memory-hungry) pivot to long format.
  mat_p <- matrix(
    stats::phyper(
      q = as.vector(mat_int) - 1,
      m = rep(prev_x_chunk, times = ncol(mat_int)),
      n = n_samples - rep(prev_x_chunk, times = ncol(mat_int)),
      k = rep(prev_y, each = nrow(mat_int)),
      lower.tail = FALSE
    ),
    nrow = nrow(mat_int), ncol = ncol(mat_int)
  )

  if (!is.null(max_p_val)) {
    keep_mask <- keep_mask & (mat_p <= max_p_val)
    if (!any(keep_mask)) return(NULL)
  }

  # --- 3. Effect-Size Matrices ---
  # Only compute secondary matrices (Jaccard/Overlap) if they are used for FILTERING.
  if (!is.null(min_jacc)) {
    # Jaccard Denominator = (SizeX + SizeY) - Intersection
    mat_denom_jacc <- outer(prev_x_chunk, prev_y, "+") - mat_int
    mat_jacc <- safe_ratio(mat_int, mat_denom_jacc)
    # An NA ratio never passes the threshold (`NA >= t` would poison the mask)
    keep_mask <- keep_mask & !is.na(mat_jacc) & (mat_jacc >= min_jacc)
  }

  if (!is.null(min_over)) {
    # Overlap Denominator = min(SizeX, SizeY)
    mat_denom_over <- outer(prev_x_chunk, prev_y, pmin)
    mat_over <- safe_ratio(mat_int, mat_denom_over)
    keep_mask <- keep_mask & !is.na(mat_over) & (mat_over >= min_over)
  }

  if (!any(keep_mask)) return(NULL)

  # --- 4. The Pivot (Melt to Long) ---
  # Convert valid matrix positions to a list of pairs
  idx_pairs <- which(keep_mask, arr.ind = TRUE)

  # Extract values using array indexing
  v_int <- mat_int[idx_pairs]
  v_pval <- mat_p[idx_pairs]
  v_sz_x <- prev_x_chunk[idx_pairs[, 1]]
  v_sz_y <- prev_y[idx_pairs[, 2]]

  # --- 5. Compute Missing Stats (Vectorized) ---
  # Calculate final stats for the passing pairs
  v_union <- v_sz_x + v_sz_y - v_int

  if (!is.null(min_jacc)) {
    v_jacc <- mat_jacc[idx_pairs]
  } else {
    v_jacc <- safe_ratio(v_int, v_union)
  }

  if (!is.null(min_over)) {
    v_over <- mat_over[idx_pairs]
  } else {
    v_over <- safe_ratio(v_int, pmin(v_sz_x, v_sz_y))
  }

  # --- 6. Return ---
  tibble::tibble(
    idx_x = as.vector(idx_x_chunk[idx_pairs[, 1]]), # Absolute index in X
    idx_y = as.vector(idx_pairs[, 2]),              # Absolute index in Y (strip 'col' name)
    intersection = as.vector(v_int),
    size_x = as.vector(v_sz_x),
    size_y = as.vector(v_sz_y),
    union = as.vector(v_union),
    jaccard = as.vector(v_jacc),
    overlap = as.vector(v_over),
    p_val = as.vector(v_pval)
  )
}

#' Ratio With Zero Denominators Mapped to NA (Internal)
#'
#' Element-wise `num / den`, returning `NA_real_` where `den == 0` instead of
#' `NaN` (0/0) or `Inf` (x/0). Works on vectors and matrices; attributes of
#' `num` are preserved.
#'
#' @param num,den Numeric vectors or matrices of the same shape.
#' @return Numeric vector or matrix of the same shape as `num`.
#' @keywords internal
safe_ratio <- function(num, den) {
  out <- num / den
  out[den == 0] <- NA_real_
  out
}
