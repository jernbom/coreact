#' Internal Worker: Compute Co-occurrence Stats for a Chunk
#'
#' Computes intersection, union, Jaccard index, overlap coefficient and a
#' one-sided hypergeometric p-value for every (row of X chunk) x (row of Y)
#' pair, then keeps the pairs that pass the structural filters in `config`.
#'
#' @param idx_x_chunk Vector of absolute indices for the rows of X in this chunk.
#' @param mat_x_chunk The actual data sub-matrix for X (dense or sparse).
#' @param mat_y The full matrix Y (dense or sparse).
#' @param prev_x_chunk Vector of prevalence (rowSums) for rows in X chunk.
#' @param prev_y Vector of prevalence (rowSums) for all rows in Y.
#' @param n_samples Integer, total number of samples (ncol).
#' @param config List of name-value pairs of filter thresholds. Current valid
#'   options are: `min_intersection`, `min_jaccard`, `min_overlap`.
#'
#' @details
#' The Jaccard index is undefined when both features are absent from every
#' sample (union of size zero), and the overlap coefficient is undefined when
#' either feature is absent from every sample. Such ratios are reported as
#' `NA` and never pass a `min_jaccard` / `min_overlap` threshold. With the
#' default `min_intersection = 1` these pairs are excluded anyway, because
#' their intersection is zero.
#'
#' @return A tibble of passing pairs with stats, or `NULL` if none pass.
#' @keywords internal
worker_chunk_calc <- function(idx_x_chunk,
                              mat_x_chunk,
                              mat_y,
                              prev_x_chunk,
                              prev_y,
                              n_samples,
                              config) {

  # --- 1. Intersection (Always Matrix) ---
  # BLAS optimized dense matrix mult (Cross-product)
  mat_int <- as.matrix(Matrix::tcrossprod(mat_x_chunk, mat_y))

  # Initialize Mask (Logical matrix of pairs to keep)
  keep_mask <- matrix(TRUE, nrow = nrow(mat_int), ncol = ncol(mat_int))

  # Filter A: Intersection
  if (!is.null(config$min_intersection)) {
    keep_mask <- keep_mask & (mat_int >= config$min_intersection)
    # Early exit if nothing passes
    if (!any(keep_mask)) return(NULL)
  }

  # --- 2. Adaptive Matrix Calculation ---
  # Only compute secondary matrices (Jaccard/Overlap) if they are used for FILTERING.
  calc_jacc_filter <- !is.null(config$min_jaccard)
  calc_over_filter <- !is.null(config$min_overlap)

  if (calc_jacc_filter) {
    # Jaccard Denominator = (SizeX + SizeY) - Intersection
    mat_denom_jacc <- outer(prev_x_chunk, prev_y, "+") - mat_int
    mat_jacc <- safe_ratio(mat_int, mat_denom_jacc)
    # An NA ratio never passes the threshold (`NA >= t` would poison the mask)
    keep_mask <- keep_mask & !is.na(mat_jacc) & (mat_jacc >= config$min_jaccard)
  }

  if (calc_over_filter) {
    # Overlap Denominator = min(SizeX, SizeY)
    mat_denom_over <- outer(prev_x_chunk, prev_y, pmin)
    mat_over <- safe_ratio(mat_int, mat_denom_over)
    keep_mask <- keep_mask & !is.na(mat_over) & (mat_over >= config$min_overlap)
  }

  if (!any(keep_mask)) return(NULL)

  # --- 3. The Pivot (Melt to Long) ---
  # Convert valid matrix positions to a list of pairs
  idx_pairs <- which(keep_mask, arr.ind = TRUE)

  # Extract values using array indexing
  v_int <- mat_int[idx_pairs]
  v_sz_x <- prev_x_chunk[idx_pairs[, 1]]
  v_sz_y <- prev_y[idx_pairs[, 2]]

  # --- 4. Compute Missing Stats (Vectorized) ---
  # Calculate final stats for the passing pairs
  v_union <- v_sz_x + v_sz_y - v_int

  if (calc_jacc_filter) {
    v_jacc <- mat_jacc[idx_pairs]
  } else {
    v_jacc <- safe_ratio(v_int, v_union)
  }

  if (calc_over_filter) {
    v_over <- mat_over[idx_pairs]
  } else {
    v_over <- safe_ratio(v_int, pmin(v_sz_x, v_sz_y))
  }

  # --- 5. Compute P-values ---
  # One-sided hypergeometric test for enrichment: P(intersection >= observed)
  v_pval <- stats::phyper(
    q = v_int - 1,
    m = v_sz_x,
    n = n_samples - v_sz_x,
    k = v_sz_y,
    lower.tail = FALSE
  )

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
