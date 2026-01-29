#' Internal Worker: Compute Co-occurrence Stats for a Chunk
#'
#' @param idx_x_chunk Vector of absolute indices for the rows of X in this chunk.
#' @param mat_x_chunk The actual data sub-matrix for X (dense or sparse).
#' @param mat_y The full matrix Y (dense or sparse).
#' @param prev_x_chunk Vector of prevalence (rowSums) for rows in X chunk.
#' @param prev_y Vector of prevalence (rowSums) for all rows in Y.
#' @param n_samples Integer, total number of samples (ncol).
#' @param config List of name-value pairs of filter thresholds. Current valid options are: `min_intersection`, `min_jaccard`, `min_overlap`.
#'
#' @return A tibble of passing pairs with stats.
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

  if (calc_jacc_filter || calc_over_filter) {

    if (calc_jacc_filter) {
      # Jaccard Denominator = (SizeX + SizeY) - Intersection
      mat_denom_jacc <- outer(prev_x_chunk, prev_y, "+") - mat_int
      mat_jacc <- mat_int / mat_denom_jacc
      keep_mask <- keep_mask & (mat_jacc >= config$min_jaccard)
    }

    if (calc_over_filter) {
      # Overlap Denominator = min(SizeX, SizeY)
      mat_denom_over <- outer(prev_x_chunk, prev_y, pmin)
      mat_over <- mat_int / mat_denom_over
      keep_mask <- keep_mask & (mat_over >= config$min_overlap)
    }
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
    v_jacc <- v_int / v_union
  }

  if (calc_over_filter) {
    v_over <- mat_over[idx_pairs]
  } else {
    v_over <- v_int / pmin(v_sz_x, v_sz_y)
  }

  # --- 5. Compute P-values ---
  # Hypergeometric test
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
