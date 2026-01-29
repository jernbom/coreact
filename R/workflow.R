#' Run Full Coreact Pipeline (Disc-to-Disc)
#'
#' @param paths Character vector (length 2). Paths to input TSV files (X and Y).
#' @param out_path Path for result file.
#' @param names Character vector (length 2). Names for datasets X and Y.
#'   Default `c("Input_X", "Input_Y")`.
#' @param meta_cols List or Vector. Names or indices of metadata columns.
#'   If a vector, applied to both X and Y.
#'   If a list of length 2, `meta_cols[[1]]` is used for X, `meta_cols[[2]]` for Y.
#' @param feature_ids List or Vector. Names or indices of columns for constructing unique Feature IDs. If `NULL`, row numbers are used.
#'   If a vector, applied to both X and Y.
#'   If a list of length 2, `feature_ids[[1]]` is used for X, `feature_ids[[2]]` for Y.
#' @param feature_id_sep String. Separator for Feature IDs.
#' @param min_prevalence Numeric vector (length 1 or 2). Minimum prevalence (rowSum)
#'   to retain a feature.
#' @param filter_config List of structural filters (e.g., `min_intersection`).
#' @param fdr_threshold Numeric. Max FDR to retain (default 0.05). If `NULL`, no filtering.
#' @param n_cores Integer. Number of cores.
#' @param chunk_size Integer or `NULL`. Auto-calculated if `NULL`.
#' @export
coreact_pipeline <- function(paths,
                             out_path,
                             names = c("Input_X", "Input_Y"),
                             meta_cols = 1,
                             feature_ids = NULL,
                             feature_id_sep = "|",
                             min_prevalence = 0,
                             filter_config = list(min_intersection = 1),
                             fdr_threshold = 0.05,
                             n_cores = 1,
                             chunk_size = NULL) {

  # --- 1. Argument Validation & Expansion ---
  if (length(paths) != 2 && !is.character(paths)) stop("Argument 'paths' must be a character vector of length 2.")
  if (length(names) != 2 && !is.character(names)) stop("Argument 'names' must be a character vector of length 2.")

  # Helper: Resolve arguments that can be shared or specific (list vs vector)
  resolve_xy <- function(val, arg_name) {
    if (is.list(val)) {
      if (length(val) == 1) return(list(x = val[[1]], y = val[[1]]))
      if (length(val) == 2) return(list(x = val[[1]], y = val[[2]]))
      stop(sprintf("Argument '%s' provided as a list must be length 1 or 2.", arg_name))
    } else {
      # Atomic vector (or NULL): Recycle for both
      return(list(x = val, y = val))
    }
  }

  m_cols <- resolve_xy(meta_cols, "meta_cols")
  f_ids  <- resolve_xy(feature_ids, "feature_ids")

  # Prevalence recycling
  if (length(min_prevalence) > 2) stop("min_prevalence must have length 1 or 2.")
  if (length(min_prevalence) == 1) min_prevalence <- rep(min_prevalence, 2)

  # --- 2. Load Data ---
  # Load X
  obj_x <- import_coreact_tsv(
    paths[1],
    meta_cols = m_cols$x,
    feature_id = f_ids$x,
    feature_id_sep = feature_id_sep,
    name = names[1],
    n_cores = n_cores
  )

  # Load Y
  obj_y <- import_coreact_tsv(
    paths[2],
    meta_cols = m_cols$y,
    feature_id = f_ids$y,
    feature_id_sep = feature_id_sep,
    name = names[2],
    n_cores = n_cores
  )

  # --- 3. Consistency Checks ---
  if (!identical(colnames(obj_x$mat), colnames(obj_y$mat))) {
    stop("Error: Sample identifiers (columns) in matrices X and Y are not identical.")
  }

  # --- 4. Pre-Filtering ---
  # Helper to filter an object on prevalence
  prev_filter <- function(obj, thresh) {
    if (thresh <= 0) return(obj)
    prev <- Matrix::rowSums(obj$mat > 0)
    keep <- prev >= thresh
    if (sum(keep) == 0) stop(sprintf("Filter removed all features from %s.", obj$name))
    obj[keep, ] # Uses the S3 subsetting method
  }

  obj_x <- prev_filter(obj_x, min_prevalence[1])
  obj_y <- prev_filter(obj_y, min_prevalence[2])

  gc()

  # --- 5. Compute (Engine) ---
  results <- run_coreact(
    data_x = obj_x,
    data_y = obj_y,
    n_cores = n_cores,
    chunk_size = chunk_size,
    filter_config = filter_config
  )

  # --- 6. Post-Processing (FDR & Write) ---
  if (nrow(results) > 0) {
    # Calculate FDR
    results$p_adj <- stats::p.adjust(results$p_val, method = "fdr")

    # Filter by FDR
    if (!is.null(fdr_threshold)) {
      message(sprintf("Filtering results by FDR <= %s ...", fdr_threshold))
      results <- results[results$p_adj <= fdr_threshold, ]
    }
  }

  # Write
  if (nrow(results) == 0) warning("No significant results found. Writing empty file.")

  message(sprintf("Writing results to %s ...", out_path))

  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data.table::fwrite(results, file = out_path, sep = "\t", quote = FALSE)
  message("Pipeline completed successfully.")
}

#' Run Coreact Engine (Internal)
#'
#' @param data_x,data_y coreact_data objects.
#' @param n_cores Integer.
#' @param chunk_size Integer.
#' @param filter_config List of thresholds.
#'
#' @return A tibble of results.
#' @export
run_coreact <- function(data_x,
                        data_y,
                        n_cores = 1,
                        chunk_size = NULL,
                        filter_config = list(min_intersection = 1)) {

  # --- 1. Swap Optimization ---
  # Strategy: Chunk the LARGER dataset.
  # If X is small (10 rows) and Y is huge (1M rows), 1 chunk of X results in
  # a calculation of (10 x N) * (N x 1M) -> 10 x 1M dense matrix (Huge memory).
  # If we swap, we chunk the 1M rows. Each chunk is small, memory is managed.

  swapped <- FALSE
  if (nrow(data_x$mat) < nrow(data_y$mat)) {
    message("Optimizing: Swapping X and Y for load balancing...")
    tmp <- data_x; data_x <- data_y; data_y <- tmp
    swapped <- TRUE
  }

  n_features_x <- nrow(data_x$mat)

  # Calculate Stats
  prev_x <- Matrix::rowSums(data_x$mat > 0)
  prev_y <- Matrix::rowSums(data_y$mat > 0)
  n_samples <- ncol(data_x$mat)

  # --- 2. Auto-Chunking ---
  if (is.null(chunk_size)) {
    target_chunks <- max(n_cores * 20, 10)
    calc_size <- ceiling(n_features_x / target_chunks)
    chunk_size <- max(50, min(calc_size, 2000))
    if (n_features_x < 200) chunk_size <- 50
    message(sprintf("Auto-chunking: Size = %d rows (%d total chunks).",
                    chunk_size, ceiling(n_features_x/chunk_size)))
  }

  # --- 3. Execution ---
  chunks <- split(1:n_features_x, ceiling(seq_along(1:n_features_x) / chunk_size))

  run_fun <- function(idx) {
    # Note: data_x is the object being chunked.
    # idx refers to rows in the CURRENT data_x.
    worker_chunk_calc(
      idx_x_chunk = idx,
      mat_x_chunk = data_x$mat[idx, , drop=FALSE],
      mat_y = data_y$mat,
      prev_x_chunk = prev_x[idx],
      prev_y = prev_y,
      n_samples = n_samples,
      config = filter_config
    )
  }

  if (n_cores > 1) {
    if (.Platform$OS.type == "windows") {
      warning("Fork clusters not supported on Windows. Running serially.")
      res <- lapply(chunks, run_fun)
    } else {
      cl <- parallel::makeForkCluster(n_cores)
      on.exit(parallel::stopCluster(cl))
      res <- parallel::parLapply(cl, chunks, run_fun)
    }
  } else {
    res <- lapply(chunks, run_fun)
  }

  # --- 4. Assembly ---
  res <- res[!sapply(res, is.null)]
  if (length(res) == 0) return(tibble::tibble())

  final_df <- data.table::rbindlist(res) %>% tibble::as_tibble()

  # --- 5. Map Indices to IDs ---
  # At this moment:
  # 'idx_x' refers to rows of the CURRENT data_x
  # 'idx_y' refers to rows of the CURRENT data_y
  # The matrix rownames contain the Feature IDs.

  final_df$feature_x <- rownames(data_x$mat)[final_df$idx_x]
  final_df$feature_y <- rownames(data_y$mat)[final_df$idx_y]

  # --- 6. Handle Un-Swapping ---
  if (swapped) {
    # If we swapped:
    # data_x holds the Original Y data.
    # data_y holds the Original X data.
    # Therefore:
    # column 'feature_x' actually contains IDs from Original Y.
    # column 'feature_y' actually contains IDs from Original X.
    # We must rename them to restore the user's expected order (X, Y).

    final_df <- dplyr::rename(final_df,
                              idx_x_tmp = idx_x, idx_y_tmp = idx_y,
                              size_x_tmp = size_x, size_y_tmp = size_y,
                              feat_x_tmp = feature_x, feat_y_tmp = feature_y
    ) %>%
      dplyr::rename(
        # Map Y (temp X) back to Y
        idx_y = idx_x_tmp,
        size_y = size_x_tmp,
        feature_y = feat_x_tmp,

        # Map X (temp Y) back to X
        idx_x = idx_y_tmp,
        size_x = size_y_tmp,
        feature_x = feat_y_tmp
      )
  }

  # Reorder columns for readability
  final_df <- final_df %>%
    dplyr::select(feature_x, feature_y, idx_x, idx_y, everything())

  return(final_df)
}
