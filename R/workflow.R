#' Run Full Coreact Pipeline (Disc-to-Disc)
#'
#' @param path_x Path to input TSV X.
#' @param path_y Path to input TSV Y.
#' @param out_path Path for result file.
#' @param name_x Name for dataset X (default "Input_X").
#' @param name_y Name for dataset Y (default "Input_Y").
#' @param meta_cols_x Vector. Metadata columns for X.
#' @param meta_cols_y Vector. Metadata columns for Y.
#' @param feature_id_x Vector. Columns for Feature IDs in X.
#' @param feature_id_y Vector. Columns for Feature IDs in Y.
#' @param feature_id_sep String. Separator for Feature IDs.
#' @param min_prevalence Numeric vector (length 1 or 2).
#' @param filter_config List of structural filters.
#' @param fdr_threshold Numeric. Max FDR to retain (default 0.05).
#' @param n_cores Integer. Number of cores.
#' @param chunk_size Integer or NULL. Auto-calculated if NULL.
#' @export
coreact_pipeline <- function(path_x,
                             path_y,
                             out_path,
                             name_x = "Input_X",
                             name_y = "Input_Y",
                             meta_cols_x = 1,
                             meta_cols_y = 1,
                             feature_id_x = NULL,
                             feature_id_y = NULL,
                             feature_id_sep = "|",
                             min_prevalence = 0,
                             filter_config = list(min_intersection = 1),
                             fdr_threshold = 0.05,
                             n_cores = 1,
                             chunk_size = NULL) {

  # Handle Argument Recycling
  if (length(min_prevalence) == 1) min_prevalence <- rep(min_prevalence, 2)
  thresh_x <- min_prevalence[1]
  thresh_y <- min_prevalence[2]

  # --- Load Data ---
  obj_x <- import_coreact_tsv(path_x, meta_cols = meta_cols_x, feature_id = feature_id_x,
                              feature_id_sep = feature_id_sep, name = name_x, n_cores = n_cores)

  obj_y <- import_coreact_tsv(path_y, meta_cols = meta_cols_y, feature_id = feature_id_y,
                              feature_id_sep = feature_id_sep, name = name_y, n_cores = n_cores)

  # --- Check sample identity ---

  if (!all(colnames(obj_x$mat) == colnames(obj_y$mat))) stop("Sample identifiers in matrices X and Y are not identical.")

  # --- Pre-Filtering ---
  if (thresh_x > 0) {
    prev_x <- Matrix::rowSums(obj_x$mat > 0)
    keep_x <- prev_x >= thresh_x
    if (sum(keep_x) == 0) stop(sprintf("Filter removed all features from %s.", name_x))
    obj_x <- obj_x[keep_x, ]
  }

  if (thresh_y > 0) {
    prev_y <- Matrix::rowSums(obj_y$mat > 0)
    keep_y <- prev_y >= thresh_y
    if (sum(keep_y) == 0) stop(sprintf("Filter removed all features from %s.", name_y))
    obj_y <- obj_y[keep_y, ]
  }

  gc()

  # --- Compute ---
  results <- run_coreact(
    data_x = obj_x,
    data_y = obj_y,
    n_cores = n_cores,
    chunk_size = chunk_size,
    filter_config = filter_config,
    fdr_threshold = fdr_threshold
  )

  # --- Write ---
  if (nrow(results) == 0) warning("No results found. Writing empty file.")

  message(sprintf("Writing results to %s ...", out_path))

  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data.table::fwrite(results, file = out_path, sep = "\t", quote = FALSE)
  message("Pipeline completed successfully.")
}

#' Run Coreact Engine (Internal)
#' @export
run_coreact <- function(data_x,
                        data_y,
                        n_cores = 1,
                        chunk_size = NULL,
                        filter_config = list(min_intersection = 1),
                        fdr_threshold = NULL) {

  # 0. Extract Stats & Names
  # We extract these BEFORE swapping to ensure we track the names correctly.
  n_features_x <- nrow(data_x$mat)
  n_features_y <- nrow(data_y$mat)

  # Calculate prevalence
  prev_x <- Matrix::rowSums(data_x$mat > 0)
  prev_y <- Matrix::rowSums(data_y$mat > 0)

  # Get Row Names (Feature IDs)
  ids_current_x <- rownames(data_x$mat)
  ids_current_y <- rownames(data_y$mat)

  # Get N samples
  if (!all(colnames(data_x$mat) == colnames(data_y$mat))) {
    stop("Sample identifiers in matrices X and Y are not identical.")
  } else {
    n_samples <- ncol(data_x$mat) # Assumes X and Y have same samples
  }

  # 1. Swap Optimization
  swapped <- FALSE
  if (n_features_x < n_features_y) {
    message("Optimizing: Swapping X and Y for load balancing...")

    # Swap Data
    tmp_data <- data_x; data_x <- data_y; data_y <- tmp_data

    # Swap Stats
    tmp_prev <- prev_x; prev_x <- prev_y; prev_y <- tmp_prev
    n_features_x <- nrow(data_x$mat) # Update count

    # Swap Name References
    tmp_ids <- ids_current_x; ids_current_x <- ids_current_y; ids_current_y <- tmp_ids

    swapped <- TRUE
  }

  # 2. Auto-Chunking
  if (is.null(chunk_size)) {
    target_chunks <- max(n_cores * 20, 10)
    calc_size <- ceiling(n_features_x / target_chunks)
    chunk_size <- max(50, min(calc_size, 2000))
    if (n_features_x < 200) chunk_size <- 50
    message(sprintf("Auto-chunking: Size = %d rows (%d total chunks).",
                    chunk_size, ceiling(n_features_x/chunk_size)))
  }

  # 3. Execution
  chunks <- split(1:n_features_x, ceiling(seq_along(1:n_features_x) / chunk_size))

  run_fun <- function(idx) {
    # Note: We pass the SUBSET of prevalence for X, but ALL of prevalence for Y
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

  # 4. Bind & Process Results
  res <- res[!sapply(res, is.null)]
  if (length(res) == 0) return(tibble::tibble())

  final_df <- data.table::rbindlist(res) %>% tibble::as_tibble()

  # 5. Map Indices to IDs
  # idx_x from worker refers to ids_current_x
  # idx_y from worker refers to ids_current_y
  final_df$feature_x <- ids_current_x[final_df$idx_x]
  final_df$feature_y <- ids_current_y[final_df$idx_y]

  # 6. Un-Swap if needed
  if (swapped) {
    # If we swapped, "feature_x" in the df actually holds the Y IDs, and vice versa.
    # We rename columns to restore original meaning.
    final_df <- dplyr::rename(final_df,
                              idx_x_tmp = idx_x, idx_y_tmp = idx_y,
                              size_x_tmp = size_x, size_y_tmp = size_y,
                              feat_x_tmp = feature_x, feat_y_tmp = feature_y) %>%
      dplyr::rename(
        idx_x = idx_y_tmp, idx_y = idx_x_tmp,
        size_x = size_y_tmp, size_y = size_x_tmp,
        feature_x = feat_y_tmp, feature_y = feat_x_tmp
      )
  }

  # 7. FDR & Cleanup
  final_df$p_adj <- stats::p.adjust(final_df$p_val, method = "fdr")

  if (!is.null(fdr_threshold)) {
    message(sprintf("Filtering results by FDR < %s", fdr_threshold))
    final_df <- final_df[final_df$p_adj < fdr_threshold, ]
  }

  # Reorder columns for readability (IDs first)
  final_df <- final_df %>%
    dplyr::select(feature_x, feature_y, idx_x, idx_y, everything())

  return(final_df)
}
