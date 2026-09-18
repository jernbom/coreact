#' Run Full Coreact Pipeline (Disc-to-Disc)
#'
#' @param paths Character vector (length 2). Paths to input TSV files (X and Y).
#' @param out_path Path for result file.
#' @param names Character vector (length 2). Names for datasets X and Y.
#'    Default `c("Input_X", "Input_Y")`.
#' @param meta_cols List or Vector. Names or indices of metadata columns.
#'    If a vector, applied to both X and Y.
#'    If a list of length 2, `meta_cols[[1]]` is used for X, `meta_cols[[2]]` for Y.
#' @param sample_cols Optional List or Vector (default=`NULL`). Names or indices of sample columns.
#'    If a vector, applied to both X and Y.
#'    If a list of length 2, `sample_cols[[1]]` is used for X, `sample_cols[[2]]` for Y.
#'    If `NULL`, all columns not defined in `meta_cols` are used as sample columns.
#' @param feature_ids List or Vector. Names or indices of columns for constructing unique Feature IDs. If `NULL`, row numbers are used.
#'    If a vector, applied to both X and Y.
#'    If a list of length 2, `feature_ids[[1]]` is used for X, `feature_ids[[2]]` for Y.
#' @param feature_id_sep String. Separator for Feature IDs.
#' @param binarize Logical vector (length 1 or 2). Should the input data matrices (non-meta columns) be converted to binary?
#' @param min_prevalence Numeric vector (length 1 or 2). Minimum prevalence (rowSum)
#'    to retain a feature. If < 1, treated as percentage. If >= 1, treated as absolute count.
#' @param filter_config List of structural filters (e.g., `min_intersection`).
#' @param fdr_threshold Numeric. Max FDR to retain (default 0.05). If `NULL`, no filtering.
#' @param n_cores Integer. Number of cores.
#' @param chunk_size Integer or `NULL`. Auto-calculated if `NULL`.
#'
#' @details
#' P-values are adjusted with the Benjamini-Hochberg procedure
#' ([stats::p.adjust()] with `method = "BH"`, equivalent to `"fdr"`). Note that
#' the adjustment is applied to the set of pairs that survive `filter_config`
#' (by default, pairs with at least one shared sample), not to all possible
#' pairs. Use `filter_config = list(min_intersection = 0)` to adjust over the
#' complete pair space.
#'
#' @return Invisibly, the tibble of results that was written to `out_path`
#'   (including the `p_adj` column). Files are written as a side effect: the
#'   main results table at `out_path` and two metadata sidecars
#'   (`*_metadata_x.tsv`, `*_metadata_y.tsv`) next to it.
#' @export
coreact_pipeline <- function(paths,
                             out_path,
                             names = c("Input_X", "Input_Y"),
                             meta_cols = 1,
                             sample_cols = NULL,
                             feature_ids = NULL,
                             feature_id_sep = "|",
                             binarize = FALSE,
                             min_prevalence = 0,
                             filter_config = list(min_intersection = 1),
                             fdr_threshold = 0.05,
                             n_cores = 1,
                             chunk_size = NULL) {

  # --- 1. Argument Validation & Expansion ---
  if (!is.character(paths) || length(paths) != 2) {
    stop("Argument 'paths' must be a character vector of length 2.")
  }
  if (!is.character(names) || length(names) != 2) {
    stop("Argument 'names' must be a character vector of length 2.")
  }
  if (!is.character(out_path) || length(out_path) != 1) {
    stop("Argument 'out_path' must be a single file path.")
  }

  # Resolve args
  m_cols <- resolve_xy_args(meta_cols, "meta_cols")
  f_ids  <- resolve_xy_args(feature_ids, "feature_ids")

  # Resolve sample_cols if present
  s_cols <- if (!is.null(sample_cols)) resolve_xy_args(sample_cols, "sample_cols") else NULL

  # Binarize recycling
  if (!is.logical(binarize) || !length(binarize) %in% 1:2 || anyNA(binarize)) {
    stop("Argument 'binarize' must be a logical vector of length 1 or 2 without NA.")
  }
  if (length(binarize) == 1) binarize <- rep(binarize, 2)

  # Prevalence recycling
  if (!is.numeric(min_prevalence) || !length(min_prevalence) %in% 1:2 || anyNA(min_prevalence)) {
    stop("Argument 'min_prevalence' must be a numeric vector of length 1 or 2 without NA.")
  }
  if (length(min_prevalence) == 1) min_prevalence <- rep(min_prevalence, 2)

  # Fail fast on a malformed filter_config, before reading (potentially large) files
  validate_filter_config(filter_config)

  # --- 2. Load Data ---
  # Load X
  obj_x <- import_coreact_tsv(
    paths[1],
    meta_cols = m_cols$x,
    sample_cols = if (!is.null(s_cols)) s_cols$x else NULL,
    feature_id = f_ids$x,
    feature_id_sep = feature_id_sep,
    name = names[1],
    n_cores = n_cores,
    binarize = binarize[1]
  )

  # Load Y
  obj_y <- import_coreact_tsv(
    paths[2],
    meta_cols = m_cols$y,
    sample_cols = if (!is.null(s_cols)) s_cols$y else NULL,
    feature_id = f_ids$y,
    feature_id_sep = feature_id_sep,
    name = names[2],
    n_cores = n_cores,
    binarize = binarize[2]
  )

  # --- 3. Consistency Checks ---
  if (!identical(colnames(obj_x$mat), colnames(obj_y$mat))) {
    # Generate detailed error report
    stop(format_col_mismatch(obj_x$mat, obj_y$mat, names[1], names[2]))
  }

  # --- 4. Pre-Filtering (Rows/Features) ---
  obj_x <- filter_by_prevalence(obj_x, min_prevalence[1])
  obj_y <- filter_by_prevalence(obj_y, min_prevalence[2])

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
  # Benjamini-Hochberg over the pairs that survived the structural filters.
  # Always add the column so the output schema is identical for empty results.
  results$p_adj <- stats::p.adjust(results$p_val, method = "BH")

  if (!is.null(fdr_threshold) && nrow(results) > 0) {
    message(sprintf("Filtering results by FDR <= %s ...", fdr_threshold))
    results <- results[results$p_adj <= fdr_threshold, , drop = FALSE]
  }

  # Write Main Interactions File
  if (nrow(results) == 0) warning("No significant results found. Writing empty file.")

  message(sprintf("Writing %d results to %s ...", nrow(results), out_path))

  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data.table::fwrite(results, file = out_path, sep = "\t", quote = FALSE)

  # --- 7. Write Sidecar Metadata ---
  ids_x_final <- if (nrow(results) > 0) unique(results$feature_x) else character(0)
  ids_y_final <- if (nrow(results) > 0) unique(results$feature_y) else character(0)

  write_metadata_sidecar(
    obj = obj_x,
    keep_ids = ids_x_final,
    out_path = out_path,
    suffix = "x"
  )

  write_metadata_sidecar(
    obj = obj_y,
    keep_ids = ids_y_final,
    out_path = out_path,
    suffix = "y"
  )

  message("Pipeline completed successfully.")
  invisible(results)
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

  # --- 0. Input Validation ---
  if (!inherits(data_x, "coreact_data") || !inherits(data_y, "coreact_data")) {
    stop("'data_x' and 'data_y' must be 'coreact_data' objects (see new_coreact_data()).")
  }
  validate_filter_config(filter_config)
  if (!is.numeric(n_cores) || length(n_cores) != 1 || is.na(n_cores) || n_cores < 1) {
    stop("'n_cores' must be a single integer >= 1.")
  }
  if (!is.null(chunk_size) &&
        (!is.numeric(chunk_size) || length(chunk_size) != 1 || is.na(chunk_size) || chunk_size < 1)) {
    stop("'chunk_size' must be NULL or a single integer >= 1.")
  }

  # Samples must be the same set, in the same order, for tcrossprod to mean
  # "shared samples". The pipeline checks this on import; check again here so
  # that direct callers get the same protection.
  cn_x <- colnames(data_x$mat)
  cn_y <- colnames(data_y$mat)
  if (is.null(cn_x) || is.null(cn_y)) {
    if (ncol(data_x$mat) != ncol(data_y$mat)) {
      stop(sprintf("Sample count differs: '%s' has %d columns, '%s' has %d.",
                   data_x$name, ncol(data_x$mat), data_y$name, ncol(data_y$mat)))
    }
  } else if (!identical(cn_x, cn_y)) {
    stop(format_col_mismatch(data_x$mat, data_y$mat, data_x$name, data_y$name))
  }

  # --- 1. Swap Optimization ---
  # Chunk over the larger feature set so that the (smaller) Y matrix is what
  # every worker holds in full.
  swapped <- FALSE
  if (nrow(data_x$mat) < nrow(data_y$mat)) {
    message("Optimizing: Swapping X and Y for load balancing...")
    tmp <- data_x
    data_x <- data_y
    data_y <- tmp
    swapped <- TRUE
  }

  n_features_x <- nrow(data_x$mat)

  # Nothing to compare: return the empty schema rather than failing in split().
  if (n_features_x == 0 || nrow(data_y$mat) == 0) {
    return(empty_result_tibble())
  }

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
                    chunk_size, ceiling(n_features_x / chunk_size)))
  }
  chunk_size <- as.integer(chunk_size)

  # --- 3. Execution ---
  feature_idx <- seq_len(n_features_x)
  chunks <- split(feature_idx, ceiling(feature_idx / chunk_size))

  run_fun <- function(idx) {
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
  res <- Filter(Negate(is.null), res)

  # Prepare Empty Schema if no results
  if (length(res) == 0) {
    return(empty_result_tibble())
  }

  final_df <- data.table::rbindlist(res) %>% tibble::as_tibble()

  # --- 5. Map Indices to IDs ---
  final_df$feature_x <- rownames(data_x$mat)[final_df$idx_x]
  final_df$feature_y <- rownames(data_y$mat)[final_df$idx_y]

  # --- 6. Handle Un-Swapping ---
  if (swapped) {
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

#' Empty Result Schema (Internal)
#'
#' The zero-row tibble returned by [run_coreact()] when no pair passes the
#' structural filters or when an input has no features. Column names and
#' types match a non-empty result.
#'
#' @return A zero-row tibble.
#' @keywords internal
empty_result_tibble <- function() {
  tibble::tibble(
    feature_x = character(), feature_y = character(),
    idx_x = integer(), idx_y = integer(),
    intersection = numeric(), size_x = numeric(), size_y = numeric(),
    union = numeric(), jaccard = numeric(), overlap = numeric(), p_val = numeric()
  )
}

#' Validate Structural Filter Configuration (Internal)
#'
#' @param config A named list. Recognised entries are `min_intersection`,
#'   `min_jaccard` and `min_overlap`, each a single non-negative number.
#' @return `config`, invisibly. Errors on unknown names or invalid values.
#' @keywords internal
validate_filter_config <- function(config) {
  if (is.null(config)) return(invisible(config))
  if (!is.list(config)) stop("'filter_config' must be a named list.")
  if (length(config) == 0) return(invisible(config))

  allowed <- c("min_intersection", "min_jaccard", "min_overlap")
  nms <- names(config)
  if (is.null(nms) || any(nms == "")) {
    stop("All elements of 'filter_config' must be named.")
  }
  unknown <- setdiff(nms, allowed)
  if (length(unknown) > 0) {
    stop(sprintf("Unknown 'filter_config' entries: %s. Valid names are: %s.",
                 paste(unknown, collapse = ", "), paste(allowed, collapse = ", ")))
  }
  for (nm in nms) {
    v <- config[[nm]]
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0) {
      stop(sprintf("'filter_config$%s' must be a single non-negative number.", nm))
    }
  }
  invisible(config)
}
