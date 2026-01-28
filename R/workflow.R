#' Run Full Coreact Pipeline (Disc-to-Disc)
#'
#' @param path_x Path to input TSV X.
#' @param path_y Path to input TSV Y.
#' @param out_path Path for result file.
#' @param name_x Name for dataset X (default "Input_X").
#' @param name_y Name for dataset Y (default "Input_Y").
#' @param meta_cols_x Vector. Metadata columns for X.
#' @param meta_cols_y Vector. Metadata columns for Y.
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
                             min_prevalence = 0,
                             filter_config = list(min_intersection = 1),
                             fdr_threshold = 0.05,
                             n_cores = 1,
                             chunk_size = NULL) {

  # Handle Argument Recycling
  if (length(min_prevalence) == 1) min_prevalence <- rep(min_prevalence, 2)
  thresh_x <- min_prevalence[1]
  thresh_y <- min_prevalence[2]

  # --- Load & Filter ---
  obj_x <- import_coreact_tsv(path_x, meta_cols = meta_cols_x, name = name_x, n_cores = n_cores)
  obj_y <- import_coreact_tsv(path_y, meta_cols = meta_cols_y, name = name_y, n_cores = n_cores)

  if (thresh_x > 0) {
    obj_x <- filter_prevalence(obj_x, thresh_x)
    if (is.null(obj_x)) stop(sprintf("Filter removed all features from %s.", name_x))
  }

  if (thresh_y > 0) {
    obj_y <- filter_prevalence(obj_y, thresh_y)
    if (is.null(obj_y)) stop(sprintf("Filter removed all features from %s.", name_y))
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

  # 1. Swap Optimization
  swapped <- FALSE
  if (data_x$n_features < data_y$n_features) {
    message("Optimizing: Swapping X and Y for load balancing...")
    tmp <- data_x; data_x <- data_y; data_y <- tmp; swapped <- TRUE
  }

  # 2. Auto-Chunking
  if (is.null(chunk_size)) {
    target_chunks <- max(n_cores * 20, 10)
    calc_size <- ceiling(data_x$n_features / target_chunks)
    chunk_size <- max(50, min(calc_size, 2000))
    if (data_x$n_features < 200) chunk_size <- 50
    message(sprintf("Auto-chunking: Size = %d rows (%d total chunks).",
                    chunk_size, ceiling(data_x$n_features/chunk_size)))
  }

  # 3. Execution
  chunks <- split(1:data_x$n_features, ceiling(seq_along(1:data_x$n_features) / chunk_size))

  run_fun <- function(idx) {
    worker_chunk_calc(idx, data_x$mat[idx,,drop=FALSE], data_y$mat,
                      data_x$prevalence[idx], data_y$prevalence,
                      data_x$n_samples, filter_config)
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

  # 4. Bind & FDR
  res <- res[!sapply(res, is.null)]
  if (length(res) == 0) return(tibble::tibble())

  final_df <- data.table::rbindlist(res) %>% tibble::as_tibble()

  if (swapped) {
    final_df <- dplyr::rename(final_df,
                              idx_x_tmp = idx_x, idx_y_tmp = idx_y, size_x_tmp = size_x, size_y_tmp = size_y) %>%
      dplyr::rename(
        idx_x = idx_y_tmp, idx_y = idx_x_tmp, size_x = size_y_tmp, size_y = size_x_tmp)
  }

  final_df$p_adj <- stats::p.adjust(final_df$p_val, method = "fdr")

  if (!is.null(fdr_threshold)) {
    message(sprintf("Filtering results by FDR < %s", fdr_threshold))
    final_df <- final_df[final_df$p_adj < fdr_threshold, ]
  }

  return(final_df)
}
