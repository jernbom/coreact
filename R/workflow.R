#' Run Full Coreact Pipeline (Disc-to-Disc)
#'
#' Reads one or two TSV files and runs the three-step co-reactivity procedure:
#' feature selection on marginal prevalence, significance testing of every
#' retained feature pair with FDR control, and effect-size filtering of the
#' significant pairs. Results and metadata sidecars are written to disk.
#'
#' @param paths Character vector (length 1 or 2). Paths to input TSV files.
#'    Two paths compare X with Y; one path compares a library with itself, and
#'    every other per-dataset argument must then be a single value. See
#'    Within-library comparison below.
#' @param out_path Path for result file.
#' @param names Character vector. Names for the datasets: length 2 for two
#'    inputs (default `c("Input_X", "Input_Y")`), length 1 for one input
#'    (default `"Input"`).
#' @param meta_cols List or Vector. Names or indices of metadata columns.
#'    If a vector, applied to both X and Y.
#'    If a list of length 2, `meta_cols[[1]]` is used for X, `meta_cols[[2]]` for Y.
#' @param sample_cols Optional List or Vector (default=`NULL`). Names or indices of sample columns.
#'    If a vector, applied to both X and Y.
#'    If a list of length 2, `sample_cols[[1]]` is used for X, `sample_cols[[2]]` for Y.
#'    If `NULL`, all columns not defined in `meta_cols` are used as sample columns.
#' @param feature_ids List or Vector. Names or indices of columns for constructing unique
#'    Feature IDs. If `NULL`, row numbers are used.
#'    If a vector, applied to both X and Y.
#'    If a list of length 2, `feature_ids[[1]]` is used for X, `feature_ids[[2]]` for Y.
#' @param feature_id_sep String. Separator for Feature IDs.
#' @param binarize Logical vector (length 1 or 2). Should the input data matrices (non-meta
#'    columns) be converted to binary?
#' @param min_prevalence Numeric vector (length 1 or 2). Minimum prevalence (number of
#'    positive samples) to retain a feature. If < 1, treated as a fraction of the samples.
#'    If >= 1, treated as an absolute count. `0` disables the filter.
#' @param max_prevalence Numeric vector (length 1 or 2) or `NULL`. Maximum prevalence to
#'    retain a feature, with the same fraction/count convention as `min_prevalence`.
#'    `NULL` (default) disables the filter.
#' @param fdr_threshold Numeric. Maximum Benjamini-Hochberg adjusted p-value to retain
#'    (default 0.05). If `NULL`, every pair is retained and only the effect-size
#'    filters apply; note that this returns all `n_x * n_y` pairs.
#' @param effect_filters Named list of minimum effect sizes applied to the significant
#'    pairs. Recognised entries: `min_intersection`, `min_jaccard`, `min_overlap`.
#'    The default applies no effect-size filter.
#' @param n_cores Integer. Number of cores.
#' @param chunk_size Integer or `NULL`. Auto-calculated if `NULL`.
#'
#' @details
#' The procedure has three steps, in this order:
#'
#' 1. **Feature selection.** Features are retained by their marginal prevalence
#'    (`min_prevalence`, `max_prevalence`). The hypergeometric test conditions
#'    on the prevalences, so this selection does not affect the validity of the
#'    p-values of the retained pairs.
#' 2. **Significance testing.** A one-sided hypergeometric p-value is computed
#'    for every pair of retained features, adjusted with the Benjamini-Hochberg
#'    procedure over all `n_x * n_y` pairs, and pairs with
#'    `p_adj <= fdr_threshold` are kept.
#' 3. **Effect-size filtering.** The significant pairs are filtered by
#'    `effect_filters`.
#'
#' Steps 2 and 3 are implemented so that pairs which cannot end up in the
#' output are discarded inside the engine before they are materialised: a
#' pair with a raw p-value above `fdr_threshold` can never have an adjusted
#' p-value below it, and the effect-size thresholds commute with the FDR
#' filter. The output is therefore identical to computing all pairs first,
#' but memory and disk usage scale with the number of reported pairs rather
#' than with `n_x * n_y`. Effect-size filters are a convenience for keeping
#' output sizes manageable, and with a moderate `fdr_threshold` they are
#' usually unnecessary.
#'
#' @section Within-library comparison:
#' A single `paths` entry compares the features of one library with each other.
#' Each unordered pair is reported once and no feature is paired with itself, so
#' step 2 adjusts over `n * (n - 1) / 2` pairs rather than `n^2`; see
#' [run_coreact()] for why. With one input, `names`, `meta_cols`,
#' `sample_cols`, `feature_ids`, `binarize`, `min_prevalence` and
#' `max_prevalence` each take a single value.
#'
#' Passing the same file twice is refused, whether as the same path or as a copy
#' under another name that imports to the same matrix.
#'
#' @return Invisibly, the tibble of results that was written to `out_path`
#'   (including the `p_adj` column), with the attributes `n_pairs` and
#'   `provenance` described in [run_coreact()]. The pipeline's `provenance` also
#'   records `fdr_threshold`, `feature_id_sep` and, per side, the input path
#'   (made absolute), name, `meta_cols`, `sample_cols`, `feature_ids`,
#'   `binarize` and prevalence bounds as given.
#'
#'   Files are written as a side effect: the main results table at `out_path`,
#'   two metadata sidecars (`*_metadata_x.tsv`, `*_metadata_y.tsv`) and a
#'   provenance sidecar (`*_provenance.tsv`) next to it. Within a library the
#'   two metadata sidecars hold the features appearing as `feature_x` and as
#'   `feature_y` respectively, which differ at the ends of the index order.
#'   [read_coreact_results()] reads the results back together with their
#'   provenance, so that [adjust_fdr()] can be applied to them again.
#' @export
#' @examples
#' path_x <- system.file("extdata", "input_x.tsv", package = "coreact")
#' path_y <- system.file("extdata", "input_y.tsv", package = "coreact")
#' out <- file.path(tempdir(), "interactions.tsv")
#'
#' res <- coreact_pipeline(
#'   paths = c(path_x, path_y),
#'   out_path = out,
#'   meta_cols = list(c("u_pep_id", "pep_nr", "protein_id", "gene_name"),
#'                    c("protein_id", "gene_name")),
#'   feature_ids = list("u_pep_id", c("protein_id", "gene_name")),
#'   min_prevalence = c(0.1, 2),
#'   fdr_threshold = 1,
#'   effect_filters = list(min_overlap = 0.5)
#' )
#' head(res)
#' list.files(tempdir(), pattern = "^interactions")
#'
#' # Within one library: a single path, single-valued arguments
#' out_within <- file.path(tempdir(), "within.tsv")
#' res_within <- coreact_pipeline(
#'   paths = path_x,
#'   out_path = out_within,
#'   meta_cols = c("u_pep_id", "pep_nr", "protein_id", "gene_name"),
#'   feature_ids = "u_pep_id",
#'   fdr_threshold = 1
#' )
#' all(res_within$idx_x < res_within$idx_y)
coreact_pipeline <- function(paths,
                             out_path,
                             names = c("Input_X", "Input_Y"),
                             meta_cols = 1,
                             sample_cols = NULL,
                             feature_ids = NULL,
                             feature_id_sep = "|",
                             binarize = FALSE,
                             min_prevalence = 0,
                             max_prevalence = NULL,
                             fdr_threshold = 0.05,
                             effect_filters = list(),
                             n_cores = 1,
                             chunk_size = NULL) {

  # --- 1. Argument Validation & Expansion ---
  if (!is.character(paths) || !length(paths) %in% 1:2) {
    stop("Argument 'paths' must be a character vector of length 1 or 2.")
  }
  # One path declares a within-library comparison; see run_coreact().
  within_library <- length(paths) == 1

  # Two paths naming the same file would test every pair twice and pair every
  # feature with itself. Compare resolved paths before reading anything.
  if (!within_library && file.exists(paths[1]) && file.exists(paths[2]) &&
        identical(normalizePath(paths[1]), normalizePath(paths[2]))) {
    stop("Both entries of 'paths' resolve to the same file. For a ",
         "within-library comparison, pass that path alone.")
  }

  # The length-2 default cannot stand for a single input.
  if (missing(names) && within_library) names <- "Input"
  if (!is.character(names) || length(names) != length(paths)) {
    stop(sprintf("Argument 'names' must be a character vector of length %d, matching 'paths'.",
                 length(paths)))
  }
  if (within_library) {
    # How "one per dataset" is written differs by argument, as resolve_xy_args()
    # defines it: meta_cols / sample_cols / feature_ids take a LIST of two,
    # since an atomic vector of length two names two columns of one dataset,
    # while binarize and the prevalence bounds take a vector of two.
    as_list <- list(meta_cols = meta_cols, sample_cols = sample_cols,
                    feature_ids = feature_ids)
    as_vec <- list(binarize = binarize, min_prevalence = min_prevalence,
                   max_prevalence = max_prevalence)
    two_valued <- c(
      vapply(as_list, function(v) is.list(v) && length(v) > 1, logical(1)),
      vapply(as_vec, function(v) !is.null(v) && length(v) > 1, logical(1))
    )
    if (any(two_valued)) {
      stop(sprintf(
        paste0("With one entry in 'paths' there is one dataset, so %s must be a single ",
               "value, not one per dataset."),
        paste(sprintf("'%s'", base::names(two_valued)[two_valued]), collapse = ", ")
      ))
    }
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

  if (is.null(max_prevalence)) {
    max_prevalence <- c(NA_real_, NA_real_)
  } else {
    if (!is.numeric(max_prevalence) || !length(max_prevalence) %in% 1:2 || anyNA(max_prevalence)) {
      stop("Argument 'max_prevalence' must be NULL or a numeric vector of length 1 or 2 without NA.")
    }
    if (length(max_prevalence) == 1) max_prevalence <- rep(max_prevalence, 2)
  }

  if (!is.null(fdr_threshold) &&
        (!is.numeric(fdr_threshold) || length(fdr_threshold) != 1 || is.na(fdr_threshold) ||
           fdr_threshold < 0 || fdr_threshold > 1)) {
    stop("Argument 'fdr_threshold' must be NULL or a single number in [0, 1].")
  }

  # Fail fast on malformed effect filters, before reading (potentially large) files
  validate_effect_filters(effect_filters)

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

  # Load Y, unless this is a within-library comparison
  obj_y <- if (within_library) NULL else import_coreact_tsv(
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
  # One input is trivially consistent with itself.
  if (!within_library) {
    if (!identical(colnames(obj_x$mat), colnames(obj_y$mat))) {
      # Generate detailed error report
      stop(format_col_mismatch(obj_x$mat, obj_y$mat, names[1], names[2]))
    }
    # Two different paths can still hold the same data: a copy under another
    # name imports to the same matrix, and would test every pair twice.
    if (identical(obj_x$mat, obj_y$mat)) {
      stop("The two entries of 'paths' hold the same data. For a within-library ",
           "comparison, pass one path alone; each pair is then reported once.")
    }
  }

  # --- 4. Step 1: Feature Selection (Marginal Prevalence) ---
  obj_x <- filter_by_prevalence(obj_x, min_prevalence[1],
                                if (is.na(max_prevalence[1])) NULL else max_prevalence[1])
  if (!within_library) {
    obj_y <- filter_by_prevalence(obj_y, min_prevalence[2],
                                  if (is.na(max_prevalence[2])) NULL else max_prevalence[2])
  }

  gc()

  # --- 5. Step 2 + 3: Significance Testing and Effect-Size Filtering (Engine) ---
  # The engine prunes on raw p-value (exact, see ?coreact_pipeline) and on
  # effect size; the FDR adjustment below is over ALL pairs of retained features.
  results <- run_coreact(
    data_x = obj_x,
    data_y = obj_y,
    n_cores = n_cores,
    chunk_size = chunk_size,
    effect_filters = effect_filters,
    max_p_val = fdr_threshold
  )

  # Add the pipeline's own arguments to the engine's provenance: with them the
  # inputs can be re-imported and re-filtered into the very objects the engine
  # ran on. Each side's values are recorded as resolved, so a per-dataset list
  # is stored per side and never has to be re-interpreted.
  import_side <- function(i, tag) {
    pick <- function(resolved) if (is.null(resolved)) NULL else resolved[[c("x", "y")[i]]]
    stats::setNames(
      list(normalizePath(paths[i], winslash = "/"), names[i], pick(m_cols), pick(s_cols),
           pick(f_ids), binarize[i], min_prevalence[i],
           if (is.na(max_prevalence[i])) NULL else max_prevalence[i]),
      paste0(tag, c("_path", "_input_name", "_meta_cols", "_sample_cols", "_feature_ids",
                    "_binarize", "_min_prevalence", "_max_prevalence"))
    )
  }
  attr(results, "provenance") <- c(
    attr(results, "provenance"),
    list(fdr_threshold = fdr_threshold, feature_id_sep = feature_id_sep),
    import_side(1, "x"),
    if (!within_library) import_side(2, "y")
  )

  results <- adjust_fdr(results)

  if (!is.null(fdr_threshold)) {
    message(sprintf("Filtering results by FDR <= %s ...", fdr_threshold))
    results <- results[results$p_adj <= fdr_threshold, , drop = FALSE]
  }

  # Write Main Interactions File
  if (nrow(results) == 0) warning("No significant results found. Writing empty file.")

  message(sprintf("Writing %d results to %s ...", nrow(results), out_path))

  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data.table::fwrite(results, file = out_path, sep = "\t", quote = FALSE)

  # --- 6. Write Sidecar Metadata ---
  ids_x_final <- if (nrow(results) > 0) unique(results$feature_x) else character(0)
  ids_y_final <- if (nrow(results) > 0) unique(results$feature_y) else character(0)

  write_metadata_sidecar(
    obj = obj_x,
    keep_ids = ids_x_final,
    out_path = out_path,
    suffix = "x"
  )

  # Within a library both sidecars come from the one object. They are not
  # identical: under the kept triangle the last feature in index order never
  # appears as feature_x, and the first never as feature_y.
  write_metadata_sidecar(
    obj = if (within_library) obj_x else obj_y,
    keep_ids = ids_y_final,
    out_path = out_path,
    suffix = "y"
  )

  # --- 7. Write Provenance ---
  # Restores what fwrite() cannot carry: the family size adjust_fdr() needs,
  # and the settings that reproduce the filtered inputs.
  write_provenance_sidecar(attr(results, "provenance"), out_path)

  message("Pipeline completed successfully.")
  invisible(results)
}

#' Run Coreact Engine on In-Memory Objects
#'
#' The computational core of the package: computes co-occurrence statistics
#' for every pair of features between two `coreact_data` objects that share
#' the same samples, or between the features of a single object.
#' [coreact_pipeline()] wraps this with file import,
#' feature selection, FDR adjustment and export; call `run_coreact()`
#' directly when the data are already in memory, then [adjust_fdr()] to add
#' adjusted p-values.
#'
#' @param data_x,data_y `coreact_data` objects (see [new_coreact_data()] and
#'   [import_coreact_tsv()]). Their sample columns must be identical and in
#'   the same order. Omit `data_y` to compare a library with itself; see
#'   Within-library comparison below. Supplying the same data as both
#'   arguments is an error.
#' @param n_cores Integer. Number of worker processes. Values above 1 use a
#'   fork cluster and are therefore only effective on Linux/macOS; on Windows
#'   the computation falls back to serial with a warning.
#' @param chunk_size Integer or `NULL`. Number of X features per work chunk.
#'   Auto-calculated from the data size and `n_cores` if `NULL`.
#' @param effect_filters Named list of minimum effect sizes. Recognised
#'   entries: `min_intersection` (minimum number of shared samples),
#'   `min_jaccard`, `min_overlap`. A missing entry or `0` disables that
#'   filter; the default `list()` keeps every pair.
#' @param max_p_val Numeric or `NULL`. Pairs with a raw p-value above this
#'   are not returned. Set it to the FDR threshold you intend to apply after
#'   [adjust_fdr()]: since an adjusted p-value is never below the raw one,
#'   the set of significant pairs is unchanged and memory is saved.
#'
#' @details
#' For every pair `(x, y)` the following are reported: `intersection`
#' (samples positive for both), `size_x` / `size_y` (prevalence of each
#' feature), `union`, `jaccard` (`intersection / union`), `overlap`
#' (`intersection / min(size_x, size_y)`) and `p_val`, the one-sided
#' hypergeometric probability of observing at least `intersection` shared
#' samples given the two prevalences and the total number of samples
#' (`phyper(intersection - 1, size_x, n_samples - size_x, size_y,
#' lower.tail = FALSE)`). Jaccard and overlap are `NA` when their
#' denominator is zero (features absent from every sample).
#'
#' The larger of the two feature sets is chunked and processed in parallel;
#' if `data_y` is larger the two are swapped internally and the result columns
#' are mapped back, so the output always refers to the inputs as supplied.
#'
#' @section Within-library comparison:
#' Omitting `data_y` compares the features of `data_x` with each other. Each
#' unordered pair is then reported once, as `idx_x < idx_y`, and no feature is
#' paired with itself: the family is `n * (n - 1) / 2` rather than `n^2`. Which
#' feature of a pair is `x` and which is `y` follows index order and carries no
#' other meaning.
#'
#' Reporting both orientations would enter the same hypothesis in the family
#' twice, and a feature paired with itself co-occurs perfectly by construction,
#' so the self-pairs would occupy the top of the ranking without being
#' findings. Both inflate the family that [adjust_fdr()] corrects over.
#'
#' Supplying one library as both `data_x` and `data_y` is refused, because the
#' duplicate pairs cannot be removed reliably once two objects are involved. A
#' copy of the same data under a different name cannot always be recognised,
#' so when the two inputs share features that are identical across all samples
#' a warning says how many and what it costs.
#'
#' @return A tibble with one row per pair passing `effect_filters` and
#'   `max_p_val`, with columns `feature_x`, `feature_y`, `idx_x`, `idx_y`,
#'   `intersection`, `size_x`, `size_y`, `union`, `jaccard`, `overlap`,
#'   `p_val`. `idx_*` are row indices into the corresponding input object.
#'   The attribute `n_pairs` holds the number of tested pairs
#'   (`nrow(data_x) * nrow(data_y)`, or `n * (n - 1) / 2` when `data_y` is
#'   omitted), which [adjust_fdr()] uses as the size
#'   of the family. No multiple-testing adjustment is applied here.
#'
#'   The attribute `provenance` is a named list recording how the result was
#'   computed: the coreact version, `within_library`, `n_samples`, `n_pairs`,
#'   `max_p_val` and the effect-size thresholds (`effect_min_*`, `0` meaning
#'   off). For each side (`x_*`, and `y_*` unless within a library) it records
#'   the object's `name`, `n_features` and `n_positive` (positive cells), as
#'   received, and the prevalence cutoffs it was filtered with as counts.
#'   Cutoffs are `NA` when the object does not carry them, for example because
#'   it was built directly rather than filtered by [coreact_pipeline()].
#' @export
#' @examples
#' path_x <- system.file("extdata", "input_x.tsv", package = "coreact")
#' path_y <- system.file("extdata", "input_y.tsv", package = "coreact")
#'
#' x <- import_coreact_tsv(path_x, meta_cols = 1:4, feature_id = "u_pep_id", name = "X")
#' y <- import_coreact_tsv(path_y, meta_cols = 1:2, feature_id = c("protein_id", "gene_name"),
#'                         name = "Y")
#'
#' res <- run_coreact(x, y)
#' attr(res, "n_pairs")
#' res <- adjust_fdr(res)
#' head(res[order(res$p_adj), ])
#'
#' # Within a single library: each pair once, no feature against itself
#' res_within <- run_coreact(x)
#' attr(res_within, "n_pairs")           # nrow(x) * (nrow(x) - 1) / 2
#' all(res_within$idx_x < res_within$idx_y)
run_coreact <- function(data_x,
                        data_y = NULL,
                        n_cores = 1,
                        chunk_size = NULL,
                        effect_filters = list(),
                        max_p_val = NULL) {

  # --- 0. Input Validation ---
  if (!inherits(data_x, "coreact_data")) {
    stop("'data_x' must be a 'coreact_data' object (see new_coreact_data()).")
  }

  # Omitting data_y declares a within-library comparison. The mode is never
  # inferred from the data: supplying one library as both arguments is refused,
  # because each unordered pair would then be tested twice (once in each
  # orientation) and every feature would also be paired with itself.
  within_library <- is.null(data_y)
  if (within_library) {
    data_y <- data_x
  } else {
    if (!inherits(data_y, "coreact_data")) {
      stop("'data_y' must be a 'coreact_data' object, or omitted for a ",
           "within-library comparison (see new_coreact_data()).")
    }
    if (identical(data_x$mat, data_y$mat)) {
      stop("'data_x' and 'data_y' hold the same data. For a within-library ",
           "comparison, supply 'data_x' only; each pair is then reported once.")
    }
    warn_shared_identical_features(data_x, data_y)
  }
  validate_effect_filters(effect_filters)
  if (!is.numeric(n_cores) || length(n_cores) != 1 || is.na(n_cores) || n_cores < 1) {
    stop("'n_cores' must be a single integer >= 1.")
  }
  if (!is.null(chunk_size) &&
        (!is.numeric(chunk_size) || length(chunk_size) != 1 || is.na(chunk_size) || chunk_size < 1)) {
    stop("'chunk_size' must be NULL or a single integer >= 1.")
  }
  if (!is.null(max_p_val) &&
        (!is.numeric(max_p_val) || length(max_p_val) != 1 || is.na(max_p_val) ||
           max_p_val < 0 || max_p_val > 1)) {
    stop("'max_p_val' must be NULL or a single number in [0, 1].")
  }

  # Samples must be the same set, in the same order, for tcrossprod to mean
  # "shared samples". The pipeline checks this on import; check again here so
  # that direct callers get the same protection.
  # One input is trivially consistent with itself.
  if (!within_library) {
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
  }

  # The family of tests. Within a library each unordered pair is tested once
  # and no feature is paired with itself, so the family is n * (n - 1) / 2.
  n_pairs <- if (within_library) {
    n_feat <- as.numeric(nrow(data_x$mat))
    n_feat * (n_feat - 1) / 2
  } else {
    as.numeric(nrow(data_x$mat)) * as.numeric(nrow(data_y$mat))
  }

  # Provenance is recorded here, BEFORE the swap below: after it, data_x may be
  # the caller's Y, and a record taken then would describe each side under the
  # other's name.
  provenance <- record_provenance(data_x, if (within_library) NULL else data_y,
                                  n_pairs = n_pairs, effect_filters = effect_filters,
                                  max_p_val = max_p_val)
  finish <- function(res) {
    attr(res, "provenance") <- provenance
    with_n_pairs(res, n_pairs)
  }

  # --- 1. Swap Optimization ---
  # Chunk over the larger feature set so that the (smaller) Y matrix is what
  # every worker holds in full. Never within a library, where the two sides are
  # one matrix and the kept triangle is defined on its index order.
  swapped <- FALSE
  if (!within_library && nrow(data_x$mat) < nrow(data_y$mat)) {
    message("Optimizing: Swapping X and Y for load balancing...")
    tmp <- data_x
    data_x <- data_y
    data_y <- tmp
    swapped <- TRUE
  }

  n_features_x <- nrow(data_x$mat)

  # Nothing to compare: return the empty schema rather than failing in split().
  if (n_features_x == 0 || nrow(data_y$mat) == 0) {
    return(finish(empty_result_tibble()))
  }

  # Calculate Stats (as double, so column types do not depend on Matrix version)
  prev_x <- as.numeric(Matrix::rowSums(data_x$mat > 0))
  prev_y <- as.numeric(Matrix::rowSums(data_y$mat > 0))
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
      mat_x_chunk = data_x$mat[idx, , drop = FALSE],
      mat_y = data_y$mat,
      prev_x_chunk = prev_x[idx],
      prev_y = prev_y,
      n_samples = n_samples,
      config = effect_filters,
      max_p_val = max_p_val,
      upper_triangle = within_library
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
    return(finish(empty_result_tibble()))
  }

  final_df <- tibble::as_tibble(data.table::rbindlist(res))

  # --- 5. Map Indices to IDs ---
  final_df$feature_x <- rownames(data_x$mat)[final_df$idx_x]
  final_df$feature_y <- rownames(data_y$mat)[final_df$idx_y]

  # --- 6. Handle Un-Swapping ---
  # Everything computed above refers to the swapped roles; swap the x/y
  # suffixed column names back in one pass.
  if (swapped) {
    swap_map <- c(idx_x = "idx_y", idx_y = "idx_x",
                  size_x = "size_y", size_y = "size_x",
                  feature_x = "feature_y", feature_y = "feature_x")
    nms <- names(final_df)
    hit <- nms %in% names(swap_map)
    nms[hit] <- swap_map[nms[hit]]
    names(final_df) <- nms
  }

  # Reorder columns for readability (same order as empty_result_tibble())
  finish(final_df[, names(empty_result_tibble())])
}

#' Adjust P-values Over the Full Pair Family
#'
#' Adds a `p_adj` column to a [run_coreact()] result using the
#' Benjamini-Hochberg procedure (or another [stats::p.adjust()] method) with
#' the number of tested pairs as the size of the family, so that pairs
#' pruned by `max_p_val` or `effect_filters` still count as tests.
#'
#' @param res A tibble returned by [run_coreact()].
#' @param n_pairs Number of tested pairs. Defaults to the `n_pairs` attribute
#'   set by [run_coreact()]; supply it explicitly if the attribute was lost.
#'   [read_coreact_results()] restores it for results read from disk.
#' @param method Adjustment method passed to [stats::p.adjust()].
#'
#' @details
#' Pairs absent from `res` are treated as having p-value 1. For `"BH"` this
#' reproduces the adjusted p-values of the complete computation for every
#' pair whose raw p-value is at most the `max_p_val` that was used, and it is
#' conservative for the rest.
#'
#' @return `res` with a `p_adj` column (replaced if already present).
#' @export
#' @examples
#' mat <- matrix(c(1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0), nrow = 3, byrow = TRUE,
#'               dimnames = list(c("A", "B", "C"), paste0("S", 1:4)))
#' obj <- new_coreact_data(mat, data.frame(row.names = rownames(mat)))
#' res <- adjust_fdr(run_coreact(obj))
#' res[, c("feature_x", "feature_y", "p_val", "p_adj")]
adjust_fdr <- function(res, n_pairs = attr(res, "n_pairs"), method = "BH") {
  if (!is.data.frame(res) || !"p_val" %in% names(res)) {
    stop("'res' must be a run_coreact() result with a 'p_val' column.")
  }
  if (is.null(n_pairs)) {
    stop("'n_pairs' is missing: 'res' has no 'n_pairs' attribute. ",
         "Pass the number of tested pairs explicitly.")
  }
  if (!is.numeric(n_pairs) || length(n_pairs) != 1 || is.na(n_pairs) || n_pairs < nrow(res)) {
    stop("'n_pairs' must be a single number no smaller than nrow(res).")
  }
  res$p_adj <- stats::p.adjust(res$p_val, method = method, n = n_pairs)
  res
}

#' Record How a Result Was Computed (Internal)
#'
#' Builds the provenance that [run_coreact()] attaches to its result. It is a
#' flat named list so that it writes to, and reads back from, a key-value
#' sidecar without restructuring.
#'
#' @param data_x,data_y `coreact_data` objects exactly as `run_coreact()`
#'   received them; `data_y = NULL` for a within-library comparison.
#' @param n_pairs Size of the tested family.
#' @param effect_filters,max_p_val As passed to `run_coreact()`.
#'
#' @details
#' Each side gets a fingerprint — object name, number of features, number of
#' positive cells — that later steps can check an object against, and the
#' prevalence cutoffs it was filtered with, as counts. Cutoffs are `NA` when
#' the object does not carry them (it was never filtered, or was subset or
#' combined since), which records them as unknown rather than as no filter.
#' Effect-size filters are recorded under all three names, `0` meaning off,
#' which is how [run_coreact()] treats a missing entry.
#'
#' @return A named list.
#' @keywords internal
record_provenance <- function(data_x, data_y, n_pairs, effect_filters, max_p_val) {
  side <- function(obj, tag) {
    cut <- attr(obj, "prevalence_cutoffs")
    if (is.null(cut)) cut <- c(min = NA_real_, max = NA_real_)
    stats::setNames(
      list(obj$name, as.numeric(nrow(obj$mat)), as.numeric(sum(obj$mat)),
           unname(cut[["min"]]), unname(cut[["max"]])),
      paste0(tag, c("_name", "_n_features", "_n_positive",
                    "_min_prevalence_count", "_max_prevalence_count"))
    )
  }
  effect <- function(name) {
    v <- effect_filters[[name]]
    if (is.null(v)) 0 else as.numeric(v)
  }

  c(
    list(
      coreact_version = as.character(utils::packageVersion("coreact")),
      within_library = is.null(data_y),
      n_samples = as.numeric(ncol(data_x$mat)),
      n_pairs = as.numeric(n_pairs),
      max_p_val = max_p_val,
      effect_min_intersection = effect("min_intersection"),
      effect_min_jaccard = effect("min_jaccard"),
      effect_min_overlap = effect("min_overlap")
    ),
    side(data_x, "x"),
    if (!is.null(data_y)) side(data_y, "y")
  )
}

#' Attach the Family Size to a Result (Internal)
#'
#' @param res Result tibble.
#' @param n_pairs Number of tested pairs.
#' @return `res` with the `n_pairs` attribute set.
#' @keywords internal
with_n_pairs <- function(res, n_pairs) {
  attr(res, "n_pairs") <- n_pairs
  res
}

#' Empty Result Schema (Internal)
#'
#' The zero-row tibble returned by [run_coreact()] when no pair passes the
#' filters or when an input has no features. Column names and types match a
#' non-empty result.
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

#' Validate Effect-Size Filter Configuration (Internal)
#'
#' @param config A named list. Recognised entries are `min_intersection`,
#'   `min_jaccard` and `min_overlap`, each a single non-negative number.
#' @return `config`, invisibly. Errors on unknown names or invalid values.
#' @keywords internal
validate_effect_filters <- function(config) {
  if (is.null(config)) return(invisible(config))
  if (!is.list(config)) stop("'effect_filters' must be a named list.")
  if (length(config) == 0) return(invisible(config))

  allowed <- c("min_intersection", "min_jaccard", "min_overlap")
  nms <- names(config)
  if (is.null(nms) || any(nms == "")) {
    stop("All elements of 'effect_filters' must be named.")
  }
  unknown <- setdiff(nms, allowed)
  if (length(unknown) > 0) {
    stop(sprintf("Unknown 'effect_filters' entries: %s. Valid names are: %s.",
                 paste(unknown, collapse = ", "), paste(allowed, collapse = ", ")))
  }
  for (nm in nms) {
    v <- config[[nm]]
    if (!is.numeric(v) || length(v) != 1 || is.na(v) || v < 0) {
      stop(sprintf("'effect_filters$%s' must be a single non-negative number.", nm))
    }
  }
  invisible(config)
}
