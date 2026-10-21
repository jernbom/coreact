#' Run the Outcome Association From Files (Disc-to-Disc)
#'
#' The disk-to-disk counterpart of [run_coreact_outcome()], as
#' [coreact_pipeline()] is of [run_coreact()]. Reads the results that
#' [coreact_pipeline()] wrote, rebuilds the inputs they were computed on from
#' the provenance sidecar beside them, reads the outcome from a TSV, runs the
#' association and writes the table.
#'
#' @param results_path String. The `out_path` given to [coreact_pipeline()].
#'   The provenance sidecar written next to it (`*_provenance.tsv`) is
#'   required: it records how to rebuild the inputs.
#' @param outcome_path String. A TSV with one row per sample: a sample-ID column
#'   and the outcome column, and any other columns, which are ignored.
#' @param out_path String. Where to write the result. A provenance sidecar
#'   recording the settings used is written next to it.
#' @param outcome_col Name or index of the outcome column in `outcome_path`.
#' @param sample_id_col Name or index of the sample-ID column in `outcome_path`.
#'   Default: the first column. IDs are read as text, so leading zeros survive.
#' @param input_paths `NULL`, or the input files' current paths when they have
#'   moved since [coreact_pipeline()] read them: X then Y, or one path for a
#'   within-library result. By default the paths recorded in the provenance are
#'   used. Either way the files must hold the data that was analysed, which is
#'   checked.
#' @inheritParams run_coreact_outcome
#'
#' @details
#' The inputs are re-imported from the paths recorded in the provenance, with
#' the recorded import settings, and filtered with the recorded prevalence
#' bounds, which reproduces the objects [coreact_pipeline()] analysed. That is
#' checked rather than assumed: each object is compared with its recorded
#' fingerprint — name, size, positive cells and a hash of its contents — and
#' every pair's indices with its feature IDs. So a file edited since, or the
#' wrong file given in `input_paths`, is refused. One recorded input path means
#' a within-library result.
#'
#' The outcome file is one row per sample, transposed relative to the X and Y
#' files, because a sample annotation mixes types and [data.table::fread()]
#' types each column: so `0`/`1` and `TRUE`/`FALSE` columns stay binary and
#' numeric ones stay numeric. A text column is refused rather than guessed at,
#' since which of two words is the event cannot be read off the data. Its sample
#' IDs must be unique and must match the X and Y sample columns exactly; empty
#' outcome cells are missing values.
#'
#' @return Invisibly, the tibble that was written to `out_path`, as returned by
#'   [run_coreact_outcome()]. Files written: the table at `out_path`, and
#'   `*_provenance.tsv` next to it, recording the input files, the column
#'   choices and the settings, including the resolved `near_nested_max` that
#'   sets the size of the tested family.
#' @export
#' @examples
#' path_x <- system.file("extdata", "input_x.tsv", package = "coreact")
#' path_y <- system.file("extdata", "input_y.tsv", package = "coreact")
#' results <- file.path(tempdir(), "cooccurrence.tsv")
#' coreact_pipeline(
#'   paths = c(path_x, path_y),
#'   out_path = results,
#'   meta_cols = list(c("u_pep_id", "pep_nr", "protein_id", "gene_name"),
#'                    c("protein_id", "gene_name")),
#'   feature_ids = list("u_pep_id", c("protein_id", "gene_name")),
#'   fdr_threshold = 1
#' )
#'
#' # One row per sample: an ID column and the outcome
#' samples <- setdiff(names(data.table::fread(path_x, nrows = 0)),
#'                    c("u_pep_id", "pep_nr", "protein_id", "gene_name"))
#' set.seed(1)
#' outcome_file <- file.path(tempdir(), "outcome.tsv")
#' utils::write.table(data.frame(sample = samples, case = stats::rbinom(length(samples), 1, 0.5)),
#'                    outcome_file, sep = "\t", row.names = FALSE, quote = FALSE)
#'
#' out <- coreact_outcome_pipeline(results, outcome_file,
#'                                 out_path = file.path(tempdir(), "association.tsv"),
#'                                 outcome_col = "case")
#' head(out[out$role == "tested", c("feature_x", "feature_y", "contrast", "beta", "p_adj")])
coreact_outcome_pipeline <- function(results_path,
                                     outcome_path,
                                     out_path,
                                     outcome_col,
                                     sample_id_col = 1,
                                     input_paths = NULL,
                                     outcome_type = c("binary", "continuous"),
                                     near_nested_max = 10,
                                     contrasts = c("interaction", "11_vs_10", "11_vs_01", "11_vs_00",
                                                   "10_vs_00", "01_vs_00", "10_vs_01"),
                                     ci = c("tested", "all", "none"),
                                     fdr_method = "BH",
                                     n_cores = 1,
                                     chunk_size = NULL) {

  # --- 1. Arguments, before any file is read ---
  for (arg in c("results_path", "outcome_path", "out_path")) {
    v <- get(arg)
    if (!is.character(v) || length(v) != 1 || is.na(v)) stop(sprintf("'%s' must be a single file path.", arg))
  }
  check_column_arg(outcome_col, "outcome_col")
  check_column_arg(sample_id_col, "sample_id_col")
  if (!is.null(input_paths)) {
    if (!is.character(input_paths) || !length(input_paths) %in% 1:2 || anyNA(input_paths)) {
      stop("'input_paths' must be NULL or one or two file paths.")
    }
    for (p in input_paths) if (!file.exists(p)) stop(sprintf("File not found: %s", p))
  }
  outcome_type <- match.arg(outcome_type)
  ci <- match.arg(ci)
  contrasts <- validate_contrasts(contrasts)
  validate_fdr_method(fdr_method)
  validate_parallel_args(n_cores, chunk_size)
  if (!is.function(near_nested_max)) check_near_nested_max(near_nested_max)

  if (!file.exists(results_path)) stop(sprintf("File not found: %s", results_path))
  sidecar <- sidecar_path(results_path, "provenance")
  if (!file.exists(sidecar)) {
    stop(sprintf(paste0("No provenance sidecar at %s, so the inputs of these results cannot be ",
                        "rebuilt. It is written by coreact_pipeline() from coreact 1.5.0 on."),
                 sidecar))
  }
  if (!file.exists(outcome_path)) stop(sprintf("File not found: %s", outcome_path))
  if (identical(normalizePath(out_path, winslash = "/", mustWork = FALSE),
                normalizePath(results_path, winslash = "/"))) {
    stop("'out_path' is the results file itself; choose another path.")
  }

  # --- 2. The results and the inputs they were computed on ---
  res <- read_coreact_results(results_path)
  prov <- attr(res, "provenance")
  n_inputs <- if (isTRUE(prov$within_library)) 1 else 2
  if (!is.null(input_paths) && length(input_paths) != n_inputs) {
    stop(sprintf("These results were computed from %d input(s), but 'input_paths' gives %d.",
                 n_inputs, length(input_paths)))
  }
  inputs <- rebuild_inputs(prov, n_cores = n_cores, input_paths = input_paths)

  # Check the rebuilt objects here, where the message can speak of files
  remedy <- if (is.null(input_paths)) {
    "The file has changed since coreact_pipeline() read it."
  } else {
    "Check 'input_paths': they must be the files coreact_pipeline() read, X first."
  }
  if (!is.null(inputs$y) && identical(inputs$x$mat, inputs$y$mat)) {
    stop("The two inputs hold the same data. ", remedy)
  }
  for (side in c("x", "y")[seq_len(n_inputs)]) {
    check_fingerprint(inputs[[side]], prov, side, remedy = remedy,
                      what = sprintf("The %s input (%s)", toupper(side), inputs$paths[[side]]))
  }

  # --- 3. The outcome, one row per sample ---
  outcome <- read_outcome_tsv(outcome_path, outcome_col, sample_id_col, n_cores = n_cores)

  # --- 4. The association ---
  out <- run_coreact_outcome(res, inputs$x, inputs$y, outcome = outcome,
                             outcome_type = outcome_type, near_nested_max = near_nested_max,
                             contrasts = contrasts, ci = ci, fdr_method = fdr_method,
                             n_cores = n_cores, chunk_size = chunk_size)

  # --- 5. Write the table and what it was computed with ---
  message(sprintf("Writing %d rows to %s ...", nrow(out), out_path))
  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  data.table::fwrite(out, file = out_path, sep = "\t", quote = FALSE)
  write_provenance_sidecar(
    outcome_provenance(out, results_path, outcome_path, outcome_col, sample_id_col,
                       input_paths = unname(inputs$paths)),
    out_path
  )

  message("Outcome pipeline completed successfully.")
  invisible(out)
}

#' Validate a Column Chosen by Name or Index (Internal)
#'
#' @param col The argument's value.
#' @param arg Its name, for the error message.
#' @return `NULL`, invisibly.
#' @keywords internal
check_column_arg <- function(col, arg) {
  ok <- length(col) == 1 && !is.na(col) &&
    (is.character(col) || (is.numeric(col) && col >= 1 && col == round(col)))
  if (!ok) stop(sprintf("'%s' must be a single column name or a column index.", arg))
  invisible(NULL)
}

#' Rebuild the Inputs of a Coreact Result From Its Provenance (Internal)
#'
#' Re-imports each recorded input file with its recorded import settings and
#' filters it with its recorded prevalence bounds, as [coreact_pipeline()] did.
#'
#' @param prov Provenance list of a [coreact_pipeline()] result.
#' @param n_cores Threads for reading.
#' @param input_paths `NULL`, or paths to read instead of the recorded ones,
#'   X then Y.
#' @return List with `x` and `y` (`NULL` for a within-library result), and
#'   `paths`, the files read.
#' @keywords internal
rebuild_inputs <- function(prov, n_cores = 1, input_paths = NULL) {
  fields <- c("_path", "_input_name", "_meta_cols", "_sample_cols", "_feature_ids", "_binarize",
              "_min_prevalence", "_max_prevalence")
  needed <- c("within_library", "feature_id_sep", paste0("x", fields))
  if (!isTRUE(prov$within_library)) needed <- c(needed, paste0("y", fields))
  absent <- setdiff(needed, names(prov))
  if (length(absent) > 0) {
    stop(sprintf(paste0("The provenance does not record how its inputs were imported (%s missing); ",
                        "only results written by coreact_pipeline() can be rebuilt from disk."),
                 paste(utils::head(absent, 3), collapse = ", ")))
  }

  paths <- c(x = prov$x_path, y = if (isTRUE(prov$within_library)) NULL else prov$y_path)
  if (!is.null(input_paths)) paths[] <- input_paths

  side <- function(tag) {
    get_field <- function(field) prov[[paste0(tag, "_", field)]]
    path <- paths[[tag]]
    if (!file.exists(path)) {
      stop(sprintf(paste0("The %s input these results were computed from is no longer at %s. ",
                          "If it has moved, give its new path in 'input_paths'."),
                   toupper(tag), path))
    }
    obj <- import_coreact_tsv(path, meta_cols = get_field("meta_cols"),
                              sample_cols = get_field("sample_cols"),
                              feature_id = get_field("feature_ids"),
                              feature_id_sep = prov$feature_id_sep,
                              name = get_field("input_name"), n_cores = n_cores,
                              binarize = get_field("binarize"))
    filter_by_prevalence(obj, get_field("min_prevalence"), get_field("max_prevalence"))
  }

  list(x = side("x"), y = if (isTRUE(prov$within_library)) NULL else side("y"), paths = paths)
}

#' Read an Outcome From a One-Row-per-Sample TSV (Internal)
#'
#' @param path String.
#' @param outcome_col,sample_id_col Name or index of each column.
#' @param n_cores Threads for reading.
#' @return The outcome vector, named by sample ID.
#' @keywords internal
read_outcome_tsv <- function(path, outcome_col, sample_id_col, n_cores = 1) {
  header <- names(data.table::fread(path, sep = "\t", header = TRUE, nrows = 0))
  resolve <- function(col, arg) {
    if (is.character(col)) {
      if (!col %in% header) {
        stop(sprintf("'%s': column '%s' not found in %s.", arg, col, basename(path)))
      }
      return(col)
    }
    if (col > length(header)) {
      stop(sprintf("'%s' = %d, but %s has %d columns.", arg, as.integer(col), basename(path), length(header)))
    }
    header[col]
  }
  id_name <- resolve(sample_id_col, "sample_id_col")
  outcome_name <- resolve(outcome_col, "outcome_col")
  if (identical(id_name, outcome_name)) {
    stop("'outcome_col' and 'sample_id_col' name the same column.")
  }

  # Sample IDs as text, or "007" would come back as 7
  dt <- data.table::fread(path, sep = "\t", header = TRUE, select = c(id_name, outcome_name),
                          colClasses = list(character = id_name), nThread = n_cores)
  data.table::setDF(dt)
  ids <- dt[[id_name]]
  values <- dt[[outcome_name]]
  if (anyNA(ids) || any(ids == "")) {
    stop(sprintf("The sample-ID column '%s' of %s has empty entries.", id_name, basename(path)))
  }
  if (is.character(values)) {
    stop(sprintf(paste0("The outcome column '%s' of %s holds text. Code a binary outcome as 0/1 or ",
                        "TRUE/FALSE, and a continuous one as numbers."),
                 outcome_name, basename(path)))
  }
  stats::setNames(values, ids)
}

#' Record How an Outcome Table Was Computed (Internal)
#'
#' @param out The table from [run_coreact_outcome()].
#' @param results_path,outcome_path,outcome_col,sample_id_col As given to
#'   [coreact_outcome_pipeline()].
#' @param input_paths The input files actually read, X then Y.
#' @return A flat named list for [write_provenance_sidecar()].
#' @keywords internal
outcome_provenance <- function(out, results_path, outcome_path, outcome_col, sample_id_col,
                               input_paths) {
  list(
    coreact_version = as.character(utils::packageVersion("coreact")),
    results_path = normalizePath(results_path, winslash = "/"),
    outcome_path = normalizePath(outcome_path, winslash = "/"),
    input_paths = normalizePath(input_paths, winslash = "/"),
    outcome_col = outcome_col,
    sample_id_col = sample_id_col,
    outcome_type = attr(out, "outcome_type"),
    near_nested_max = as.numeric(attr(out, "near_nested_max")),
    contrasts = attr(out, "contrasts"),
    ci = attr(out, "ci"),
    fdr_method = attr(out, "fdr_method"),
    within_library = attr(out, "within_library"),
    n_used = if (nrow(out) > 0) out$n_used[1] else NA_integer_,
    frac_unbounded_tested = as.numeric(attr(out, "frac_unbounded_tested")),
    n_pairs_below_selection = attr(out, "n_pairs_below_selection")
  )
}
