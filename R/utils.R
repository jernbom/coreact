#' Resolve X/Y Pipeline Arguments (Internal)
#'
#' Handles the logic for arguments that can be provided as a single vector (applied to both)
#' or a list of length 2 (specific to X and Y).
#'
#' @param val The input argument (list or vector).
#' @param arg_name String. Name of the argument for error reporting.
#' @return A list with elements `$x` and `$y`.
#' @keywords internal
resolve_xy_args <- function(val, arg_name) {
  if (is.list(val)) {
    if (length(val) == 1) return(list(x = val[[1]], y = val[[1]]))
    if (length(val) == 2) return(list(x = val[[1]], y = val[[2]]))
    stop(sprintf("Argument '%s' provided as a list must be length 1 or 2.", arg_name))
  } else {
    # If atomic vector (or NULL), recycle the same value for both X and Y
    list(x = val, y = val)
  }
}

#' Filter Coreact Data by Marginal Prevalence (Internal)
#'
#' Retains the features of a coreact_data object whose prevalence (number of
#' positive samples) lies within `[min_thresh, max_thresh]`. Each threshold is
#' a fraction of the samples if `< 1` and an absolute count if `>= 1`.
#'
#' @param obj A `coreact_data` object.
#' @param min_thresh Numeric. Minimum prevalence; `0` disables the lower bound.
#' @param max_thresh Numeric or `NULL`. Maximum prevalence; `NULL` disables
#'   the upper bound.
#' @return A subsetted `coreact_data` object carrying the attribute
#'   `prevalence_cutoffs`: the bounds applied, resolved to counts, as
#'   `c(min = , max = )`, with `max = Inf` when there is no upper bound. Fractions
#'   are resolved with `ceiling()` for the lower bound and `floor()` for the
#'   upper, so the counts reproduce the filter exactly on the same samples.
#'
#'   The attribute describes this object only. Any later subsetting or `cbind()`
#'   builds a new object without it, and [run_coreact()] then records the
#'   cutoffs as unknown rather than as no filter.
#' @keywords internal
filter_by_prevalence <- function(obj, min_thresh = 0, max_thresh = NULL) {
  if (!is.numeric(min_thresh) || length(min_thresh) != 1 || is.na(min_thresh)) {
    stop("Minimum prevalence threshold must be a single non-NA number.")
  }
  if (!is.null(max_thresh) &&
        (!is.numeric(max_thresh) || length(max_thresh) != 1 || is.na(max_thresh))) {
    stop("Maximum prevalence threshold must be NULL or a single non-NA number.")
  }
  if (min_thresh <= 0 && is.null(max_thresh)) {
    return(with_prevalence_cutoffs(obj, 0, Inf))
  }

  n_samp <- ncol(obj$mat)

  # Lower cutoff: fraction -> smallest count reaching it; count -> as is
  lo <- if (min_thresh < 1) ceiling(min_thresh * n_samp) else min_thresh
  # Upper cutoff: fraction -> largest count not exceeding it; count -> as is
  hi <- if (is.null(max_thresh)) Inf else if (max_thresh < 1) floor(max_thresh * n_samp) else max_thresh

  if (hi < lo) {
    stop(sprintf("[%s] Prevalence bounds are empty: min %s > max %s samples.",
                 obj$name, format(lo), format(hi)))
  }

  prev <- Matrix::rowSums(obj$mat > 0)
  keep <- prev >= lo & prev <= hi
  n_keep <- sum(keep)

  bounds <- if (is.finite(hi)) sprintf("%s to %s", format(lo), format(hi)) else sprintf(">= %s", format(lo))

  if (n_keep == 0) {
    stop(sprintf("[%s] Filter removed all features (prevalence %s samples).", obj$name, bounds))
  }

  message(sprintf("[%s] Prevalence filter (%s samples): keeping %d of %d features.",
                  obj$name, bounds, n_keep, length(keep)))

  # Uses the S3 subsetting method '[.coreact_data' defined in class_coreact.R.
  # That method builds a fresh object, so the cutoffs are attached afterwards.
  with_prevalence_cutoffs(obj[keep, ], lo, hi)
}

#' Attach Resolved Prevalence Cutoffs (Internal)
#'
#' @param obj A `coreact_data` object.
#' @param lo,hi Lower and upper bounds as counts; `hi = Inf` for no upper bound.
#' @return `obj` with the `prevalence_cutoffs` attribute set.
#' @keywords internal
with_prevalence_cutoffs <- function(obj, lo, hi) {
  attr(obj, "prevalence_cutoffs") <- c(min = as.numeric(lo), max = as.numeric(hi))
  obj
}

#' Path of a Sidecar File Next to a Results File (Internal)
#'
#' `results.tsv` with tag `metadata_x` gives `results_metadata_x.tsv`. An
#' `out_path` without an extension gets the tag and `.tsv` appended.
#'
#' @param out_path String. The path of the main results file.
#' @param tag String. What the sidecar holds, e.g. `"metadata_x"`.
#' @return String.
#' @keywords internal
sidecar_path <- function(out_path, tag) {
  ext <- tools::file_ext(out_path)
  if (nzchar(ext)) {
    sprintf("%s_%s.%s", tools::file_path_sans_ext(out_path), tag, ext)
  } else {
    sprintf("%s_%s.tsv", out_path, tag)
  }
}

#' Write Sidecar Metadata (Internal)
#'
#' Writes a metadata TSV for features present in the final results.
#' Appends `_metadata_\{suffix\}` to the base filename.
#'
#' @param obj A `coreact_data` object (containing `$meta` and `$mat`).
#' @param keep_ids Character vector. Feature IDs appearing in the final interaction list.
#' @param out_path String. The path of the main results file.
#' @param suffix String. Suffix to distinguish X and Y (e.g., "x" or "y").
#' @return NULL (Writes file to disk).
#' @keywords internal
write_metadata_sidecar <- function(obj, keep_ids, out_path, suffix) {
  # Construct filename: results.tsv -> results_metadata_x.tsv
  meta_file <- sidecar_path(out_path, paste0("metadata_", suffix))

  # Extract IDs from the object (assuming rownames of matrix are the master IDs)
  all_ids <- rownames(obj$mat)

  # Check which features are in the final results
  to_keep <- all_ids %in% keep_ids

  # Create the export table: Feature ID first
  export_df <- data.frame(feature_id = all_ids[to_keep], stringsAsFactors = FALSE)

  # Bind metadata columns if they exist
  if (!is.null(obj$meta) && ncol(obj$meta) > 0) {
    # Note: obj$meta is expected to be row-aligned with obj$mat
    meta_subset <- obj$meta[to_keep, , drop = FALSE]
    rownames(meta_subset) <- NULL
    export_df <- cbind(export_df, meta_subset)
  }

  message(sprintf("Writing %s metadata (%d features) to %s ...", suffix, nrow(export_df), meta_file))

  # Ensure directory exists (redundant check but safe)
  out_dir <- dirname(meta_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  data.table::fwrite(export_df, file = meta_file, sep = "\t", quote = FALSE)
}

#' Write the Provenance Sidecar (Internal)
#'
#' Writes a result's provenance next to it as `<out>_provenance.tsv`, a
#' key-value table read back by [read_provenance_sidecar()].
#'
#' @details
#' The columns are `key`, `type`, `index` and `value`. A scalar is one row; a
#' vector is one row per element with `index` counting from 1; `NULL` and
#' empty vectors are a single row with `index` 0. `type` is one of
#' `character`, `integer`, `double`, `logical` or `null`, so that every value
#' reads back with its type. That matters for the import arguments: column 2
#' given by index and a column named `"2"` are different requests, and only
#' the type tells them apart once both are text.
#'
#' Doubles are written with the fewest significant digits that read back to
#' the identical number, so values such as `0.1` stay legible while `1/3`
#' round-trips exactly. A missing character value is refused rather than
#' written, because it could not be told apart from the string `"NA"`.
#'
#' @param provenance A flat named list of atomic vectors or `NULL`s.
#' @param out_path String. The path of the main results file.
#' @return The sidecar's path, invisibly.
#' @keywords internal
write_provenance_sidecar <- function(provenance, out_path) {
  file <- sidecar_path(out_path, "provenance")

  rows <- lapply(names(provenance), function(key) {
    v <- provenance[[key]]
    if (is.null(v)) {
      return(data.frame(key = key, type = "null", index = 0L, value = ""))
    }
    type <- typeof(v)
    if (!type %in% c("character", "integer", "double", "logical")) {
      stop(sprintf("Provenance entry '%s' has unsupported type '%s'.", key, type))
    }
    if (type == "character" && anyNA(v)) {
      stop(sprintf("Provenance entry '%s' contains a missing string.", key))
    }
    if (length(v) == 0) {
      return(data.frame(key = key, type = type, index = 0L, value = ""))
    }
    data.frame(key = key, type = type, index = seq_along(v), value = encode_value(v))
  })
  table <- do.call(rbind, rows)

  message(sprintf("Writing provenance (%d entries) to %s ...", length(provenance), file))
  out_dir <- dirname(file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  data.table::fwrite(table, file = file, sep = "\t")
  invisible(file)
}

#' Encode Atomic Values as Exact Text (Internal)
#'
#' @param v An atomic vector.
#' @return A character vector of the same length.
#' @keywords internal
encode_value <- function(v) {
  if (!is.double(v)) return(as.character(v))
  out <- trimws(formatC(v, digits = 15, format = "g"))
  # 15 significant digits are not always enough to identify a double; fall
  # back to 17, which always is. Only finite values need checking: the
  # tokens for NA and the infinities are exact by construction.
  fin <- is.finite(v)
  inexact <- fin
  inexact[fin] <- as.numeric(out[fin]) != v[fin]
  out[inexact] <- trimws(formatC(v[inexact], digits = 17, format = "g"))
  out
}

#' Decode Exact Text Back to Atomic Values (Internal)
#'
#' Inverse of [encode_value()]. The token `"NA"` is a missing value for every
#' type except character, where the writer never emits a missing value and
#' `"NA"` is therefore the literal string.
#'
#' @param x Character vector.
#' @param type One of `"character"`, `"integer"`, `"double"`, `"logical"`.
#' @return A vector of that type.
#' @keywords internal
decode_value <- function(x, type) {
  if (type == "character") return(x)
  x[x == "NA"] <- NA_character_
  switch(type, integer = as.integer(x), double = as.numeric(x), logical = as.logical(x))
}

#' Read the Provenance Sidecar (Internal)
#'
#' Inverse of [write_provenance_sidecar()].
#'
#' @param file String. Path of the sidecar itself.
#' @return A named list, with `NULL` entries preserved.
#' @keywords internal
read_provenance_sidecar <- function(file) {
  table <- data.table::fread(file, sep = "\t", colClasses = "character", na.strings = NULL)
  data.table::setDF(table)
  if (!identical(names(table), c("key", "type", "index", "value"))) {
    stop(sprintf("'%s' is not a coreact provenance sidecar.", file))
  }
  types <- c("character", "integer", "double", "logical")

  out <- list()
  for (key in unique(table$key)) {
    part <- table[table$key == key, ]
    type <- unique(part$type)
    if (length(type) != 1 || !type %in% c(types, "null")) {
      stop(sprintf("Provenance entry '%s' in '%s' has an invalid type.", key, file))
    }
    out[key] <- list(
      if (type == "null") NULL
      else if (all(part$index == "0")) decode_value(character(0), type)
      else decode_value(part$value, type)
    )
  }
  out
}

#' Check if a Matrix is Binary (0/1)
#'
#' @description
#' Efficiently determines whether a matrix (dense or sparse) contains only
#' values of 0 and 1.
#'
#' @details
#' This function is optimized for performance depending on the input type:
#' \itemize{
#'   \item **Sparse Matrices (`dgCMatrix`):** It inspects only the stored values
#'   in the `@x` slot. Since implicit (non-stored) values are always 0, verifying
#'   that explicit values are either 0 or 1 is sufficient. This avoids expanding
#'   the matrix.
#'   \item **Dense Matrices:** It uses `any()` with short-circuiting logic. The
#'   scan stops immediately upon finding a value that is not 0 or 1.
#' }
#'
#' @param x A matrix-like object (base `matrix` or `Matrix::sparseMatrix`).
#'
#' @return `TRUE` if all elements are 0 or 1; otherwise `FALSE`. `NA` values
#'   count as non-binary.
#'
#' @keywords internal
is_binary_matrix <- function(x) {
  if (inherits(x, "nsparseMatrix")) {
    # Pattern matrices store no values: every entry is TRUE/FALSE by construction
    return(TRUE)
  }
  if (inherits(x, "sparseMatrix")) {
    # Check only stored values (@x slot)
    # If a value is stored, it must be 0 or 1.
    # (Implicit values are always 0, which is valid).
    return(all(x@x %in% c(0, 1)))
  }
  # Dense (base matrix or dense Matrix): NA is neither 0 nor 1
  vals <- as.vector(x)
  !anyNA(vals) && !any(vals != 0 & vals != 1)
}

#' Warn When Two Inputs Share Identical Features (Internal)
#'
#' @description
#' Detects the case where two inputs are really one library, which
#' [run_coreact()] cannot refuse because the two objects are not identical:
#' a copy of the same file under another name, or one library filtered twice
#' with different prevalence bounds.
#'
#' @details
#' Every feature present and identical on both sides is paired with itself,
#' and every two of them are paired twice, once in each orientation. Both
#' enter the family that [adjust_fdr()] corrects over, so the adjusted
#' p-values are computed over more hypotheses than exist.
#'
#' A feature counts only if its prevalence is at least 2. Under row-number
#' feature IDs two unrelated libraries can match on a single positive sample
#' by chance, and features absent from every sample are identical trivially.
#' Only the shared IDs are compared, so the cost is proportional to the
#' overlap and nil when the two inputs share no IDs.
#'
#' @param data_x,data_y `coreact_data` objects.
#' @return `NULL`, invisibly. Emits a warning when the count is positive.
#' @keywords internal
warn_shared_identical_features <- function(data_x, data_y) {
  shared <- intersect(rownames(data_x$mat), rownames(data_y$mat))
  if (length(shared) == 0) return(invisible(NULL))
  if (ncol(data_x$mat) != ncol(data_y$mat)) return(invisible(NULL))

  mat_x <- data_x$mat[shared, , drop = FALSE]
  mat_y <- data_y$mat[shared, , drop = FALSE]

  # Sparse difference: a row is identical when every entry cancels
  same <- Matrix::rowSums(abs(mat_x - mat_y)) == 0
  informative <- Matrix::rowSums(mat_x > 0) >= 2
  hits <- shared[same & informative]
  n_hits <- length(hits)
  if (n_hits == 0) return(invisible(NULL))

  examples <- paste(utils::head(hits, 3), collapse = ", ")
  if (n_hits > 3) examples <- paste0(examples, ", ...")

  # Subsetting a coreact_data object appends "_subset" to its name, so the two
  # inputs can carry the same name; name them once when they do.
  subject <- if (identical(data_x$name, data_y$name)) {
    sprintf("Both inputs ('%s')", data_x$name)
  } else {
    sprintf("'%s' and '%s'", data_x$name, data_y$name)
  }
  counted <- if (n_hits == 1) "1 feature that is" else sprintf("%d features that are", n_hits)

  warning(sprintf(
    paste0("%s share %s identical across all samples (%s). Each is paired with itself, ",
           "and each two of them are paired twice, inside a family of %s pairs. If these ",
           "inputs are one library, pass it as 'data_x' alone and every pair is reported once."),
    subject, counted, examples,
    format(as.numeric(nrow(data_x$mat)) * as.numeric(nrow(data_y$mat)), scientific = FALSE)
  ), call. = FALSE)
  invisible(NULL)
}

#' Format Column Mismatch Error (Internal Helper)
#'
#' @param mat_x Matrix X
#' @param mat_y Matrix Y
#' @param name_x Name of dataset X
#' @param name_y Name of dataset Y
#' @importFrom utils head
#' @return A formatted error string detailing specific mismatches.
#' @keywords internal
format_col_mismatch <- function(mat_x, mat_y, name_x, name_y) {
  cx <- colnames(mat_x)
  cy <- colnames(mat_y)

  # 1. Identify discrepancies
  in_x_not_y <- setdiff(cx, cy)
  in_y_not_x <- setdiff(cy, cx)

  # Base error message
  msg <- sprintf("Error: Sample identifiers (columns) in '%s' and '%s' are not identical.\n",
                 name_x, name_y)

  # 2. Case: Identical content, wrong order
  if (length(in_x_not_y) == 0 && length(in_y_not_x) == 0) {
    return(paste0(msg, "The datasets contain the same sample names, but in a different order. ",
                  "Please sort or align the columns to match exactly."))
  }

  # 3. Case: Content Mismatch - Helper to format list
  # Returns: "ColName (Idx 5), OtherCol (Idx 9)..."
  fmt_list <- function(diffs, all_cols) {
    idxs <- match(diffs, all_cols)
    items <- paste0("'", diffs, "' (Idx ", idxs, ")")
    if (length(items) > 5) {
      paste0(paste(head(items, 5), collapse = ", "),
             sprintf(", ... (+%d more)", length(items) - 5))
    } else {
      paste(items, collapse = ", ")
    }
  }

  # 4. Append details
  msg <- paste0(msg, "If using 'sample_cols', ensure both datasets resolve to the exact same sample set.\n",
                "Details:\n")

  if (length(in_x_not_y) > 0) {
    msg <- paste0(msg, sprintf("Found in '%s' but missing in '%s' (%d cols): %s\n",
                               name_x, name_y, length(in_x_not_y), fmt_list(in_x_not_y, cx)))
  }
  if (length(in_y_not_x) > 0) {
    msg <- paste0(msg, sprintf("Found in '%s' but missing in '%s' (%d cols): %s\n",
                               name_y, name_x, length(in_y_not_x), fmt_list(in_y_not_x, cy)))
  }

  msg
}

#' Preview Column Names and Indices
#'
#' A helper utility to inspect the headers of two TSV files side-by-side.
#' This assists in identifying the correct indices for `meta_cols`, `sample_cols`,
#' or `feature_ids` before running the pipeline.
#'
#' @param path_x Path to the first TSV file (dataset X).
#' @param path_y Path to the second TSV file (dataset Y).
#' @importFrom utils head
#' @importFrom data.table fread
#' @return A data.frame showing indices and column names for both files side-by-side.
#'   Returns invisibly; prints to console by default.
#' @export
#' @examples
#' preview_cols(system.file("extdata", "input_x.tsv", package = "coreact"),
#'              system.file("extdata", "input_y.tsv", package = "coreact"))
preview_cols <- function(path_x, path_y) {

  if (!file.exists(path_x)) stop("File not found: ", path_x)
  if (!file.exists(path_y)) stop("File not found: ", path_y)

  # Read only headers efficiently
  cols_x <- names(data.table::fread(path_x, nrows = 0))
  cols_y <- names(data.table::fread(path_y, nrows = 0))

  len_x <- length(cols_x)
  len_y <- length(cols_y)
  max_len <- max(len_x, len_y)

  # Pad shorter vector with NA
  pad_vec <- function(v, n) {
    if (length(v) < n) c(v, rep(NA, n - length(v))) else v
  }

  # Pad indices for display
  idx_x <- seq_len(len_x)
  idx_y <- seq_len(len_y)

  # Construct display table
  df <- data.frame(
    Idx_X  = pad_vec(idx_x, max_len),
    Name_X = pad_vec(cols_x, max_len),
    " | "  = rep("|", max_len), # Visual separator
    Idx_Y  = pad_vec(idx_y, max_len),
    Name_Y = pad_vec(cols_y, max_len),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # Handle NA indices to look cleaner (optional, depends on preference)
  df$Idx_X[is.na(df$Idx_X)] <- ""
  df$Idx_Y[is.na(df$Idx_Y)] <- ""
  df$Name_X[is.na(df$Name_X)] <- ""
  df$Name_Y[is.na(df$Name_Y)] <- ""

  print(df, row.names = FALSE, right = FALSE)
  invisible(df)
}
