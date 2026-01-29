#' Coreact Data Constructor
#'
#' Internal constructor to create a `coreact_data` object.
#' Strict validation: Row names must be present, unique, and identical in both inputs.
#'
#' @param mat A sparse matrix (dgCMatrix) containing counts/values.
#'   Rows are features, columns are samples.
#' @param meta A data.frame containing feature metadata (e.g., Gene Symbols).
#'   Number of rows must match `mat`.
#' @param name String. An identifier for the dataset.
#'
#' @return An S3 object of class `coreact_data`.
#' @importFrom methods is
#' @export
new_coreact_data <- function(mat, meta, name = "unknown") {

  # 1. Type Validation
  if (!methods::is(mat, "Matrix")) {
    stop("Input 'mat' must be a Matrix (preferably sparse).")
  }
  if (!is.data.frame(meta)) {
    stop("Input 'meta' must be a data.frame.")
  }

  # 2. Dimension Validation
  if (nrow(mat) != nrow(meta)) {
    stop(sprintf("Dimension mismatch: Matrix has %d rows but metadata has %d rows.",
                 nrow(mat), nrow(meta)))
  }

  # 3. Row Name Validation (Strict)
  mat_names <- rownames(mat)
  meta_names <- rownames(meta)

  if (is.null(mat_names)) {
    stop("Input 'mat' must have row names (Feature IDs).")
  }

  # Check Identity
  if (!identical(mat_names, meta_names)) {
    stop("Row names (Feature IDs) in Matrix and Metadata must be identical and in the same order.")
  }

  if (any(duplicated(rownames(mat))) || any(duplicated(rownames(meta)))) stop("Row names must be unique.")

  # 4. Create Object
  structure(
    list(
      mat = mat,
      meta = meta,
      name = name
    ),
    class = "coreact_data"
  )
}

#' Print Method for Coreact Data
#'
#' Displays a summary of the Coreact object.
#'
#' @param x A `coreact_data` object.
#' @param ... Additional arguments (ignored).
#' @export
print.coreact_data <- function(x, ...) {
  cat(sprintf("=== Coreact Data Object: [%s] ===\n", x$name))
  cat(sprintf("Dimensions: %d Features x %d Samples\n", nrow(x$mat), ncol(x$mat)))

  mem_size <- format(object.size(x), units = "auto")
  cat(sprintf("Size: %s\n", mem_size))

  if (!is.null(rownames(x$mat))) {
    n_show <- min(3, nrow(x$mat))
    ids <- rownames(x$mat)[1:n_show]
    cat(sprintf("Features: %s%s\n", paste(ids, collapse = ", "),
                if (nrow(x$mat) > 3) " ..." else ""))
  }
  invisible(x)
}

#' Subset Coreact Data
#'
#' Subsetting operator.
#'
#' @param x A `coreact_data` object.
#' @param i Row indices (features).
#' @param j Column indices (samples).
#' @param ... Additional arguments (ignored).
#'
#' @export
`[.coreact_data` <- function(x, i, j, ...) {
  if (missing(i)) i <- seq_len(nrow(x$mat))
  if (missing(j)) j <- seq_len(ncol(x$mat))

  new_mat <- x$mat[i, j, drop = FALSE]
  new_meta <- x$meta[i, , drop = FALSE]
  new_name <- paste0(x$name, "_subset")

  # Constructor will validate that the subsetted names still match
  new_coreact_data(new_mat, new_meta, name = new_name)
}

#' Combine Coreact Data Objects (Column-bind)
#'
#' Combines multiple coreact_data objects by adding samples.
#'
#' @param ... `coreact_data` objects to combine.
#' @param deparse.level Integer. Ignored.
#'
#' @return A combined `coreact_data` object.
#' @export
cbind.coreact_data <- function(..., deparse.level = 1) {
  objects <- list(...)

  if (length(objects) == 0) return(NULL)
  if (length(objects) == 1) return(objects[[1]])

  ref_ids <- rownames(objects[[1]]$mat)

  for (idx in 2:length(objects)) {
    curr_obj <- objects[[idx]]
    curr_ids <- rownames(curr_obj$mat)

    if (!identical(ref_ids, curr_ids)) {
      stop(sprintf(
        "Cannot cbind objects: Feature IDs (row names) do not match between object 1 and object %d.",
        idx
      ))
    }
  }

  all_mats <- lapply(objects, function(x) x$mat)
  combined_mat <- do.call(cbind, all_mats)
  combined_meta <- objects[[1]]$meta
  new_name <- paste(sapply(objects, function(x) x$name), collapse = "+")

  new_coreact_data(combined_mat, combined_meta, name = new_name)
}
