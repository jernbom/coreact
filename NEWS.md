# coreact 1.5.0

## New features

* `run_coreact_outcome()` associates each pair in a coreact result with a
  binary or continuous sample-level outcome. Each pair is classified by which
  cells of its 2x2 table are occupied, and tested on the contrast that
  occupancy supports: the interaction where all four cells are occupied, a
  conditional contrast for nested pairs, double-positive against
  double-negative for identical features. The choice never looks at the
  outcome. Tests are likelihood ratio (binary) or exact F (continuous) tests,
  adjusted once over all tested rows; every estimable contrast is reported.
  Separated logistic fits report infinite estimates with their p-values kept.
  The continuous path assumes the four cells are equally variable. Tested rows
  get 95% intervals that agree with their tests — profile likelihood for a
  binary outcome, exact t for a continuous one — and `ci = "all"` extends them
  to every estimable row.
* `coreact_outcome_pipeline()` runs the outcome association from files: it
  reads the results `coreact_pipeline()` wrote, rebuilds the inputs from their
  provenance sidecar, reads the outcome from a one-row-per-sample TSV (columns
  by name or index), and writes the table with a sidecar of the settings used.

* Results record how they were computed. `run_coreact()` attaches a
  `provenance` attribute holding the family size, the number of samples, the
  effect-size thresholds and `max_p_val`, whether the run was within-library,
  the coreact version, and for each input its name, feature count, number of
  positive cells and the prevalence cutoffs it was filtered with, as counts.
  `coreact_pipeline()` adds `fdr_threshold` and the import settings of each
  input, and writes the lot next to the results as `*_provenance.tsv`.
* `read_coreact_results()` reads a results file back with the column types it
  had in memory and with its `n_pairs` and `provenance` attributes, so
  `adjust_fdr()` can be applied to results from disk. It could not before,
  because `n_pairs` did not survive being written.

## Changes to existing behaviour

* `coreact_pipeline()` writes a third sidecar, `*_provenance.tsv`, beside the
  results and the two metadata files.

# coreact 1.4.0

## New features

* `run_coreact()` and `coreact_pipeline()` can compare a library with itself.
  Omit `data_y`, or pass a single entry in `paths`, and every pair of the
  library's features is tested once: the result satisfies `idx_x < idx_y`, no
  feature is paired with itself, and the family `adjust_fdr()` corrects over is
  `n * (n - 1) / 2` rather than `n^2`. Reporting both orientations would enter
  each hypothesis in the family twice, and a feature paired with itself
  co-occurs perfectly by construction, so those rows would top the ranking
  without being findings. Which feature of a pair is `x` and which is `y`
  follows index order and carries no other meaning.
* `run_coreact()` warns when its two inputs share features that are identical
  across all samples, with the count, examples and what it costs. This catches
  one library supplied as two objects — a copy of a file under another name, or
  the same library filtered twice — which cannot be refused outright because
  the two objects differ.

## Changes to existing behaviour

* `run_coreact(x, y)` and `coreact_pipeline()` now refuse inputs that hold the
  same data, directing the caller to the single-input form. For `paths` this is
  checked both on the resolved paths and on the imported matrices. Code that
  relied on `run_coreact(obj, obj)` returning all `n^2` ordered pairs must now
  call `run_coreact(obj)`, which returns each pair once.
* `coreact_pipeline()` accepts `paths` of length 1. With one input, `names`,
  `meta_cols`, `sample_cols`, `feature_ids`, `binarize`, `min_prevalence` and
  `max_prevalence` must each be a single value, and `names` defaults to
  `"Input"`.
* `import_coreact_tsv()`'s error for feature-ID columns outside `meta_cols` now
  also names `feature_ids`, the argument that sets them when the call came
  through `coreact_pipeline()`.

# coreact 1.3.0

## Procedure

* `coreact_pipeline()` is now explicitly structured as three steps: feature
  selection on marginal prevalence (`min_prevalence` and the new
  `max_prevalence`), significance testing of every pair of retained features
  with Benjamini-Hochberg adjustment over all `n_x * n_y` pairs and filtering
  on `fdr_threshold`, and effect-size filtering of the significant pairs with
  `effect_filters` (renamed from `filter_config`; entries `min_intersection`,
  `min_jaccard`, `min_overlap`). The effect-size filters now default to none.
* Pairs that cannot end up in the output are pruned inside the engine (raw
  p-value above `fdr_threshold`, effect size below threshold), so memory and
  disk usage scale with the number of reported pairs while the result is
  identical to computing all pairs.
* `run_coreact()` gains `max_p_val` and records the number of tested pairs in
  the `n_pairs` attribute of its result; the new `adjust_fdr()` adds `p_adj`
  over that family for in-memory use.

## Bug fixes

* `coreact_pipeline()` argument validation used `&&` where `||` was intended,
  so malformed `paths` / `names` were not caught. Now validated, together with
  `out_path`, `binarize`, `min_prevalence` and `effect_filters`, before any
  file is read.
* Structural filters (`min_jaccard`, `min_overlap`) failed with
  "missing value where TRUE/FALSE needed" when a feature was absent from every
  sample and `min_intersection` was not set. Ratios with a zero denominator
  are now `NA` (previously `NaN`) and never pass a threshold.
* An empty result was written without the `p_adj` column; the output schema
  is now identical whether or not any pair passes.
* `import_coreact_tsv()` always returns a `dgCMatrix`. `Matrix::Matrix()`
  returned `ddiMatrix` / `dsCMatrix` / `dtCMatrix` for square inputs that
  happened to be diagonal, symmetric or triangular.
* `run_coreact()` now checks that the sample columns of `data_x` and `data_y`
  are identical (previously only `coreact_pipeline()` did), returns the empty
  schema for zero-feature inputs instead of failing on `1:0`, and validates
  `chunk_size`, `n_cores` and `effect_filters` (a misspelt filter name was
  silently ignored, i.e. no filter was applied).
* `cbind()` on `coreact_data` objects refuses to combine objects with
  overlapping sample names.
* Zero-row `coreact_data` objects (e.g. `obj[integer(0), ]`) are now valid.
* When X had fewer features than Y (internal swap), the results file had
  `size_y` before `size_x`; the column order is now always the documented one.

## Improvements

* `coreact_pipeline()` returns the written results tibble invisibly.
* `new_coreact_data()` accepts a base `matrix` (converted to `dgCMatrix`).
* `dim()` method for `coreact_data`, so `nrow()` / `ncol()` work directly.
* Clearer import errors: non-numeric sample columns, `NA` values and
  non-binary values are each reported with the offending columns / counts and
  a hint to use `binarize = TRUE`.
* The prevalence filter logs the effective cutoff and the number of features
  kept.
* `run_coreact()` is now documented as the in-memory entry point, with every
  output column and the test described; runnable examples for all exported
  functions.
* Roxygen markdown enabled, so help pages render code and emphasis correctly.

## Internal

* `dplyr` is no longer a dependency (three trivial uses replaced by base R).
  `utils` and `tools` are declared in `Imports`, and `renv.lock` no longer
  pins `dplyr` or its exclusive dependencies.
* Test suite grown from 126 to 340+ expectations: seeded random data, a
  brute-force reference implementation as oracle, and invariance tests for
  parallel vs serial, X/Y swap and chunk size.
* GitHub Actions workflows for `R CMD check` (Linux, macOS, Windows) and lintr.

# coreact 1.2.0

* Strict binary validation, improved logging, and `R CMD check` fixes.

# coreact 1.0.0

* Core engine, pipeline, and documentation.
