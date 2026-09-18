# coreact (development version)

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
  `out_path`, `binarize`, `min_prevalence` and `filter_config`, before any
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
  `chunk_size`, `n_cores` and `filter_config` (a misspelt filter name was
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
  `utils` and `tools` are declared in `Imports`. Run `renv::snapshot()` to
  reconcile `renv.lock`.
* Test suite grown from 126 to 300+ expectations: seeded random data, a
  brute-force reference implementation as oracle, and invariance tests for
  parallel vs serial, X/Y swap and chunk size.
* GitHub Actions `R CMD check` workflow (Linux, macOS, Windows).

# coreact 1.2.0

* Strict binary validation, improved logging, and `R CMD check` fixes.

# coreact 1.0.0

* Core engine, pipeline, and documentation.
