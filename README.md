# coreact: High-Performance Co-reactivity Analysis

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
**coreact** is a high-performance R pipeline designed for the computationally efficient analysis of co-occurrence statistics between two large-scale reactome datasets (e.g., anti-viral vs. anti-human peptide/protein reactivity).

It addresses the dimensionality problem where both feature sets are large. By utilizing **adaptive matrix algebra**, **broadcasting**, and **sparse matrix representations**, `coreact` can calculate intersection, union, Jaccard indices, Overlap coefficients, and Hypergeometric P-values for millions of pairs on standard HPC nodes without memory overflow.

## ⚡ Key Features

* **Chunked Matrix Engine:** Intersections for a block of X features against all of Y are computed in one sparse cross-product; the secondary matrices (Jaccard, overlap) are only materialised when they are needed for filtering, and the final statistics are computed vectorised on the surviving pairs.
* **Memory Efficiency:** Built on the `Matrix` package for sparse data handling. Uses **Fork-based parallelization** (Copy-On-Write) on Linux/macOS systems to share memory across workers without duplication. On Windows the computation runs serially.
* **Disc-to-Disc Pipeline:** A single wrapper function handles data loading, feature selection, computation, FDR correction, effect-size filtering and result export. The engine (`run_coreact()`) is also available for objects already in memory.
* **Statistical Rigor:** Computes one-tailed Hypergeometric P-values (enrichment) for every pair of retained features and controls the FDR over that full family with the Benjamini-Hochberg procedure; pairs that cannot reach significance are discarded inside the engine before they take up memory.

## 📦 Installation

You can install the development version of coreact directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("jernbom/coreact")

```

*(Note: `coreact` contains no compiled code of its own; it depends on `Matrix`, `data.table` and `tibble`, which are available as binaries on most HPC R environments.)*

## 🚀 Quick Start

The primary entry point is `coreact_pipeline`. The package includes small example datasets (`input_x.tsv` and `input_y.tsv`) so you can run the following code immediately to see how it works.

```r
library(coreact)

# 1. Locate the example files included in the package
path_x <- system.file("extdata", "input_x.tsv", package = "coreact")
path_y <- system.file("extdata", "input_y.tsv", package = "coreact")

# 2. Run the Pipeline
coreact_pipeline(
  # Input paths (Vector of length 2)
  paths = c(path_x, path_y),
  out_path = "results/interactions.tsv",
  
  # Metadata Configuration (List of length 2)
  # Specify which columns are metadata (indices or names).
  # The pipeline separates these from the numeric sample data.
  meta_cols = list(
    c("u_pep_id", "pep_nr", "protein_id", "gene_name"), # Metadata for X
    c("protein_id", "gene_name")                        # Metadata for Y
  ),
  
  # Feature ID Configuration (List of length 2)
  # Which column(s) uniquely identifies a row?
  feature_ids = list(
    "u_pep_id",    # Unique ID for X
    c("protein_id", "gene_name")   # Unique ID for Y
  ),
  
  # Step 1: Feature selection on marginal prevalence
  # (fraction of samples if < 1, absolute count if >= 1; max_prevalence is optional)
  min_prevalence = c(0.1, 2), # 10% prevalence in X, N>=2 in Y
  
  # Step 2: Significance filter (BH-adjusted p-value over all retained feature pairs)
  # For the purpose of this example, keep all hits regardless of significance.
  # NB! I strongly recommend lowering this threshold to, e.g., 0.05 for an actual use case!
  fdr_threshold = 1,
  
  # Step 3: Effect-size filters on the significant pairs (none by default)
  effect_filters = list(min_overlap = 0.5),
  
  # HPC Configuration. Use 1 if running on a Windows machine or running interactively in RStudio.
  n_cores = 1
)

```

`coreact_pipeline()` writes the results to disk and also returns them invisibly, so `res <- coreact_pipeline(...)` gives you the same table in-session.

### Working in memory

If the data are already loaded, or you want the unadjusted statistics as a tibble, use the importer and the engine directly:

```r
x <- import_coreact_tsv(path_x, meta_cols = 1:4, feature_id = "u_pep_id", name = "X")
y <- import_coreact_tsv(path_y, meta_cols = 1:2, feature_id = c("protein_id", "gene_name"), name = "Y")

res <- run_coreact(x, y)           # every pair; add max_p_val = 0.05 to prune on raw p
res <- adjust_fdr(res)             # BH over attr(res, "n_pairs") = nrow(x) * nrow(y) tests
sig <- res[res$p_adj <= 0.05, ]
```

`adjust_fdr()` uses the number of tested pairs as the size of the family, so pruning inside the engine (`max_p_val`, `effect_filters`) never changes the adjusted p-values of the significant pairs.

### Within a single library

To compare the features of one library with each other, supply it once: `run_coreact(x)`, or a single entry in `paths`. Each pair is then tested once and no feature is paired with itself, so the family is `n * (n - 1) / 2` rather than `n^2`. Supplying the same library twice is an error.

```r
res <- adjust_fdr(run_coreact(x))  # all(res$idx_x < res$idx_y) is TRUE
```

`preview_cols(path_x, path_y)` prints the column headers of both files side by side, which helps when choosing `meta_cols`, `sample_cols` and `feature_ids` by index.

## 📄 Input Data Format

`coreact` expects two TSV files (X and Y), or one for a within-library comparison. The structure must be consistent:

1. **Metadata Columns:** Text or IDs (must be specified in `meta_cols`).
2. **Sample Columns:** Binary data (0 or 1), numeric and without missing values. Count data can be thresholded at `>= 1` with `binarize = TRUE`.
3. **Consistency:** Sample column names must match exactly between X and Y (use `sample_cols` to select a shared subset).

**Example Input X (`input_x.tsv`):**

| u_pep_id | pep_nr | protein_id | gene_name | Sample_1 | Sample_2 | ... |
|----------|--------|------------|-----------|----------|----------|-----|
| PEP_1001 | 1      | PROT_1     | GENE_E    | 1        | 0        | ... |
| PEP_1002 | 2      | PROT_5     | GENE_M    | 0        | 1        | ... |

**Example Input Y (`input_y.tsv`):**

| protein_id | gene_name | Sample_1 | Sample_2 | ... |
|------------|-----------|----------|----------|-----|
| PROT_1     | GENE_A    | 1        | 1        | ... |
| PROT_2     | GENE_B    | 0        | 0        | ... |


## 📄 Output Data Format

The pipeline produces tab-separated files. Below is a preview of the main results file (`results/interactions.tsv`), showing significant co-occurrences between features from **Input X** (e.g., Peptides) and **Input Y** (e.g., Proteins).

**`results/interactions.tsv`**

| feature_x | feature_y | idx_x | idx_y | intersection | size_x | size_y | union | jaccard | overlap | p_val | p_adj |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PEP_1001 | PROT_1\|GENE_A | 1 | 1 | 1 | 1 | 5 | 5 | 0.200 | 1.000 | 0.500 | 1.000 |
| PEP_1002 | PROT_1\|GENE_A | 2 | 1 | 2 | 3 | 5 | 6 | 0.333 | 0.667 | 0.500 | 1.000 |
| PEP_1003 | PROT_1\|GENE_A | 3 | 1 | 1 | 2 | 5 | 6 | 0.167 | 0.500 | 0.778 | 1.000 |
| PEP_1004 | PROT_1\|GENE_A | 4 | 1 | 2 | 2 | 5 | 5 | 0.400 | 1.000 | 0.222 | 1.000 |

`idx_x` / `idx_y` are the row numbers of the features in the (prevalence-filtered) X and Y inputs. `jaccard` and `overlap` are `NA` when their denominator is zero, i.e. for features that are absent from every sample. (With ten samples and `fdr_threshold = 1`, nothing in the toy data reaches significance, hence `p_adj` of 1.)


Corresponding metadata is saved in sidecar files (linked by `feature_x` / `feature_y`):

**`results/interactions_metadata_x.tsv`**

| feature_id | u_pep_id | pep_nr | protein_id | gene_name |
| --- | --- | --- | --- | --- |
| PEP_1001 | PEP_1001 | 1 | PROT_1 | GENE_E |
| PEP_1002 | PEP_1002 | 2 | PROT_5 | GENE_M |
| PEP_1003 | PEP_1003 | 3 | PROT_1 | GENE_E |
| PEP_1004 | PEP_1004 | 4 | PROT_1 | GENE_T |

A provenance sidecar, **`results/interactions_provenance.tsv`**, records how the results were computed: the family size `n_pairs`, the prevalence cutoffs as counts, the FDR and effect-size thresholds, the import settings for each input, and the coreact version. `read_coreact_results()` reads the results back with it, so they can be adjusted again from disk:

```r
res <- read_coreact_results("results/interactions.tsv")
res <- adjust_fdr(res, method = "holm")   # same family, another method
```


## 📐 Methodology

The pipeline runs three steps in a fixed order:

1. **Feature selection.** Features are retained by their marginal prevalence (`min_prevalence`, `max_prevalence`). The hypergeometric test conditions on the prevalences, so selecting on them leaves the p-values of the retained pairs valid.
2. **Significance testing.** Every pair of retained features is tested, the p-values are adjusted over all $n_X \times n_Y$ pairs, and pairs with `p_adj <= fdr_threshold` are kept.
3. **Effect-size filtering.** The significant pairs are filtered by `effect_filters` (`min_intersection`, `min_jaccard`, `min_overlap`; none by default). With a moderate `fdr_threshold` these are mostly a convenience for keeping output sizes manageable.

Steps 2 and 3 are executed inside the engine so that pairs which cannot appear in the output are never materialised (a raw p-value above `fdr_threshold` cannot become significant after adjustment, and the effect-size thresholds commute with the FDR filter); the result is identical to computing all pairs first.

For every pair of features $(A, B)$:

1. **Intersection ($A \cap B$):** Number of samples where both features are present.
2. **Union ($A \cup B$):** Number of samples where $A$ or $B$ is present ($|A| + |B| - |A \cap B|$).
3. **Jaccard Index:**
  $$J = \frac{|A \cap B|}{|A \cup B|}$$
4. **Overlap Coefficient:**
  $$O = \frac{|A \cap B|}{\min(|A|, |B|)}$$
5. **P-value (Hypergeometric):** The one-sided probability of observing at least $|A \cap B|$ shared samples, calculated using `phyper(q, m, n, k, lower.tail = FALSE)` where:
* $q = |A \cap B| - 1$
* $m = |A|$ (Prevalence of Feature A)
* $n = N_{samples} - |A|$
* $k = |B|$ (Prevalence of Feature B)
6. **FDR (`p_adj`):** Benjamini-Hochberg adjustment (`p.adjust(method = "BH", n = n_X * n_Y)`) over all pairs of retained features.

## 🛠️ Development

```r
devtools::load_all()
devtools::test()          # testthat suite, including a brute-force reference implementation
devtools::document()      # regenerate NAMESPACE and man/ after editing roxygen comments
devtools::check()         # R CMD check
lintr::lint_package()     # style; CI fails on any lint
```

Continuous integration runs `R CMD check` on Linux (devel, release, oldrel), macOS and Windows for every push and pull request. See `NEWS.md` for the changelog.

## ⚠️ Status & Attribution

**Status:** Active Development (Alpha).

**Author:** August F. Jernbom (Johns Hopkins Medicine)

This package is intended for research use. The API may change as the methods are refined for publication.

If you encounter issues or have feature requests, please open an issue in this repository.
