# coreact: High-Performance Co-reactivity Analysis

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
**coreact** is a high-performance R pipeline designed for the computationally efficient analysis of co-occurrence statistics between two large-scale reactome datasets (e.g., anti-viral vs. anti-human peptide/protein reactivity).

It addresses the dimensionality problem where both feature sets are large. By utilizing **adaptive matrix algebra**, **broadcasting**, and **sparse matrix representations**, `coreact` can calculate intersection, union, Jaccard indices, Overlap coefficients, and Hypergeometric P-values for millions of pairs on standard HPC nodes without memory overflow.

## ⚡ Key Features

* **Adaptive Engine:** Automatically switches between dense matrix operations (for fast block processing) and vectorized operations (for memory safety) based on data density.
* **Memory Efficiency:** Built on the `Matrix` package for sparse data handling. Uses **Fork-based parallelization** (Copy-On-Write) on Linux/macOS systems to share memory across workers without duplication.
* **Disc-to-Disc Pipeline:** A single wrapper function handles data loading, prevalence filtering, computation, FDR correction, and result export.
* **Statistical Rigor:** Computes one-tailed Hypergeometric P-values (enrichment) and supports standard FDR (Benjamini-Hochberg) correction.

## 📦 Installation

You can install the development version of coreact directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("jernbom/coreact")

```

*(Note: Requires a C++ compiler and the `Matrix` package, which are available on many HPC R environments.)*

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
  
  # Prevalence Filter (Important for speed/noise reduction)
  min_prevalence = c(0.1, 2), # 10% prevalence in X, N>=2 in Y
  
  # Structural Filters
  # Only keep pairs with at least 1 shared sample and at least 50% overlap
  filter_config = list(min_intersection = 1, min_overlap = 0.5),
  
  # Statistical Filter
  # For the purpose of this example, keep all hits regardless of significance.
  # NB! I strongly recommend lowering this threshold to, e.g., 0.05 for an actual use case!
  fdr_threshold = 1,
  
  # HPC Configuration. Use 1 if running on a Windows machine or running interactively in RStudio.
  n_cores = 1
)

```

## 📄 Input Data Format

`coreact` expects two TSV files (X and Y). The structure must be consistent:

1. **Metadata Columns:** Text or IDs (must be specified in `meta_cols`).
2. **Sample Columns:** Binary data (0 or 1).
3. **Consistency:** Sample column names must match exactly between X and Y.

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
| PEP_1001 | PROT_1|GENE_A | 1 | 1 | 1 | 1 | 5 | 5 | 0.200 | 1.000 | 0.500 | 0.829 |
| PEP_1002 | PROT_1|GENE_A | 2 | 1 | 2 | 3 | 5 | 6 | 0.333 | 0.667 | 0.500 | 0.829 |
| PEP_1003 | PROT_1|GENE_A | 3 | 1 | 1 | 2 | 5 | 6 | 0.167 | 0.500 | 0.778 | 0.829 |
| PEP_1004 | PROT_1|GENE_A | 4 | 1 | 2 | 2 | 5 | 5 | 0.400 | 1.000 | 0.222 | 0.829 |


Corresponding metadata is saved in sidecar files (linked by `feature_x` / `feature_y`):

**`results/interactions_metadata_x.tsv`**

| feature_id | u_pep_id | pep_nr | protein_id | gene_name |
| --- | --- | --- | --- | --- |
| PEP_1001 | PEP_1001 | 1 | PROT_1 | GENE_E |
| PEP_1002 | PEP_1002 | 2 | PROT_5 | GENE_M |
| PEP_1003 | PEP_1003 | 3 | PROT_1 | GENE_E |
| PEP_1004 | PEP_1004 | 4 | PROT_1 | GENE_T |


## 📐 Methodology

For every pair of features $(A, B)$:

1. **Intersection ($A \cap B$):** Number of samples where both features are present.
2. **Union ($A \cup B$):** Number of samples where $A$ or $B$ is present ($Size_A + Size_B - Intersection$
3. **Jaccard Index:**
  $$J = \frac{|A \cap B|}{|A \cup B|}$$
4. **Overlap Coefficient:**
  $$J = \frac{|A \cap B|}{min(|A|, |B|)}$$
5. **P-value (Hypergeometric):** Calculated using `phyper(q, m, n, k, lower.tail=FALSE)` where:
* $q = |A \cap B| - 1$
* $m = |A|$ (Prevalence of Feature A)
* $n = N_{samples} - |A|$
* $k = |B|$ (Prevalence of Feature B)

## ⚠️ Status & Attribution

**Status:** Active Development (Alpha).

**Author:** August F. Jernbom (Johns Hopkins Medicine)

This package is intended for research use. The API may change as the methods are refined for publication.

If you encounter issues or have feature requests, please open an issue in this repository.
