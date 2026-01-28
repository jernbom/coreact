# coreact: High-Performance Co-reactivity Analysis

[](https://www.google.com/search?q=https://lifecycle.r-lib.org/articles/stages.html%23experimental)
**coreact** is an R package designed for the computationally efficient analysis of co-occurrence statistics between two large-scale reactome datasets (e.g., anti-viral vs. anti-human peptide/protein reactivity).

It addresses the  dimensionality problem where both feature sets are large (e.g., viral peptides vs  human proteins). By utilizing **adaptive matrix algebra**, **broadcasting**, and **sparse matrix representations**, `coreact` can calculate intersection, union, Jaccard indices, Overlap coefficients, and Hypergeometric P-values for millions of pairs on standard HPC nodes without memory overflow.

## ⚡ Key Features

* **Adaptive Engine:** Automatically switches between dense matrix operations (for fast block processing) and vectorized operations (for memory safety) based on data density.
* **Memory Efficiency:** Built on the `Matrix` package for sparse data handling. Uses **Fork-based parallelization** (Copy-On-Write) on Linux/macOS systems to share memory across workers without duplication.
* **Disc-to-Disc Pipeline:** A single wrapper function handles data loading, prevalence filtering, computation, FDR correction, and result export.
* **Statistical Rigor:** Computes one-tailed Hypergeometric P-values (enrichment) and supports standard FDR (Benjamini-Hochberg) correction.

## 📦 Installation

This package is currently in **development** and is not yet on CRAN. You can install it directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("jernbom/coreact")

```

*(Note: Requires a C++ compiler and the `Matrix` package, which are available on many HPC R environments.)*

## 🚀 Quick Start: The Pipeline

The primary workflow is `coreact_pipeline`, which takes two TSV files and outputs a TSV result file.

```r
library(coreact)

# Define file paths
viral_file <- "data/raw/viral_counts.tsv"
human_file <- "data/raw/human_counts.tsv"
output_file <- "results/cooccurrence_hits.tsv"

# Run the pipeline
coreact_pipeline(
  path_x = viral_file,
  path_y = human_file,
  out_path = output_file,
  
  # Metadata columns (Index 1:10 are metadata in the standard Larman Lab PhIP-Seq file format)
  meta_cols_x = 1:10,
  meta_cols_y = 1:10,
  
  # Prevalence Filters (Important for speed/noise reduction)
  # Remove features present in fewer than 10% of samples
  min_prevalence = 0.1, 
  
  # Structural Filters
  # Only keep pairs with at least 10 shared samples and an Overlap coefficient of 50%.
  filter_config = list(min_intersection = 10, min_overlap = 0.5),
  
  # Statistical Filters
  # Keep only significant hits (FDR < 0.05)
  fdr_threshold = 0.05,
  
  # HPC Configuration
  n_cores = 16  # Adjust based on your SLURM allocation
)

```

## 📄 Input Data Format

The input files should be standard Tab-Separated Values (TSV).

* **Rows:** Features (Genes, Peptides, OTUs)
* **Columns:** Samples
* **Content:** Binary

**Example:**

| Feature_ID | Sample_01 | Sample_02 | Sample_03 | ... |
| :--- | :--- | :--- | :--- | :--- |
| **Gene_A** | 0 | 1 | 1 | ... |
| **Gene_B** | 0 | 0 | 0 | ... |
| **Gene_C** | 1 | 1 | 0 | ... |

*Note: The pipeline requires that the metadata columns are specified by `meta_cols` and treats all other columns as samples.*

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
