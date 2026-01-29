# data-raw/generate_readme_data.R

# Ensure the output directory exists
out_dir <- file.path(rprojroot::find_package_root_file(), "inst", "extdata")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

set.seed(42)

# --- Generate input_x.tsv ---
# 20 Features, 10 Samples
# Metadata: u_pep_id, pep_nr, protein_id, gene_name
n_feat_x <- 20
n_samp <- 10
samples <- paste0("Sample_", 1:n_samp)

df_x <- data.frame(
  u_pep_id   = paste0("PEP_", 1001:(1000+n_feat_x)),
  pep_nr     = 1:n_feat_x,
  protein_id = paste0("PROT_", sample(1:5, n_feat_x, replace=TRUE)),
  gene_name  = paste0("GENE_", sample(LETTERS, n_feat_x, replace=TRUE)),
  stringsAsFactors = FALSE
)

# FIXED: Use rbinom(size=1) to generate strict binary (0/1) data
# prob = 0.3 means roughly 30% of entries will be 1 (sparse)
mat_x <- matrix(stats::rbinom(n_feat_x * n_samp, size = 1, prob = 0.3),
                nrow = n_feat_x)
colnames(mat_x) <- samples
df_x <- cbind(df_x, mat_x)

write.table(df_x, file.path(out_dir, "input_x.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Created inst/extdata/input_x.tsv")

# --- Generate input_y.tsv ---
# 5 Features, 10 Samples
# Metadata: protein_id, gene_name
n_feat_y <- 5

df_y <- data.frame(
  protein_id = paste0("PROT_", 1:n_feat_y),
  gene_name  = paste0("GENE_", LETTERS[1:n_feat_y]),
  stringsAsFactors = FALSE
)

# FIXED: Use rbinom(size=1)
# prob = 0.5 means roughly 50% density
mat_y <- matrix(stats::rbinom(n_feat_y * n_samp, size = 1, prob = 0.5),
                nrow = n_feat_y)
colnames(mat_y) <- samples
df_y <- cbind(df_y, mat_y)

write.table(df_y, file.path(out_dir, "input_y.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Created inst/extdata/input_y.tsv")
