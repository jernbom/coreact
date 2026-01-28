test_that("import_coreact_tsv reads and parses correctly", {
  # 1. Create a dummy TSV
  tf <- tempfile(fileext = ".tsv")
  df <- data.frame(
    gene_id = c("g1", "g2"),
    protein_name = c("p1", "p2"),
    S1 = c(1, 0),
    S2 = c(0, 1),
    S3 = c(1, 1)
  )
  write.table(df, tf, sep="\t", quote=FALSE, row.names=FALSE)

  # 2. Test Import with Name
  obj <- import_coreact_tsv(tf, meta_cols = c("gene_id", "protein_name"), name = "TestImport")

  # 3. Validation
  expect_s3_class(obj, "coreact_data")
  expect_equal(obj$n_features, 2)
  expect_equal(obj$n_samples, 3) # S1, S2, S3
  expect_equal(colnames(obj$meta), c("gene_id", "protein_name"))

  # Cleanup
  unlink(tf)
})

test_that("import_coreact_tsv handles index-based meta columns", {
  tf <- tempfile(fileext = ".tsv")
  df <- data.frame(ID=1:5, Val=c(1,1,0,0,0), note=5:1)
  write.table(df, tf, sep="\t", row.names=FALSE)

  obj <- import_coreact_tsv(tf, meta_cols = c(1,3))
  expect_equal(obj$n_samples, 1) # Only 'Val' is sample
  expect_equal(colnames(obj$meta), c("ID", "note"))

  unlink(tf)
})
