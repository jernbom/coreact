test_that("run_coreact handles swapping for efficiency", {
  # Create Large X and Small Y
  # Logic dictates X should be the larger one to split chunks effectively
  # But if user inputs Small X and Large Y, function should swap internally
  # and output correct IDs.

  # X (Small): 1 feature
  mat_x <- Matrix::Matrix(matrix(c(1,1), nrow=1), sparse=TRUE)
  meta_x <- data.frame(id="X1")
  obj_x <- new_coreact_data(mat_x, meta_x) # ID idx_x=1

  # Y (Large): 2 features
  mat_y <- Matrix::Matrix(matrix(c(1,1, 0,0), nrow=2, byrow=TRUE), sparse=TRUE)
  meta_y <- data.frame(id=c("Y1", "Y2"))
  obj_y <- new_coreact_data(mat_y, meta_y) # IDs idx_y=1,2

  # Run
  res <- run_coreact(obj_x, obj_y, n_cores = 1, chunk_size=10, filter_config=list(min_intersection=1))

  # Check results
  # X1 (1,1) matches Y1 (1,1). Should be idx_x=1, idx_y=1.
  expect_equal(nrow(res), 1)
  expect_equal(res$idx_x, 1)
  expect_equal(res$idx_y, 1)
})

test_that("coreact_pipeline runs end-to-end", {
  # Setup paths
  in_x <- tempfile()
  in_y <- tempfile()
  out_res <- tempfile()

  # Create Data
  # X: Feature A (1,1,0), B (0,0,1)
  write.table(data.frame(ID=c("A","B"), S1=c(1,0), S2=c(1,0), S3=c(0,1)),
              in_x, sep="\t", quote=FALSE, row.names=FALSE)

  # Y: Feature C (1,1,0) -> Perfect match with A
  write.table(data.frame(ID=c("C"), S1=c(1), S2=c(1), S3=c(0)),
              in_y, sep="\t", quote=FALSE, row.names=FALSE)

  # Run Pipeline
  # We expect A-C to match.
  coreact_pipeline(
    path_x = in_x,
    path_y = in_y,
    out_path = out_res,
    min_prevalence = 0, # Keep everything
    n_cores = 1
  )

  # Verify Output Exists
  expect_true(file.exists(out_res))

  # Verify Content
  res_df <- read.table(out_res, header=TRUE)
  expect_equal(nrow(res_df), 1)
  expect_equal(res_df$intersection, 2) # S1 and S2 match

  # Cleanup
  unlink(c(in_x, in_y, out_res))
})
