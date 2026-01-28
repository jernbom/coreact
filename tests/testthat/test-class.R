test_that("new_coreact_data creates valid object", {
  mat <- Matrix::Matrix(matrix(c(1,0,1,0), nrow=2), sparse=TRUE)
  meta <- data.frame(id = c("A", "B"))

  obj <- new_coreact_data(mat, meta, name = "Test")

  expect_s3_class(obj, "coreact_data")
  expect_equal(obj$n_features, 2)
  expect_equal(obj$n_samples, 2)
  expect_equal(as.numeric(obj$prevalence), c(1, 1)) # Row sums
})

test_that("new_coreact_data catches errors", {
  mat <- matrix(1:4, nrow=2)
  meta <- data.frame(id = "A") # Mismatch rows

  expect_error(new_coreact_data(mat, meta), "Dimension mismatch")
})

test_that("filter_prevalence works with counts and percents", {
  # 5 samples
  # Feat 1: 5/5
  # Feat 2: 2/5
  # Feat 3: 1/5
  mat <- Matrix::Matrix(matrix(
    c(1,1,1,1,1,
      1,1,0,0,0,
      1,0,0,0,0), nrow=3, byrow=TRUE), sparse=TRUE)

  meta <- data.frame(id = 1:3)
  obj <- new_coreact_data(mat, meta)

  # Filter 1: Count >= 2 (Should keep feat 1 and 2)
  res_count <- filter_prevalence(obj, threshold = 2)
  expect_equal(res_count$n_features, 2)

  # Filter 2: Percent >= 0.5 (50%) -> Needs >= 2.5 samples -> 3 samples
  # Should keep only feat 1
  res_perc <- filter_prevalence(obj, threshold = 0.5)
  expect_equal(res_perc$n_features, 1)

  # Filter 3: Too high (Keep none) -> returns NULL
  expect_warning(res_null <- filter_prevalence(obj, threshold = 10), "removed all features")
  expect_null(res_null)
})
