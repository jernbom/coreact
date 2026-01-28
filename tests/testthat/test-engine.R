test_that("Engine calculates stats correctly (Math check)", {
  # Scenario: 5 samples total
  # X: Present in [1, 2, 3] (Size = 3)
  # Y: Present in [3, 4, 5] (Size = 3)
  # Intersection: Sample 3 (Count = 1)
  # Union: 1,2,3,4,5 (Count = 5)

  mat_x <- Matrix::Matrix(matrix(c(1,1,1,0,0), nrow=1), sparse=TRUE)
  mat_y <- Matrix::Matrix(matrix(c(0,0,1,1,1), nrow=1), sparse=TRUE)

  prev_x <- 3
  prev_y <- 3
  n_samples <- 5

  # Config: Relaxed filters to ensure we get output
  config <- list(min_intersection = 0)

  # Note: accessing internal function with :::
  res <- coreact:::worker_chunk_calc(
    idx_x_chunk = 1,
    mat_x_chunk = mat_x,
    mat_y = mat_y,
    prev_x_chunk = prev_x,
    prev_y = prev_y,
    n_samples = n_samples,
    config = config
  )

  # Structural checks
  expect_equal(nrow(res), 1)
  expect_equal(res$intersection, 1)
  expect_equal(res$size_x, 3)
  expect_equal(res$size_y, 3)

  # Derived stats
  # Jaccard = Int / Union = 1 / 5 = 0.2
  expect_equal(res$jaccard, 0.2)

  # Overlap = Int / min(size_x, size_y) = 1 / 3
  expect_equal(res$overlap, 1/3)

  # P-value Check
  # phyper(q, m, n, k) -> q=int-1, m=size_x, n=total-size_x, k=size_y
  expected_p <- phyper(0, 3, 2, 3, lower.tail = FALSE)
  expect_equal(res$p_val, expected_p)
})

test_that("Engine respects filters", {
  # Case where intersection is 1, but we require 2
  mat_x <- Matrix::Matrix(matrix(c(1,0), nrow=1), sparse=TRUE)
  mat_y <- Matrix::Matrix(matrix(c(1,0), nrow=1), sparse=TRUE)

  config <- list(min_intersection = 2)

  res <- coreact:::worker_chunk_calc(1, mat_x, mat_y, 1, 1, 2, config)
  expect_null(res)
})
