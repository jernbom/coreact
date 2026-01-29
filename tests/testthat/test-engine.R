library(testthat)
library(Matrix)
library(tibble)

# --- Helper: Create Manual Test Data ---
# X (2 rows, 4 cols)
# R1: 1 1 1 0 (Sum 3)
# R2: 0 0 0 1 (Sum 1)
mat_x <- Matrix::Matrix(
  c(1, 1, 1, 0,
    0, 0, 0, 1),
  nrow = 2, byrow = TRUE, sparse = TRUE
)

# Y (2 rows, 4 cols)
# C1: 1 1 0 0 (Sum 2)
# C2: 0 0 1 1 (Sum 2)
mat_y <- Matrix::Matrix(
  c(1, 1, 0, 0,
    0, 0, 1, 1),
  nrow = 2, byrow = TRUE, sparse = TRUE
)

# Pre-calculate sums
prev_x <- Matrix::rowSums(mat_x)
prev_y <- Matrix::rowSums(mat_y)
n_samples <- 4

test_that("worker_chunk_calc computes correct stats (No Filters)", {

  # Run worker with NO filters
  res <- worker_chunk_calc(
    idx_x_chunk = c(1, 2), # Absolute indices
    mat_x_chunk = mat_x,
    mat_y = mat_y,
    prev_x_chunk = prev_x,
    prev_y = prev_y,
    n_samples = n_samples,
    config = list() # Empty config = no filters
  )

  # We expect 4 pairs (2x2)
  expect_equal(nrow(res), 4)

  # Sort for easy checking: X1-Y1, X1-Y2, X2-Y1, X2-Y2
  res <- res[order(res$idx_x, res$idx_y), ]

  # --- Pair 1: X1 vs Y1 ---
  # Int: 2 (cols 1,2) | Union: 3+2-2=3 | Jacc: 2/3 | Over: 2/2=1
  p1 <- res[1, ]
  expect_equal(p1$intersection, 2)
  expect_equal(p1$jaccard, 2/3)
  expect_equal(p1$overlap, 1)

  # --- Pair 2: X1 vs Y2 ---
  # Int: 1 (col 3) | Union: 3+2-1=4 | Jacc: 1/4 | Over: 1/2=0.5
  p2 <- res[2, ]
  expect_equal(p2$intersection, 1)
  expect_equal(p2$jaccard, 0.25)
  expect_equal(p2$overlap, 0.5)

  # --- Pair 3: X2 vs Y1 ---
  # Int: 0
  p3 <- res[3, ]
  expect_equal(p3$intersection, 0)
  expect_equal(p3$jaccard, 0)

  # --- Pair 4: X2 vs Y2 ---
  # Int: 1 (col 4) | Union: 1+2-1=2 | Jacc: 1/2=0.5 | Over: 1/1=1
  p4 <- res[4, ]
  expect_equal(p4$intersection, 1)
  expect_equal(p4$jaccard, 0.5)
  expect_equal(p4$overlap, 1)
})

test_that("worker_chunk_calc applies filters correctly", {

  # 1. Filter Intersection >= 2
  # Should keep only Pair 1 (X1-Y1) which has Int=2
  res_int <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list(min_intersection = 2)
  )
  expect_equal(nrow(res_int), 1)
  expect_equal(res_int$idx_x, 1)
  expect_equal(res_int$idx_y, 1)

  # 2. Filter Jaccard >= 0.5
  # Should keep:
  # Pair 1 (0.66)
  # Pair 4 (0.5)
  res_jacc <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list(min_jaccard = 0.5)
  )
  expect_equal(nrow(res_jacc), 2)
  expect_equal(res_jacc$idx_x, c(1, 2))
  expect_equal(res_jacc$idx_y, c(1, 2)) # Y indices

  # 3. Filter Overlap >= 0.8
  # Should keep:
  # Pair 1 (1.0)
  # Pair 4 (1.0)
  res_over <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list(min_overlap = 0.8)
  )
  expect_equal(nrow(res_over), 2)
})

test_that("worker_chunk_calc returns NULL when nothing passes", {
  # Filter Intersection >= 10 (Impossible)
  res <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list(min_intersection = 10)
  )
  expect_null(res)
})

test_that("worker_chunk_calc computes p-values correctly", {
  # Use Pair 1 (X1-Y1)
  # q = 2-1 = 1
  # m = 3 (Size X)
  # n = 4-3 = 1 (Samples - Size X)
  # k = 2 (Size Y)

  expected_pval <- stats::phyper(q = 1, m = 3, n = 1, k = 2, lower.tail = FALSE)

  res <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list()
  )

  # Get Pair 1
  actual_pval <- res$p_val[res$idx_x == 1 & res$idx_y == 1]

  expect_equal(actual_pval, expected_pval)
})
