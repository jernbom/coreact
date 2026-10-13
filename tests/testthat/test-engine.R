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
  expect_equal(p1$union, 3)
  expect_equal(p1$jaccard, 2 / 3)
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
  expect_equal(p3$overlap, 0)

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

  # 4. Combined filters intersect
  res_both <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list(min_intersection = 2, min_jaccard = 0.5, min_overlap = 0.8)
  )
  expect_equal(nrow(res_both), 1)
  expect_equal(res_both$idx_x, 1)

  # 5. Stats computed via the filter path equal those computed afterwards
  res_none <- worker_chunk_calc(c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples, list())
  res_none <- res_none[order(res_none$idx_x, res_none$idx_y), ]
  res_filt <- worker_chunk_calc(c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
                                list(min_jaccard = 0, min_overlap = 0))
  res_filt <- res_filt[order(res_filt$idx_x, res_filt$idx_y), ]
  expect_equal(res_none, res_filt)
})

test_that("worker_chunk_calc returns NULL when nothing passes", {
  # Filter Intersection >= 10 (Impossible)
  res <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list(min_intersection = 10)
  )
  expect_null(res)

  # Same via the secondary filters (early-exit path not taken)
  res2 <- worker_chunk_calc(
    c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples,
    config = list(min_jaccard = 1.5)
  )
  expect_null(res2)
})

test_that("worker_chunk_calc computes p-values correctly", {
  # Use Pair 1 (X1-Y1)
  # q is the intersection minus one (1)
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

test_that("worker_chunk_calc reports absolute X indices for a chunk", {
  # Chunk covering absolute rows 11 and 12 of some larger X
  res <- worker_chunk_calc(
    c(11, 12), mat_x, mat_y, prev_x, prev_y, n_samples, config = list()
  )
  expect_setequal(res$idx_x, c(11, 12))
  expect_setequal(res$idx_y, c(1, 2))
  expect_type(res$idx_x, "double")
})

test_that("zero-prevalence features give NA ratios and never pass ratio filters", {
  mx <- Matrix::Matrix(c(0, 0, 0,
                         1, 1, 0), 2, 3, byrow = TRUE, sparse = TRUE)
  my <- Matrix::Matrix(c(0, 0, 0,
                         0, 1, 1), 2, 3, byrow = TRUE, sparse = TRUE)
  px <- Matrix::rowSums(mx)
  py <- Matrix::rowSums(my)

  # No filters: NA where undefined, never NaN
  res <- worker_chunk_calc(1:2, mx, my, px, py, 3, list())
  res <- res[order(res$idx_x, res$idx_y), ]
  expect_equal(nrow(res), 4)
  expect_true(is.na(res$jaccard[1]))         # 0/0 union
  expect_false(is.nan(res$jaccard[1]))
  expect_equal(res$jaccard[2:3], c(0, 0))    # 0 / positive union
  expect_true(all(is.na(res$overlap[1:3])))  # min size is 0
  expect_equal(res$overlap[4], 0.5)
  expect_equal(res$p_val[1:3], c(1, 1, 1))

  # Ratio filters without min_intersection used to fail with
  # the classic missing-value-where-TRUE/FALSE-needed error
  expect_no_error(r_j <- worker_chunk_calc(1:2, mx, my, px, py, 3, list(min_jaccard = 0.3)))
  expect_equal(nrow(r_j), 1)
  expect_equal(r_j$idx_x, 2)
  expect_equal(r_j$idx_y, 2)

  expect_no_error(r_o <- worker_chunk_calc(1:2, mx, my, px, py, 3, list(min_overlap = 0.5)))
  expect_equal(nrow(r_o), 1)

  # A threshold that nothing can satisfy returns NULL, not an error
  expect_null(worker_chunk_calc(1:2, mx, my, px, py, 3, list(min_jaccard = 0.9)))
  expect_null(worker_chunk_calc(1:2, mx, my, px, py, 3, list(min_overlap = 0.9)))
})

test_that("safe_ratio maps zero denominators to NA and keeps shape", {
  expect_equal(safe_ratio(c(1, 0, 2), c(2, 0, 0)), c(0.5, NA, NA))
  m <- safe_ratio(matrix(1:4, 2), matrix(c(1, 0, 2, 0), 2))
  expect_equal(dim(m), c(2, 2))
  expect_equal(m[, 1], c(1, NA))
  expect_equal(m[, 2], c(1.5, NA))
})

test_that("worker_chunk_calc matches a brute-force reference on random data", {
  withr::local_seed(1234)
  n <- 15
  mx <- Matrix::rsparsematrix(12, n, density = 0.35, rand.x = function(k) rep(1, k))
  my <- Matrix::rsparsematrix(9, n, density = 0.45, rand.x = function(k) rep(1, k))
  # Force one all-zero row in each to cover the NA paths
  mx[3, ] <- 0
  my[5, ] <- 0
  mx <- methods::as(mx, "CsparseMatrix")
  my <- methods::as(my, "CsparseMatrix")

  ref <- brute_force_coreact(mx, my)

  res <- worker_chunk_calc(
    seq_len(nrow(mx)), mx, my,
    as.numeric(Matrix::rowSums(mx)), as.numeric(Matrix::rowSums(my)),
    n, config = list()
  )
  res <- as.data.frame(res[order(res$idx_x, res$idx_y), ])
  rownames(res) <- NULL
  rownames(ref) <- NULL

  expect_equal(nrow(res), nrow(ref))
  for (col in c("intersection", "size_x", "size_y", "union", "jaccard", "overlap", "p_val")) {
    expect_equal(res[[col]], ref[[col]], info = col)
  }
})

test_that("worker_chunk_calc prunes on raw p-value without changing surviving stats", {
  full <- worker_chunk_calc(c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples, list())
  full <- full[order(full$idx_x, full$idx_y), ]
  thr <- 0.6
  pruned <- worker_chunk_calc(c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples, list(),
                              max_p_val = thr)
  pruned <- pruned[order(pruned$idx_x, pruned$idx_y), ]

  expect_true(all(pruned$p_val <= thr))
  expect_equal(nrow(pruned), sum(full$p_val <= thr))
  expect_equal(pruned, full[full$p_val <= thr, ])

  # Nothing below an impossible ceiling -> NULL
  expect_null(worker_chunk_calc(c(1, 2), mat_x, mat_y, prev_x, prev_y, n_samples, list(),
                                max_p_val = 0))
})

test_that("effect thresholds of 0 are no-ops, including for NA ratios", {
  mx <- Matrix::Matrix(c(0, 0, 0,
                         1, 1, 0), 2, 3, byrow = TRUE, sparse = TRUE)
  my <- Matrix::Matrix(c(0, 0, 0,
                         0, 1, 1), 2, 3, byrow = TRUE, sparse = TRUE)
  px <- Matrix::rowSums(mx)
  py <- Matrix::rowSums(my)
  none <- worker_chunk_calc(1:2, mx, my, px, py, 3, list())
  zero <- worker_chunk_calc(1:2, mx, my, px, py, 3,
                            list(min_intersection = 0, min_jaccard = 0, min_overlap = 0))
  expect_equal(nrow(zero), 4)
  expect_equal(zero[order(zero$idx_x, zero$idx_y), ], none[order(none$idx_x, none$idx_y), ])
})
