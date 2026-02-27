# Tests for C++ distance functions (row_euclidean_distance, row_mahalanobis_distance)

# ============================================================
# row_euclidean_distance()
# ============================================================

test_that("row_euclidean_distance: basic correctness", {
	X <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, ncol = 2, byrow = TRUE)
	v <- c(0, 0)
	result <- row_euclidean_distance(X, v)
	expected <- c(1, 1, sqrt(2))
	expect_equal(result, expected, tolerance = 1e-10)
})

test_that("row_euclidean_distance: zero distance", {
	X <- matrix(c(1, 2, 1, 2), nrow = 2, ncol = 2, byrow = TRUE)
	v <- c(1, 2)
	result <- row_euclidean_distance(X, v)
	expect_equal(result, c(0, 0))
})

test_that("row_euclidean_distance: empty matrix returns empty vector", {
	X <- matrix(numeric(0), nrow = 0, ncol = 2)
	v <- c(1, 2)
	result <- row_euclidean_distance(X, v)
	expect_length(result, 0)
})

test_that("row_euclidean_distance: errors on dimension mismatch", {
	X <- matrix(1:6, nrow = 3, ncol = 2)
	v <- c(1, 2, 3) # 3 elements vs 2 columns
	expect_error(
		row_euclidean_distance(X, v),
		"length.*must match.*columns"
	)
})

test_that("row_euclidean_distance: NA in v produces error", {
	X <- matrix(1:4, nrow = 2, ncol = 2)
	v <- c(1, NA)
	expect_error(
		row_euclidean_distance(X, v),
		"NA"
	)
})

test_that("row_euclidean_distance: NA in X produces NA result", {
	X <- matrix(c(1, NA, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
	v <- c(0, 0)
	result <- row_euclidean_distance(X, v)
	expect_true(is.na(result[1]))
	expect_false(is.na(result[2]))
})

test_that("row_euclidean_distance: single row", {
	X <- matrix(c(3, 4), nrow = 1)
	v <- c(0, 0)
	result <- row_euclidean_distance(X, v)
	expect_equal(result, 5)
})

# ============================================================
# row_mahalanobis_distance()
# ============================================================

test_that("row_mahalanobis_distance: identity S_inv = euclidean", {
	X <- matrix(c(1, 0, 0, 1, 1, 1), nrow = 3, ncol = 2, byrow = TRUE)
	v <- c(0, 0)
	S_inv <- diag(2)
	result <- row_mahalanobis_distance(X, v, S_inv)
	expected <- c(1, 1, sqrt(2))
	expect_equal(result, expected, tolerance = 1e-10)
})

test_that("row_mahalanobis_distance: empty matrix returns empty", {
	X <- matrix(numeric(0), nrow = 0, ncol = 2)
	v <- c(1, 2)
	S_inv <- diag(2)
	result <- row_mahalanobis_distance(X, v, S_inv)
	expect_length(result, 0)
})

test_that("row_mahalanobis_distance: errors on v length mismatch", {
	X <- matrix(1:6, nrow = 3, ncol = 2)
	v <- c(1, 2, 3)
	S_inv <- diag(2)
	expect_error(
		row_mahalanobis_distance(X, v, S_inv),
		"length.*must match"
	)
})

test_that("row_mahalanobis_distance: errors on S_inv dimension mismatch", {
	X <- matrix(1:6, nrow = 3, ncol = 2)
	v <- c(1, 2)
	S_inv <- diag(3) # 3x3 instead of 2x2
	expect_error(
		row_mahalanobis_distance(X, v, S_inv),
		"square matrix"
	)
})

test_that("row_mahalanobis_distance: NA in v produces error", {
	X <- matrix(1:4, nrow = 2, ncol = 2)
	v <- c(NA, 2)
	S_inv <- diag(2)
	expect_error(
		row_mahalanobis_distance(X, v, S_inv),
		"NA"
	)
})

test_that("row_mahalanobis_distance: NA in X produces NA result", {
	X <- matrix(c(1, NA, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
	v <- c(0, 0)
	S_inv <- diag(2)
	result <- row_mahalanobis_distance(X, v, S_inv)
	expect_true(is.na(result[1]))
	expect_false(is.na(result[2]))
})

test_that("row_mahalanobis_distance: scaled S_inv gives scaled distances", {
	X <- matrix(c(1, 0), nrow = 1, ncol = 2)
	v <- c(0, 0)
	# With scale=4 on first dimension, mahalanobis distance should be 2*1 = 2
	S_inv <- matrix(c(4, 0, 0, 1), nrow = 2)
	result <- row_mahalanobis_distance(X, v, S_inv)
	expect_equal(result, 2.0, tolerance = 1e-10)
})

test_that("row_mahalanobis_distance: result is non-negative", {
	set.seed(42)
	X <- matrix(rnorm(200), nrow = 100, ncol = 2)
	v <- c(0, 0)
	S_inv <- solve(cov(X))
	result <- row_mahalanobis_distance(X, v, S_inv)
	expect_true(all(result >= 0, na.rm = TRUE))
})
