# Tests for internal helper functions

# ============================================================
# ncs_compute() and NCS functions
# ============================================================

test_that("ncs_compute: absolute_error works", {
	pred <- c(1, 2, 3)
	truth <- c(1.1, 2.2, 2.7)
	result <- pintervals:::ncs_compute("absolute_error", pred, truth)
	expect_equal(result, abs(pred - truth))
})

test_that("ncs_compute: raw_error works", {
	pred <- c(1, 2, 3)
	truth <- c(1.1, 2.2, 2.7)
	result <- pintervals:::ncs_compute("raw_error", pred, truth)
	expect_equal(result, truth - pred)
})

test_that("ncs_compute: relative_error works", {
	pred <- c(1, 2, 3)
	truth <- c(1.1, 2.2, 2.7)
	result <- pintervals:::ncs_compute("relative_error", pred, truth)
	expect_equal(result, abs((pred - truth) / pred))
})

test_that("ncs_compute: za_relative_error works", {
	pred <- c(1, 2, 3)
	truth <- c(1.1, 2.2, 2.7)
	result <- pintervals:::ncs_compute("za_relative_error", pred, truth)
	expect_equal(result, abs((pred - truth) / (1 + pred)))
})

test_that("ncs_compute: heterogeneous_error works with coefs", {
	pred <- c(1, 2, 3)
	truth <- c(1.1, 2.2, 2.7)
	coefs <- c(0.1, 0.5)  # intercept + slope
	result <- pintervals:::ncs_compute("heterogeneous_error", pred, truth, coefs = coefs)
	expected <- abs(pred - truth) / (coefs[1] + coefs[2] * pred)
	expect_equal(result, expected)
})

test_that("ncs_compute: unknown type errors", {
	expect_error(
		pintervals:::ncs_compute("nonexistent", c(1), c(1)),
		"unknown.*non-conformity"
	)
})

test_that("ncs_compute: heterogeneous_error requires coefs", {
	expect_error(
		pintervals:::ncs_compute("heterogeneous_error", c(1, 2), c(1, 2), coefs = NULL),
		"coefs.*provided"
	)
})

test_that("rel_error warns on zero prediction", {
	expect_warning(
		pintervals:::rel_error(c(0, 1), c(1, 2)),
		"zero.*Inf"
	)
})

test_that("heterogeneous_error warns on non-positive denominators", {
	expect_warning(
		pintervals:::heterogeneous_error(c(1, 2), c(1, 2), coefs = c(1, -1)),
		"non-positive"
	)
})

test_that("heterogeneous_error requires coefs of length 2", {
	expect_error(
		pintervals:::heterogeneous_error(c(1), c(1), coefs = c(1, 2, 3)),
		"length 2"
	)
})

# ============================================================
# resolve_weight_function()
# ============================================================

test_that("resolve_weight_function returns correct kernel functions", {
	gk <- pintervals:::resolve_weight_function("gaussian_kernel")
	expect_equal(gk(0), 1)
	expect_true(gk(1) < 1)

	ck <- pintervals:::resolve_weight_function("caucy_kernel")
	expect_equal(ck(0), 1)

	lk <- pintervals:::resolve_weight_function("logistic")
	expect_equal(lk(0), 0.5)

	rl <- pintervals:::resolve_weight_function("reciprocal_linear")
	expect_equal(rl(0), 1)
})

test_that("resolve_weight_function accepts a custom function", {
	custom_fn <- function(d) exp(-d)
	result <- pintervals:::resolve_weight_function(custom_fn)
	expect_true(is.function(result))
	expect_equal(result(0), 1)
})

test_that("resolve_weight_function errors on invalid string", {
	expect_error(
		pintervals:::resolve_weight_function("invalid_kernel"),
		"arg"
	)
})

# ============================================================
# validate_distance_inputs()
# ============================================================

test_that("validate_distance_inputs: errors on NULL inputs", {
	expect_error(
		pintervals:::validate_distance_inputs(NULL, matrix(1:4, 2, 2), 2, 2),
		"must be provided"
	)
	expect_error(
		pintervals:::validate_distance_inputs(matrix(1:4, 2, 2), NULL, 2, 2),
		"must be provided"
	)
})

test_that("validate_distance_inputs: errors on wrong type", {
	expect_error(
		pintervals:::validate_distance_inputs("abc", matrix(1:4, 2, 2), 2, 2),
		"must be a matrix"
	)
})

test_that("validate_distance_inputs: errors on wrong nrow", {
	expect_error(
		pintervals:::validate_distance_inputs(
			matrix(1:6, 3, 2), matrix(1:4, 2, 2), 2, 2
		),
		"rows.*match"
	)
})

test_that("validate_distance_inputs: errors on mismatched ncol", {
	expect_error(
		pintervals:::validate_distance_inputs(
			matrix(1:6, 2, 3), matrix(1:4, 2, 2), 2, 2
		),
		"same number of columns"
	)
})

test_that("validate_distance_inputs: accepts numeric vectors", {
	expect_silent(
		pintervals:::validate_distance_inputs(c(1, 2, 3), c(4, 5), 3, 2)
	)
})

# ============================================================
# bin_chopper()
# ============================================================

test_that("bin_chopper: produces correct number of bins", {
	set.seed(42)
	x <- rnorm(100)
	bins <- pintervals:::bin_chopper(x, nbins = 4)
	expect_equal(length(unique(bins)), 4)
	expect_equal(length(bins), 100)
})

test_that("bin_chopper: errors when nbins < 2", {
	expect_error(
		pintervals:::bin_chopper(1:10, nbins = 1),
		"nbins.*greater than 1"
	)
})

test_that("bin_chopper: errors when nbins > length(x)", {
	expect_error(
		pintervals:::bin_chopper(1:5, nbins = 10),
		"nbins.*less than or equal"
	)
})

test_that("bin_chopper: errors when x has single unique value", {
	expect_error(
		pintervals:::bin_chopper(rep(1, 10), nbins = 2),
		"more than one unique"
	)
})

test_that("bin_chopper: return_breaks works", {
	set.seed(42)
	x <- rnorm(100)
	brks <- pintervals:::bin_chopper(x, nbins = 4, return_breaks = TRUE)
	expect_true(is.numeric(brks))
	expect_equal(length(brks), 5)  # nbins + 1 breaks
	expect_equal(brks[1], -Inf)
	expect_equal(brks[length(brks)], Inf)
})
