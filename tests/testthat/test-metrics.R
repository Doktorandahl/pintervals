# Tests for metric functions: interval_coverage, interval_miscoverage,
# interval_score, interval_width

# --- Setup ---
truth <- c(1, 2, 3, 4, 5)
lb <- c(0.5, 1.5, 2.5, 3.5, 4.5)
ub <- c(1.5, 2.5, 3.5, 4.5, 5.5)
# All 5 observations are covered by their intervals

lb_partial <- c(0.5, 1.5, 5.0, 3.5, 4.5)
ub_partial <- c(1.5, 2.5, 5.5, 4.5, 5.5)
# Observation 3 (truth=3) is NOT covered by [5.0, 5.5]

# ============================================================
# interval_coverage()
# ============================================================

test_that("interval_coverage: perfect coverage", {
	cov <- interval_coverage(truth = truth, lower_bound = lb, upper_bound = ub)
	expect_equal(cov, 1.0)
})

test_that("interval_coverage: partial coverage", {
	cov <- interval_coverage(truth = truth, lower_bound = lb_partial, upper_bound = ub_partial)
	expect_equal(cov, 0.8)
})

test_that("interval_coverage: return_vector works", {
	cov_vec <- interval_coverage(truth = truth, lower_bound = lb, upper_bound = ub, return_vector = TRUE)
	expect_length(cov_vec, 5)
	expect_true(all(cov_vec == TRUE))
})

test_that("interval_coverage: truth must be numeric", {
	expect_error(
		interval_coverage(truth = "a", lower_bound = lb, upper_bound = ub),
		"truth.*numeric"
	)
})

test_that("interval_coverage: lower_bound must be numeric", {
	expect_error(
		interval_coverage(truth = truth, lower_bound = "a", upper_bound = ub),
		"lower_bound.*numeric"
	)
})

test_that("interval_coverage: length mismatch", {
	expect_error(
		interval_coverage(truth = truth, lower_bound = lb[1:3], upper_bound = ub),
		"same length"
	)
})

test_that("interval_coverage: requires at least one of intervals or bounds", {
	expect_error(
		interval_coverage(truth = truth),
		"Either.*intervals.*lower_bound.*upper_bound"
	)
})

test_that("interval_coverage: works with intervals list-column", {
	intervals <- list(
		list(lower_bound = 0.5, upper_bound = 1.5),
		list(lower_bound = 1.5, upper_bound = 2.5),
		list(lower_bound = 2.5, upper_bound = 3.5),
		list(lower_bound = 3.5, upper_bound = 4.5),
		list(lower_bound = 4.5, upper_bound = 5.5)
	)
	cov <- interval_coverage(truth = truth, intervals = intervals)
	expect_equal(cov, 1.0)
})

test_that("interval_coverage: non-contiguous intervals", {
	# truth = 3, covered by second segment [2.5, 3.5]
	intervals <- list(
		NULL,
		NULL,
		list(lower_bound = c(0.5, 2.5), upper_bound = c(1.0, 3.5)),
		NULL,
		NULL
	)
	cov <- interval_coverage(
		truth = truth,
		lower_bound = lb,
		upper_bound = ub,
		intervals = intervals
	)
	# observation 3 should be covered by the non-contiguous interval
	expect_equal(cov, 1.0)
})

# ============================================================
# interval_miscoverage()
# ============================================================

test_that("interval_miscoverage: zero miscoverage with perfect coverage", {
	mc <- interval_miscoverage(truth = truth, lower_bound = lb, upper_bound = ub, alpha = 0.1)
	# coverage = 1.0, expected coverage = 0.9, miscoverage = 1.0 - 0.9 = 0.1
	expect_equal(mc, 0.1)
})

test_that("interval_miscoverage: correct calculation", {
	mc <- interval_miscoverage(truth = truth, lower_bound = lb_partial, upper_bound = ub_partial, alpha = 0.1)
	# coverage = 0.8, expected = 0.9, miscoverage = 0.8 - 0.9 = -0.1
	expect_equal(mc, -0.1)
})

test_that("interval_miscoverage: alpha must be in (0,1)", {
	expect_error(
		interval_miscoverage(truth = truth, lower_bound = lb, upper_bound = ub, alpha = 0),
		"alpha"
	)
	expect_error(
		interval_miscoverage(truth = truth, lower_bound = lb, upper_bound = ub, alpha = 1),
		"alpha"
	)
})

test_that("interval_miscoverage: length mismatch", {
	expect_error(
		interval_miscoverage(truth = truth, lower_bound = lb[1:3], upper_bound = ub, alpha = 0.1),
		"same length"
	)
})

test_that("interval_miscoverage: truth must be numeric", {
	expect_error(
		interval_miscoverage(truth = "a", lower_bound = lb, upper_bound = ub, alpha = 0.1),
		"truth.*numeric"
	)
})

# ============================================================
# interval_score()
# ============================================================

test_that("interval_score: perfect coverage yields width-only score", {
	is_val <- interval_score(truth = truth, lower_bound = lb, upper_bound = ub, alpha = 0.1)
	# When all covered, interval score = mean(ub - lb) = mean(1) = 1
	expect_equal(is_val, 1.0)
})

test_that("interval_score: penalty for undercoverage", {
	is_partial <- interval_score(truth = truth, lower_bound = lb_partial, upper_bound = ub_partial, alpha = 0.1)
	is_full <- interval_score(truth = truth, lower_bound = lb, upper_bound = ub, alpha = 0.1)
	expect_true(is_partial > is_full)
})

test_that("interval_score: return_vector works", {
	is_vec <- interval_score(truth = truth, lower_bound = lb, upper_bound = ub, alpha = 0.1, return_vector = TRUE)
	expect_length(is_vec, 5)
	expect_true(all(is_vec == 1.0))  # all widths are 1
})

test_that("interval_score: alpha must be valid", {
	expect_error(
		interval_score(truth = truth, lower_bound = lb, upper_bound = ub, alpha = 0),
		"alpha"
	)
})

test_that("interval_score: truth must be numeric", {
	expect_error(
		interval_score(truth = "a", lower_bound = lb, upper_bound = ub, alpha = 0.1),
		"truth.*numeric"
	)
})

test_that("interval_score: length mismatch", {
	expect_error(
		interval_score(truth = truth, lower_bound = lb[1:3], upper_bound = ub, alpha = 0.1),
		"same length"
	)
})

# ============================================================
# interval_width()
# ============================================================

test_that("interval_width: correct mean width", {
	w <- interval_width(lower_bound = lb, upper_bound = ub)
	expect_equal(w, 1.0)
})

test_that("interval_width: return_vector works", {
	w_vec <- interval_width(lower_bound = lb, upper_bound = ub, return_vector = TRUE)
	expect_length(w_vec, 5)
	expect_true(all(w_vec == 1.0))
})

test_that("interval_width: requires bounds or intervals", {
	expect_error(
		interval_width(),
		"Either.*intervals.*lower_bound.*upper_bound"
	)
})

test_that("interval_width: lower_bound must be numeric", {
	expect_error(
		interval_width(lower_bound = "a", upper_bound = ub),
		"lower_bound.*numeric"
	)
})

test_that("interval_width: length mismatch", {
	expect_error(
		interval_width(lower_bound = lb[1:3], upper_bound = ub),
		"same length"
	)
})

test_that("interval_width: works with intervals list-column", {
	intervals <- list(
		list(lower_bound = c(0, 2), upper_bound = c(1, 3)),  # total width = 1 + 1 = 2
		list(lower_bound = 1, upper_bound = 4)                 # width = 3
	)
	w <- interval_width(intervals = intervals)
	expect_equal(w, 2.5)  # mean(2, 3)
})
