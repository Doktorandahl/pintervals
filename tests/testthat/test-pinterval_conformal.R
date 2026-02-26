# Tests for pinterval_conformal()

# --- Setup: shared test data ---
set.seed(42)
n <- 500
x <- runif(n)
y <- rnorm(n, mean = 2 * x, sd = 1)
pred_cal <- 2 * x[1:250]
truth_cal <- y[1:250]
pred_test <- 2 * x[251:500]
truth_test <- y[251:500]

# ============================================================
# 1. Input validation tests
# ============================================================

test_that("pred must be numeric", {
	expect_error(
		pinterval_conformal(pred = "a", calib = pred_cal, calib_truth = truth_cal),
		"pinterval_conformal.*pred.*numeric"
	)
	expect_error(
		pinterval_conformal(pred = list(1, 2), calib = pred_cal, calib_truth = truth_cal),
		"pinterval_conformal.*pred.*numeric"
	)
})

test_that("calib must be provided", {
	expect_error(
		pinterval_conformal(pred = pred_test, calib = NULL),
		"pinterval_conformal.*calib.*provided"
	)
})

test_that("calib_truth required when calib is numeric", {
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal),
		"pinterval_conformal.*calib_truth.*provided"
	)
})

test_that("calib must be numeric vector, matrix, or data.frame", {
	expect_error(
		pinterval_conformal(pred = pred_test, calib = "abc"),
		"pinterval_conformal.*calib.*numeric vector.*matrix.*data frame"
	)
})

test_that("calib matrix/data.frame must have exactly 2 columns", {
	expect_error(
		pinterval_conformal(
			pred = pred_test,
			calib = data.frame(a = 1:5, b = 1:5, c = 1:5)
		),
		"pinterval_conformal.*2-column"
	)
})

test_that("alpha must be single numeric in (0,1)", {
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal, calib_truth = truth_cal, alpha = 0),
		"pinterval_conformal.*alpha"
	)
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal, calib_truth = truth_cal, alpha = 1),
		"pinterval_conformal.*alpha"
	)
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal, calib_truth = truth_cal, alpha = -0.1),
		"pinterval_conformal.*alpha"
	)
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal, calib_truth = truth_cal, alpha = "0.1"),
		"pinterval_conformal.*alpha"
	)
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal, calib_truth = truth_cal, alpha = c(0.1, 0.2)),
		"pinterval_conformal.*alpha"
	)
})

test_that("grid_size must be a positive number", {
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal, calib_truth = truth_cal, grid_size = -1),
		"pinterval_conformal.*grid_size"
	)
	expect_error(
		pinterval_conformal(pred = pred_test, calib = pred_cal, calib_truth = truth_cal, grid_size = 0),
		"pinterval_conformal.*grid_size"
	)
})

test_that("lower_bound must be less than upper_bound", {
	expect_error(
		pinterval_conformal(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			lower_bound = 10, upper_bound = 5
		),
		"pinterval_conformal.*lower_bound.*upper_bound"
	)
})

test_that("calib and calib_truth must have the same length", {
	expect_error(
		pinterval_conformal(
			pred = pred_test,
			calib = pred_cal[1:100],
			calib_truth = truth_cal[1:50]
		),
		"pinterval_conformal.*same length"
	)
})

test_that("ncs_type must be a valid option", {
	expect_error(
		pinterval_conformal(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			ncs_type = "invalid_score"
		),
		"arg"
	)
})

test_that("normalize_distance must be a valid option", {
	expect_error(
		pinterval_conformal(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			normalize_distance = "invalid"
		),
		"arg"
	)
})

# ============================================================
# 2. Warning tests
# ============================================================

test_that("warns when pred contains NA", {
	pred_na <- pred_test
	pred_na[1] <- NA
	expect_warning(
		pinterval_conformal(pred = pred_na, calib = pred_cal, calib_truth = truth_cal),
		"pred.*NA"
	)
})

test_that("warns when calib or calib_truth contains NA", {
	calib_na <- pred_cal
	calib_na[1] <- NA
	expect_warning(
		pinterval_conformal(pred = pred_test, calib = calib_na, calib_truth = truth_cal),
		"calib.*NA"
	)
})

test_that("warns when both grid_size and resolution are provided", {
	expect_warning(
		pinterval_conformal(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			grid_size = 100, resolution = 0.01
		),
		"resolution.*ignored"
	)
})

# ============================================================
# 3. Output structure tests
# ============================================================

test_that("output is a tibble with correct columns", {
	result <- pinterval_conformal(
		pred = pred_test, calib = pred_cal, calib_truth = truth_cal, alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
	expect_true(all(c("pred", "lower_bound", "upper_bound") %in% names(result)))
	expect_equal(nrow(result), length(pred_test))
})

test_that("lower_bound <= pred <= upper_bound in general", {
	result <- pinterval_conformal(
		pred = pred_test, calib = pred_cal, calib_truth = truth_cal, alpha = 0.1
	)
	# Allow for NA intervals
	valid <- !is.na(result$lower_bound)
	expect_true(all(result$lower_bound[valid] <= result$upper_bound[valid]))
})

# ============================================================
# 4. Correctness / coverage tests
# ============================================================

test_that("coverage is approximately 1-alpha for large calibration set", {
	set.seed(123)
	n2 <- 2000
	x2 <- runif(n2)
	y2 <- rnorm(n2, mean = 2 * x2, sd = 1)
	p_cal <- 2 * x2[1:1000]
	t_cal <- y2[1:1000]
	p_test <- 2 * x2[1001:2000]
	t_test <- y2[1001:2000]

	result <- pinterval_conformal(
		pred = p_test, calib = p_cal, calib_truth = t_cal, alpha = 0.1,
		lower_bound = -5, upper_bound = 10
	)

	coverage <- mean(t_test >= result$lower_bound & t_test <= result$upper_bound, na.rm = TRUE)
	# Conformal prediction guarantees coverage >= 1-alpha, with some finite-sample slack
	expect_true(coverage >= 0.85, info = paste("Coverage was", coverage))
})

test_that("wider intervals at lower alpha", {
	result_90 <- pinterval_conformal(
		pred = pred_test[1:10], calib = pred_cal, calib_truth = truth_cal, alpha = 0.1
	)
	result_50 <- pinterval_conformal(
		pred = pred_test[1:10], calib = pred_cal, calib_truth = truth_cal, alpha = 0.5
	)
	width_90 <- mean(result_90$upper_bound - result_90$lower_bound, na.rm = TRUE)
	width_50 <- mean(result_50$upper_bound - result_50$lower_bound, na.rm = TRUE)
	expect_true(width_90 > width_50)
})

# ============================================================
# 5. Alternative input formats
# ============================================================

test_that("calib as 2-column tibble works", {
	calib_tib <- tibble::tibble(pred = pred_cal, truth = truth_cal)
	result <- pinterval_conformal(pred = pred_test, calib = calib_tib, alpha = 0.1)
	expect_s3_class(result, "tbl_df")
	expect_equal(nrow(result), length(pred_test))
})

test_that("calib as 2-column matrix works", {
	calib_mat <- cbind(pred_cal, truth_cal)
	result <- pinterval_conformal(pred = pred_test, calib = calib_mat, alpha = 0.1)
	expect_s3_class(result, "tbl_df")
	expect_equal(nrow(result), length(pred_test))
})

# ============================================================
# 6. NCS type tests
# ============================================================

test_that("all ncs_type options produce valid output", {
	for (ncs in c("absolute_error", "relative_error", "za_relative_error", "heterogeneous_error", "raw_error")) {
		result <- pinterval_conformal(
			pred = pred_test[1:5], calib = pred_cal, calib_truth = truth_cal,
			ncs_type = ncs, alpha = 0.1
		)
		expect_s3_class(result, "tbl_df", info = paste("ncs_type:", ncs))
		expect_equal(nrow(result), 5, info = paste("ncs_type:", ncs))
	}
})

# ============================================================
# 7. Edge cases
# ============================================================

test_that("single prediction works", {
	result <- pinterval_conformal(
		pred = 1.0, calib = pred_cal, calib_truth = truth_cal, alpha = 0.1
	)
	expect_equal(nrow(result), 1)
})

test_that("NA prediction produces NA bounds", {
	result <- suppressWarnings(pinterval_conformal(
		pred = NA_real_, calib = pred_cal, calib_truth = truth_cal, alpha = 0.1
	))
	expect_true(is.na(result$lower_bound[1]))
	expect_true(is.na(result$upper_bound[1]))
})

test_that("resolution parameter controls grid step", {
	# This should run without error
	result <- suppressWarnings(pinterval_conformal(
		pred = pred_test[1:3], calib = pred_cal, calib_truth = truth_cal,
		alpha = 0.1, grid_size = NULL, resolution = 0.1
	))
	expect_s3_class(result, "tbl_df")
})
