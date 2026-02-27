# Tests for pinterval_bootstrap()

# --- Setup ---
set.seed(42)
n <- 500
x <- runif(n)
y <- rnorm(n, mean = 2 * x, sd = 1)
pred_cal <- 2 * x[1:250]
truth_cal <- y[1:250]
pred_test <- 2 * x[251:500]
truth_test <- y[251:500]

# ============================================================
# 1. Input validation
# ============================================================

test_that("pred must be numeric", {
	expect_error(
		pinterval_bootstrap(pred = "a", calib = pred_cal, calib_truth = truth_cal),
		"pinterval_bootstrap.*pred.*numeric"
	)
})

test_that("calib must be provided", {
	expect_error(
		pinterval_bootstrap(pred = pred_test, calib = NULL),
		"pinterval_bootstrap.*calib.*provided"
	)
})

test_that("alpha must be single numeric in (0,1)", {
	expect_error(
		pinterval_bootstrap(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			alpha = 0
		),
		"pinterval_bootstrap.*alpha"
	)
	expect_error(
		pinterval_bootstrap(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			alpha = 1.5
		),
		"pinterval_bootstrap.*alpha"
	)
})

test_that("n_bootstraps must be a single positive integer", {
	expect_error(
		pinterval_bootstrap(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			n_bootstraps = -1
		),
		"pinterval_bootstrap.*n_bootstraps"
	)
	expect_error(
		pinterval_bootstrap(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			n_bootstraps = 0
		),
		"pinterval_bootstrap.*n_bootstraps"
	)
	expect_error(
		pinterval_bootstrap(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			n_bootstraps = 1.5
		),
		"pinterval_bootstrap.*n_bootstraps"
	)
})

test_that("error_type must be valid", {
	expect_error(
		pinterval_bootstrap(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			error_type = "invalid"
		),
		"arg"
	)
})

# ============================================================
# 2. Output structure
# ============================================================

test_that("output is a tibble with correct columns and rows", {
	result <- pinterval_bootstrap(
		pred = pred_test,
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
	expect_true(all(c("pred", "lower_bound", "upper_bound") %in% names(result)))
	expect_equal(nrow(result), length(pred_test))
})

test_that("lower_bound <= upper_bound", {
	result <- pinterval_bootstrap(
		pred = pred_test,
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1
	)
	expect_true(all(result$lower_bound <= result$upper_bound, na.rm = TRUE))
})

# ============================================================
# 3. Correctness / coverage
# ============================================================

test_that("coverage is reasonable", {
	set.seed(123)
	result <- pinterval_bootstrap(
		pred = pred_test,
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1,
		n_bootstraps = 5000
	)
	coverage <- mean(
		truth_test >= result$lower_bound & truth_test <= result$upper_bound,
		na.rm = TRUE
	)
	# Bootstrap doesn't have formal coverage guarantee, but should be roughly close
	expect_true(coverage >= 0.80 && coverage <= 1.0)
})

test_that("wider intervals with lower alpha", {
	result_90 <- pinterval_bootstrap(
		pred = pred_test[1:10],
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1,
		n_bootstraps = 2000
	)
	result_50 <- pinterval_bootstrap(
		pred = pred_test[1:10],
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.5,
		n_bootstraps = 2000
	)
	width_90 <- mean(result_90$upper_bound - result_90$lower_bound, na.rm = TRUE)
	width_50 <- mean(result_50$upper_bound - result_50$lower_bound, na.rm = TRUE)
	expect_true(width_90 > width_50)
})

# ============================================================
# 4. Both error types produce valid output
# ============================================================

test_that("raw and absolute error types both produce valid output", {
	for (et in c("raw", "absolute")) {
		result <- pinterval_bootstrap(
			pred = pred_test[1:5],
			calib = pred_cal,
			calib_truth = truth_cal,
			error_type = et,
			alpha = 0.1
		)
		expect_s3_class(result, "tbl_df")
	}
})

# ============================================================
# 5. Alternative input formats
# ============================================================

test_that("calib as 2-column tibble works", {
	calib_tib <- tibble::tibble(pred = pred_cal, truth = truth_cal)
	result <- pinterval_bootstrap(
		pred = pred_test[1:5],
		calib = calib_tib,
		alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
})

test_that("calib as 2-column matrix works", {
	calib_mat <- cbind(pred_cal, truth_cal)
	result <- pinterval_bootstrap(
		pred = pred_test[1:5],
		calib = calib_mat,
		alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
})

# ============================================================
# 6. Edge cases
# ============================================================

test_that("single prediction works", {
	result <- pinterval_bootstrap(
		pred = 1.0,
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1
	)
	expect_equal(nrow(result), 1)
})

test_that("warns when pred contains NA", {
	expect_warning(
		pinterval_bootstrap(
			pred = c(NA_real_, 1),
			calib = pred_cal,
			calib_truth = truth_cal
		),
		"pred.*NA"
	)
})
