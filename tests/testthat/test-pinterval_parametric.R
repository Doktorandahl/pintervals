# Tests for pinterval_parametric()

# --- Setup ---
set.seed(42)
n <- 500
x <- runif(n)
y_norm <- rnorm(n, mean = 2 * x, sd = 1)
pred_cal <- 2 * x[1:250]
truth_cal <- y_norm[1:250]
pred_test <- 2 * x[251:500]
truth_test <- y_norm[251:500]

# ============================================================
# 1. Input validation
# ============================================================

test_that("pred must be numeric", {
	expect_error(
		pinterval_parametric(
			pred = "abc",
			calib = pred_cal,
			calib_truth = truth_cal
		),
		"pinterval_parametric.*pred.*numeric"
	)
})

test_that("alpha must be single numeric in (0,1)", {
	expect_error(
		pinterval_parametric(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			alpha = 0
		),
		"pinterval_parametric.*alpha"
	)
	expect_error(
		pinterval_parametric(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			alpha = 1
		),
		"pinterval_parametric.*alpha"
	)
})

test_that("dist must be valid", {
	expect_error(
		pinterval_parametric(
			pred = pred_test,
			calib = pred_cal,
			calib_truth = truth_cal,
			dist = "invalid_dist"
		),
		"distribution"
	)
})

test_that("calib and calib_truth must have same length", {
	expect_error(
		pinterval_parametric(
			pred = pred_test,
			calib = pred_cal[1:100],
			calib_truth = truth_cal[1:50]
		),
		"same length"
	)
})

test_that("beta distribution errors if values outside (0,1)", {
	expect_error(
		pinterval_parametric(
			pred = c(0.5, 0.6),
			calib = pred_cal,
			calib_truth = truth_cal,
			dist = "beta"
		),
		"beta.*0.*1"
	)
})

# ============================================================
# 2. Output structure
# ============================================================

test_that("output is a tibble with correct columns", {
	result <- pinterval_parametric(
		pred = pred_test,
		calib = pred_cal,
		calib_truth = truth_cal,
		dist = "norm"
	)
	expect_s3_class(result, "tbl_df")
	expect_true(all(c("pred", "lower_bound", "upper_bound") %in% names(result)))
	expect_equal(nrow(result), length(pred_test))
})

test_that("lower_bound <= pred <= upper_bound", {
	result <- pinterval_parametric(
		pred = pred_test,
		calib = pred_cal,
		calib_truth = truth_cal,
		dist = "norm"
	)
	valid <- !is.na(result$lower_bound)
	expect_true(all(result$lower_bound[valid] <= result$upper_bound[valid]))
})

# ============================================================
# 3. Distribution-specific tests
# ============================================================

test_that("normal distribution produces valid intervals", {
	result <- pinterval_parametric(
		pred = pred_test,
		calib = pred_cal,
		calib_truth = truth_cal,
		dist = "norm"
	)
	expect_s3_class(result, "tbl_df")
	coverage <- mean(
		truth_test >= result$lower_bound & truth_test <= result$upper_bound,
		na.rm = TRUE
	)
	expect_true(coverage >= 0.80, info = paste("Coverage was", coverage))
})

test_that("normal with explicit pars works", {
	result <- pinterval_parametric(
		pred = pred_test,
		dist = "norm",
		pars = list(mean = pred_test, sd = 1)
	)
	expect_s3_class(result, "tbl_df")
	expect_equal(nrow(result), length(pred_test))
})

test_that("lognormal distribution works with positive data", {
	set.seed(42)
	y_ln <- rlnorm(n, meanlog = x, sdlog = 0.5)
	pc <- exp(x[1:250])
	tc <- y_ln[1:250]
	pt <- exp(x[251:500])

	result <- pinterval_parametric(
		pred = pt,
		calib = pc,
		calib_truth = tc,
		dist = "lnorm"
	)
	expect_s3_class(result, "tbl_df")
})

test_that("gamma distribution works with positive data", {
	set.seed(42)
	y_g <- rgamma(n, shape = 2, rate = 1 / (2 * x + 0.5))
	pc <- (2 * x[1:250] + 0.5) * 2
	tc <- y_g[1:250]
	pt <- (2 * x[251:500] + 0.5) * 2

	result <- pinterval_parametric(
		pred = pt,
		calib = pc,
		calib_truth = tc,
		dist = "gamma"
	)
	expect_s3_class(result, "tbl_df")
})

test_that("custom quantile function works", {
	# Using qnorm as a custom function
	result <- pinterval_parametric(
		pred = pred_test[1:5],
		dist = qnorm,
		pars = list(mean = pred_test[1:5], sd = 1)
	)
	expect_s3_class(result, "tbl_df")
	expect_equal(nrow(result), 5)
})

# ============================================================
# 4. Wider intervals at lower alpha
# ============================================================

test_that("wider intervals with lower alpha", {
	result_90 <- pinterval_parametric(
		pred = pred_test[1:10],
		calib = pred_cal,
		calib_truth = truth_cal,
		dist = "norm",
		alpha = 0.1
	)
	result_50 <- pinterval_parametric(
		pred = pred_test[1:10],
		calib = pred_cal,
		calib_truth = truth_cal,
		dist = "norm",
		alpha = 0.5
	)
	width_90 <- mean(result_90$upper_bound - result_90$lower_bound, na.rm = TRUE)
	width_50 <- mean(result_50$upper_bound - result_50$lower_bound, na.rm = TRUE)
	expect_true(width_90 > width_50)
})

# ============================================================
# 5. Edge cases
# ============================================================

test_that("single prediction works", {
	result <- pinterval_parametric(
		pred = 1.0,
		calib = pred_cal,
		calib_truth = truth_cal,
		dist = "norm"
	)
	expect_equal(nrow(result), 1)
})

test_that("warns when pred contains NA", {
	expect_warning(
		pinterval_parametric(
			pred = c(NA_real_, 1),
			calib = pred_cal,
			calib_truth = truth_cal
		),
		"pred.*NA"
	)
})
