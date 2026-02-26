# Tests for pinterval_bccp()

# --- Setup ---
set.seed(42)
n <- 800
x <- runif(n)
y <- rnorm(n, mean = 3 * x, sd = 1)
pred_cal <- 3 * x[1:400]
truth_cal <- y[1:400]
pred_test <- 3 * x[401:800]
truth_test <- y[401:800]

# Create bins based on truth quartiles
breaks <- quantile(truth_cal, probs = seq(0, 1, by = 0.25))
calib_bins <- cut(truth_cal, breaks = breaks, labels = FALSE, include.lowest = TRUE)

# ============================================================
# 1. Input validation
# ============================================================

test_that("pred must be numeric", {
	expect_error(
		pinterval_bccp(pred = "a", calib = pred_cal, calib_truth = truth_cal, calib_bins = calib_bins),
		"pinterval_bccp.*pred.*numeric"
	)
})

test_that("calib must be provided", {
	expect_error(
		pinterval_bccp(pred = pred_test, calib = NULL),
		"pinterval_bccp.*calib.*provided"
	)
})

test_that("alpha must be valid", {
	expect_error(
		pinterval_bccp(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			calib_bins = calib_bins, alpha = 0
		),
		"pinterval_bccp.*alpha"
	)
})

test_that("breaks must be sorted and numeric", {
	expect_error(
		pinterval_bccp(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			breaks = c(3, 1, 2)
		),
		"pinterval_bccp.*breaks.*sorted"
	)
	expect_error(
		pinterval_bccp(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			breaks = "abc"
		),
		"pinterval_bccp.*breaks.*numeric"
	)
})

test_that("right must be a single logical", {
	expect_error(
		pinterval_bccp(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			calib_bins = calib_bins, right = "yes"
		),
		"pinterval_bccp.*right.*logical"
	)
})

test_that("contiguize must be a single logical", {
	expect_error(
		pinterval_bccp(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			calib_bins = calib_bins, contiguize = "yes"
		),
		"pinterval_bccp.*contiguize.*logical"
	)
})

test_that("must have at least two bins", {
	one_bin <- rep(1, length(pred_cal))
	expect_error(
		pinterval_bccp(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			calib_bins = one_bin
		),
		"at least two bins"
	)
})

test_that("calib_bins must have same length as calib", {
	expect_error(
		pinterval_bccp(
			pred = pred_test, calib = pred_cal, calib_truth = truth_cal,
			calib_bins = calib_bins[1:10]
		),
		"same length"
	)
})

# ============================================================
# 2. Output structure
# ============================================================

test_that("contiguize=TRUE returns tibble with lower/upper bounds", {
	result <- suppressWarnings(pinterval_bccp(
		pred = pred_test[1:10], calib = pred_cal, calib_truth = truth_cal,
		calib_bins = calib_bins, alpha = 0.1, contiguize = TRUE
	))
	expect_s3_class(result, "tbl_df")
	expect_true(all(c("pred", "lower_bound", "upper_bound") %in% names(result)))
})

test_that("contiguize=FALSE returns tibble with intervals list-column", {
	result <- suppressWarnings(pinterval_bccp(
		pred = pred_test[1:10], calib = pred_cal, calib_truth = truth_cal,
		calib_bins = calib_bins, alpha = 0.1, contiguize = FALSE
	))
	expect_s3_class(result, "tbl_df")
	expect_true("intervals" %in% names(result))
})

# ============================================================
# 3. Breaks vs calib_bins
# ============================================================

test_that("breaks parameter works to define bins", {
	result <- pinterval_bccp(
		pred = pred_test[1:10], calib = pred_cal, calib_truth = truth_cal,
		breaks = as.numeric(breaks), alpha = 0.1, contiguize = TRUE
	)
	expect_s3_class(result, "tbl_df")
})

test_that("warns when both breaks and calib_bins provided", {
	expect_warning(
		pinterval_bccp(
			pred = pred_test[1:10], calib = pred_cal, calib_truth = truth_cal,
			calib_bins = calib_bins, breaks = as.numeric(breaks), alpha = 0.1, contiguize = TRUE
		),
		"breaks.*calib_bins.*ignored"
	)
})
