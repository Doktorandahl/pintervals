# Tests for pinterval_mondrian()

# --- Setup ---
set.seed(42)
n <- 600
x <- runif(n)
group <- sample(c("A", "B", "C"), n, replace = TRUE)
mu <- ifelse(group == "A", 1 + x, ifelse(group == "B", 2 + x, 3 + x))
y <- rnorm(n, mean = mu, sd = 0.5)

pred_cal <- mu[1:300]
truth_cal <- y[1:300]
class_cal <- group[1:300]
pred_test <- mu[301:600]
truth_test <- y[301:600]
class_test <- group[301:600]

# ============================================================
# 1. Input validation
# ============================================================

test_that("pred must be numeric", {
	expect_error(
		pinterval_mondrian(pred = "a", pred_class = class_test, calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal),
		"pinterval_mondrian.*pred.*numeric"
	)
})

test_that("calib must be provided", {
	expect_error(
		pinterval_mondrian(pred = pred_test, pred_class = class_test, calib = NULL),
		"pinterval_mondrian.*calib.*provided"
	)
})

test_that("pred_class and pred must have the same length", {
	expect_error(
		pinterval_mondrian(
			pred = pred_test, pred_class = class_test[1:10],
			calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal
		),
		"pred_class.*same length"
	)
})

test_that("calib_class and calib must have the same length", {
	expect_error(
		pinterval_mondrian(
			pred = pred_test, pred_class = class_test,
			calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal[1:10]
		),
		"calib_class.*same length"
	)
})

test_that("alpha must be valid", {
	expect_error(
		pinterval_mondrian(
			pred = pred_test, pred_class = class_test,
			calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
			alpha = 0
		),
		"alpha"
	)
})

# ============================================================
# 2. Output structure
# ============================================================

test_that("output is a tibble with correct columns", {
	result <- pinterval_mondrian(
		pred = pred_test, pred_class = class_test,
		calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
		alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
	expect_true(all(c("pred", "lower_bound", "upper_bound") %in% names(result)))
	expect_equal(nrow(result), length(pred_test))
})

# ============================================================
# 3. Correctness / coverage
# ============================================================

test_that("per-class coverage is approximately 1-alpha", {
	set.seed(123)
	result <- pinterval_mondrian(
		pred = pred_test, pred_class = class_test,
		calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
		alpha = 0.1, lower_bound = -2, upper_bound = 8
	)

	for (cls in unique(class_test)) {
		idx <- class_test == cls
		coverage_cls <- mean(
			truth_test[idx] >= result$lower_bound[idx] &
				truth_test[idx] <= result$upper_bound[idx],
			na.rm = TRUE
		)
		expect_true(
			coverage_cls >= 0.80,
			info = paste("Class", cls, "coverage was", coverage_cls)
		)
	}
})

# ============================================================
# 4. Warns about unseen classes in pred
# ============================================================

test_that("warns about classes in pred not in calib", {
	pred_class_new <- class_test
	pred_class_new[1:3] <- "D"  # class D not in calibration
	expect_warning(
		pinterval_mondrian(
			pred = pred_test, pred_class = pred_class_new,
			calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
			alpha = 0.1
		),
		"class"
	)
})
