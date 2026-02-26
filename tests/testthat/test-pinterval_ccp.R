# Tests for pinterval_ccp()

# --- Setup ---
set.seed(42)
n <- 600
x <- runif(n)
class_raw <- sample(1:6, n, replace = TRUE)
mu <- ifelse(class_raw %in% 1:2, 1 + x, ifelse(class_raw %in% 3:4, 2 + x, 3 + x))
y <- rnorm(n, mean = mu, sd = 0.5)

pred_cal <- mu[1:300]
truth_cal <- y[1:300]
class_cal <- factor(class_raw[1:300])
pred_test <- mu[301:600]
truth_test <- y[301:600]
class_test <- factor(class_raw[301:600])

# ============================================================
# 1. Input validation
# ============================================================

test_that("pred must be numeric", {
	expect_error(
		pinterval_ccp(
			pred = "a", pred_class = class_test,
			calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
			n_clusters = 3
		),
		"pinterval_ccp.*pred.*numeric"
	)
})

test_that("calib must be provided", {
	expect_error(
		pinterval_ccp(
			pred = pred_test, pred_class = class_test,
			calib = NULL, n_clusters = 3
		),
		"pinterval_ccp.*calib.*provided"
	)
})

test_that("n_clusters must be a positive integer", {
	expect_error(
		pinterval_ccp(
			pred = pred_test, pred_class = class_test,
			calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
			n_clusters = -1
		),
		"n_clusters"
	)
})

test_that("cluster_method must be valid", {
	expect_error(
		pinterval_ccp(
			pred = pred_test, pred_class = class_test,
			calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
			n_clusters = 3, cluster_method = "invalid"
		),
		"arg"
	)
})

# ============================================================
# 2. Output structure
# ============================================================

test_that("output is a tibble with expected columns", {
	result <- pinterval_ccp(
		pred = pred_test, pred_class = class_test,
		calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
		n_clusters = 3, alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
	expect_true(all(c("pred", "lower_bound", "upper_bound") %in% names(result)))
	expect_equal(nrow(result), length(pred_test))
})

# ============================================================
# 3. Clustering methods both work
# ============================================================

test_that("kmeans cluster method produces valid output", {
	result <- pinterval_ccp(
		pred = pred_test[1:20], pred_class = class_test[1:20],
		calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
		n_clusters = 3, cluster_method = "kmeans", alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
})

test_that("ks cluster method produces valid output", {
	result <- pinterval_ccp(
		pred = pred_test[1:20], pred_class = class_test[1:20],
		calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
		n_clusters = 3, cluster_method = "ks", alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
})

# ============================================================
# 4. Optimize clusters
# ============================================================

test_that("optimize_n_clusters works", {
	result <- pinterval_ccp(
		pred = pred_test[1:20], pred_class = class_test[1:20],
		calib = pred_cal, calib_truth = truth_cal, calib_class = class_cal,
		optimize_n_clusters = TRUE,
		min_n_clusters = 2, max_n_clusters = 4,
		alpha = 0.1
	)
	expect_s3_class(result, "tbl_df")
})
