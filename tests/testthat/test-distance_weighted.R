# Tests for distance-weighted conformal prediction functionality

# --- Setup ---
set.seed(42)
n <- 400
x1 <- runif(n)
x2 <- runif(n)
y <- rnorm(n, mean = x1 + x2, sd = 0.5)

pred_cal <- (x1 + x2)[1:200]
truth_cal <- y[1:200]
feat_cal <- cbind(x1[1:200], x2[1:200])
pred_test <- (x1 + x2)[201:400]
truth_test <- y[201:400]
feat_test <- cbind(x1[201:400], x2[201:400])

# ============================================================
# 1. Distance-weighted conformal
# ============================================================

test_that("distance-weighted conformal runs with euclidean", {
	result <- pinterval_conformal(
		pred = pred_test[1:5],
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1,
		distance_weighted_cp = TRUE,
		distance_features_calib = feat_cal,
		distance_features_pred = feat_test[1:5, ],
		distance_type = "euclidean"
	)
	expect_s3_class(result, "tbl_df")
	expect_equal(nrow(result), 5)
})

test_that("distance-weighted conformal runs with mahalanobis", {
	result <- pinterval_conformal(
		pred = pred_test[1:5],
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1,
		distance_weighted_cp = TRUE,
		distance_features_calib = feat_cal,
		distance_features_pred = feat_test[1:5, ],
		distance_type = "mahalanobis"
	)
	expect_s3_class(result, "tbl_df")
	expect_equal(nrow(result), 5)
})

test_that("distance-weighted conformal errors without features", {
	expect_error(
		pinterval_conformal(
			pred = pred_test[1:5],
			calib = pred_cal,
			calib_truth = truth_cal,
			distance_weighted_cp = TRUE,
			distance_features_calib = NULL,
			distance_features_pred = feat_test[1:5, ]
		),
		"must be provided"
	)
})

test_that("distance normalization options work", {
	for (norm in c("none", "minmax", "sd")) {
		result <- pinterval_conformal(
			pred = pred_test[1:3],
			calib = pred_cal,
			calib_truth = truth_cal,
			alpha = 0.1,
			distance_weighted_cp = TRUE,
			distance_features_calib = feat_cal,
			distance_features_pred = feat_test[1:3, ],
			distance_type = "euclidean",
			normalize_distance = norm
		)
		expect_s3_class(result, "tbl_df")
	}
})

test_that("all weight functions work", {
	for (wf in c(
		"gaussian_kernel",
		"caucy_kernel",
		"logistic",
		"reciprocal_linear"
	)) {
		result <- pinterval_conformal(
			pred = pred_test[1:3],
			calib = pred_cal,
			calib_truth = truth_cal,
			alpha = 0.1,
			distance_weighted_cp = TRUE,
			distance_features_calib = feat_cal,
			distance_features_pred = feat_test[1:3, ],
			distance_type = "euclidean",
			weight_function = wf
		)
		expect_s3_class(result, "tbl_df")
	}
})

# ============================================================
# 2. Distance-weighted bootstrap
# ============================================================

test_that("distance-weighted bootstrap runs", {
	result <- pinterval_bootstrap(
		pred = pred_test[1:5],
		calib = pred_cal,
		calib_truth = truth_cal,
		alpha = 0.1,
		distance_weighted_bootstrap = TRUE,
		distance_features_calib = feat_cal,
		distance_features_pred = feat_test[1:5, ],
		distance_type = "euclidean"
	)
	expect_s3_class(result, "tbl_df")
	expect_equal(nrow(result), 5)
})

# ============================================================
# 3. Dimension mismatch errors
# ============================================================

test_that("mismatched feature dimensions error", {
	expect_error(
		pinterval_conformal(
			pred = pred_test[1:3],
			calib = pred_cal,
			calib_truth = truth_cal,
			distance_weighted_cp = TRUE,
			distance_features_calib = feat_cal,
			distance_features_pred = feat_test[1:3, 1, drop = FALSE] # only 1 column vs 2
		),
		"same number of columns"
	)
})

test_that("wrong number of rows in features errors", {
	expect_error(
		pinterval_conformal(
			pred = pred_test[1:3],
			calib = pred_cal,
			calib_truth = truth_cal,
			distance_weighted_cp = TRUE,
			distance_features_calib = feat_cal[1:10, ], # 10 rows vs 200 calib
			distance_features_pred = feat_test[1:3, ]
		),
		"rows.*match"
	)
})
