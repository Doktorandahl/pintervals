#' Mondrian Conformal Prediction Intervals for Continuous Predictions
#'
#'@description
#'This function calculates Mondrian conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction on a Mondrian class-by-class basis. The intervals are computed using a calibration set with predicted and true values and their associated classes. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals. Mondrian conformal prediction intervals are useful when the prediction error is not constant across groups or classes, as they allow for locally valid coverage by ensuring that the coverage level \eqn{1 - \alpha} holds within each class—assuming exchangeability of non-conformity scores within classes.
#'
#' @inheritParams pinterval_conformal
#' @param pred_class A vector of class identifiers for the predicted values. This is used to group the predictions by class for Mondrian conformal prediction.
#' @param calib_class A vector of class identifiers for the calibration set.
#'
#' @return A tibble with predicted values, lower and upper prediction interval bounds, and class labels.
#'
#' @details
#' `pinterval_mondrian()` extends [pinterval_conformal()] to the Mondrian
#' setting, where prediction intervals are calibrated separately within
#' user-defined groups (often called "Mondrian categories"). Instead of
#' pooling all calibration residuals into a single reference distribution,
#' the method constructs a separate non-conformity distribution for each
#' subgroup defined by a grouping variable (e.g., region, regime type, or
#' income category). This allows the intervals to adapt to systematic
#' differences in error magnitude or variance across groups and targets
#' coverage conditional on group membership. It is especially useful when prediction error varies systematically across known categories, allowing for class-conditional validity by ensuring that the prediction intervals attain the desired coverage level \eqn{1 - \alpha} within each class—under the assumption of exchangeability within classes.
#'
#' Conceptually, the underlying  inductive conformal machinery is the same as in [pinterval_conformal()], but applied within groups rather than globally. For a detailed description of non-conformity scores, distance-weighting, and the general conformal  prediction framework, see [pinterval_conformal()].
#'
#' For `pinterval_mondrian()`, the calibration set must include predicted values, true values, and corresponding class labels. These can be supplied as separate vectors (`calib`, `calib_truth`, and `calib_class`) or as a single three-column matrix or tibble.
#'
#' @seealso \code{\link[pintervals]{pinterval_conformal}}
#'
#' @export
#'
#' @examples
#'
#' # Generate synthetic data
#' library(dplyr)
#' library(tibble)
#' set.seed(123)
#' x1 <- runif(1000)
#' x2 <- runif(1000)
#' group <- sample(c("A", "B", "C"), size = 1000, replace = TRUE)
#' mu <- ifelse(group == "A", 1 + x1 + x2,
#'       ifelse(group == "B", 2 + x1 + x2,
#'                         3 + x1 + x2))
#' y <- rlnorm(1000, meanlog = mu, sdlog = 0.4)
#'
#' df <- tibble(x1, x2, group, y)
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model to the training data
#' mod <- lm(log(y) ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' calib <- exp(predict(mod, newdata = df_cal))
#' calib_truth <- df_cal$y
#' calib_class <- df_cal$group
#'
#' pred_test <- exp(predict(mod, newdata = df_test))
#' pred_test_class <- df_test$group
#'
#' # Apply Mondrian conformal prediction
#' pinterval_mondrian(pred = pred_test,
#' pred_class = pred_test_class,
#' calib = calib,
#' calib_truth = calib_truth,
#' calib_class = calib_class,
#' alpha = 0.1)
#'
pinterval_mondrian = function(
	pred,
	pred_class = NULL,
	calib = NULL,
	calib_truth = NULL,
	calib_class = NULL,
	alpha = 0.1,
	ncs_type = c(
		'absolute_error',
		'relative_error',
		'za_relative_error',
		'heterogeneous_error',
		'raw_error'
	),
	lower_bound = NULL,
	upper_bound = NULL,
	grid_size = 10000,
	resolution = NULL,
	distance_weighted_cp = FALSE,
	distance_features_calib = NULL,
	distance_features_pred = NULL,
	distance_type = c('mahalanobis', 'euclidean'),
	normalize_distance = 'none',
	weight_function = c(
		'gaussian_kernel',
		'caucy_kernel',
		'logistic',
		'reciprocal_linear'
	)
) {
	i <- NA

	# Validate pred
	if (!is.numeric(pred) && !is.matrix(pred) && !is.data.frame(pred)) {
		stop(
			'pinterval_mondrian: pred must be a numeric scalar or vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the class labels',
			call. = FALSE
		)
	}

	if (!is.numeric(pred) && ncol(pred) != 2) {
		stop(
			'pinterval_mondrian: pred must be a numeric scalar or vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the class labels',
			call. = FALSE
		)
	}

	if (any(is.na(pred))) {
		warning(
			'pinterval_mondrian: pred contains NA values',
			call. = FALSE
		)
	}

	if (is.numeric(pred) && is.null(pred_class)) {
		stop(
			'pinterval_mondrian: If pred is numeric, pred_class must be provided',
			call. = FALSE
		)
	}

	if (!is.numeric(pred)) {
		pred_class <- as.numeric(pred[[2]])
		pred <- as.numeric(pred[[1]])
	}

	# Validate calib
	if (is.null(calib)) {
		stop(
			'pinterval_mondrian: calib must be provided',
			call. = FALSE
		)
	}

	if (!is.numeric(calib) && !is.matrix(calib) && !is.data.frame(calib)) {
		stop(
			'pinterval_mondrian: calib must be a numeric vector or a 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and the third column being the class labels',
			call. = FALSE
		)
	}

	if (!is.vector(calib) && ncol(calib) < 3) {
		stop(
			'pinterval_mondrian: calib must be a numeric vector or a 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and the third column being the class labels',
			call. = FALSE
		)
	}

	if (any(is.na(calib))) {
		warning(
			'pinterval_mondrian: calib contains NA values',
			call. = FALSE
		)
	}

	if (is.numeric(calib) && is.null(calib_truth)) {
		stop(
			'pinterval_mondrian: If calib is numeric, calib_truth must be provided',
			call. = FALSE
		)
	}

	if (is.numeric(calib) && is.null(calib_class)) {
		stop(
			'pinterval_mondrian: If calib is numeric, calib_class must be provided',
			call. = FALSE
		)
	}

	if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha) != 1) {
		stop(
			'pinterval_mondrian: alpha must be a single numeric value between 0 and 1',
			call. = FALSE
		)
	}

	ncs_type <- match.arg(
		ncs_type,
		c(
			'absolute_error',
			'relative_error',
			'za_relative_error',
			'heterogeneous_error',
			'raw_error'
		)
	)

	# Parse calib if it's a matrix or data frame
	if (!is.vector(calib)) {
		calib_org <- calib
		if (is.matrix(calib)) {
			calib <- as.numeric(calib_org[, 1])
			calib_truth <- as.numeric(calib_org[, 2])
			if (is.null(calib_class) && ncol(calib_org) == 3) {
				calib_class <- as.numeric(calib_org[, 3])
			}
		} else {
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
			if (is.null(calib_class) && ncol(calib_org) == 3) {
				calib_class <- as.numeric(calib_org[[3]])
			}
		}
	}

	# Validate lengths after parsing
	if (length(pred_class) != length(pred)) {
		stop(
			'pinterval_mondrian: pred_class must have the same length as pred',
			call. = FALSE
		)
	}

	if (length(calib_class) != length(calib)) {
		stop(
			'pinterval_mondrian: calib_class must have the same length as calib',
			call. = FALSE
		)
	}

	if (length(calib_truth) != length(calib)) {
		stop(
			'pinterval_mondrian: calib_truth must have the same length as calib',
			call. = FALSE
		)
	}

	# Check for classes in pred that are not in calib (AFTER validation)
	if (setdiff(unique(pred_class), unique(calib_class)) %>% length() > 0) {
		warning(
			'pinterval_mondrian: Some classes in pred_class are not present in calib_class. These will result in NA prediction intervals for those classes.',
			call. = FALSE
		)
	}

	if (ncs_type == 'heterogeneous_error') {
		coefs <- stats::coef(stats::lm(abs(calib - calib_truth) ~ calib))
	} else {
		coefs <- NULL
	}

	# Validate and normalize distance parameters
	normalize_distance <- match.arg(normalize_distance, c('none', 'minmax', 'sd'))

	if (distance_weighted_cp) {
		validate_distance_inputs(
			distance_features_calib,
			distance_features_pred,
			length(calib),
			length(pred),
			fn_name = "pinterval_mondrian"
		)
		distance_features_calib <- as.matrix(distance_features_calib)
		distance_features_pred <- as.matrix(distance_features_pred)
		distance_type <- match.arg(distance_type, c('mahalanobis', 'euclidean'))
		weight_function <- resolve_weight_function(weight_function)
	}

	nobs_class <- as.numeric(table(calib_class))

	if (any(nobs_class * alpha / 2 < 1)) {
		warning(
			'pinterval_mondrian: Some classes have too few observations to calculate prediction intervals at the specified alpha level. Consider using a larger calibration set or a higher alpha level',
			call. = FALSE
		)
	}

	class_labels <- sort(unique(pred_class))
	class_labels_calib <- sort(unique(calib_class))
	if (length(class_labels) < 2 || length(class_labels_calib) < 2) {
		stop(
			"pinterval_mondrian: calibration set must have at least two classes. For continuous prediction intervals without classes, use pinterval_conformal() instead.",
			call. = FALSE
		)
	}

	cp_intervals <- foreach::foreach(
		i = 1:length(class_labels),
		.final = dplyr::bind_rows
	) %do%
		{
			indices <- which(pred_class == class_labels[i])
			if (length(calib[calib_class == class_labels[i]]) == 0) {
				res <- tibble::tibble(
					pred = pred[pred_class == class_labels[i]],
					lower_bound = NA_real_,
					upper_bound = NA_real_,
					indices = which(pred_class == class_labels[i])
				)
			} else {
				res <- suppressWarnings(pinterval_conformal(
					pred = pred[pred_class == class_labels[i]],
					lower_bound = lower_bound,
					upper_bound = upper_bound,
					ncs_type = ncs_type,
					calib = calib[calib_class == class_labels[i]],
					calib_truth = calib_truth[calib_class == class_labels[i]],
					distance_weighted_cp = distance_weighted_cp,
					distance_features_calib = distance_features_calib[
						calib_class == class_labels[i],
					],
					distance_features_pred = distance_features_pred[
						pred_class == class_labels[i],
					],
					distance_type = distance_type,
					normalize_distance = normalize_distance,
					weight_function = weight_function,
					alpha = alpha,
					resolution = resolution,
					grid_size = grid_size
				))
			}
			res$indices <- indices
			res
		}

	cp_intervals2 <- cp_intervals %>%
		dplyr::arrange(indices) %>%
		dplyr::select(-indices) %>%
		dplyr::mutate(class = pred_class)

	return(cp_intervals2)
}
