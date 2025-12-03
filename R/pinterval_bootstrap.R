#' Bootstrap prediction intervals
#'
#' @description
#' This function computes bootstrapped prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using bootstrapped prediction errors. The prediction errors to bootstrap from are computed using either a calibration set with predicted and true values or a set of pre-computed prediction errors from a calibration dataset or other data which the model was not trained on (e.g. OOB errors from a model using bagging). The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals.
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values
#' @param calib_truth A numeric vector of true values in the calibration partition. Only required if calib is a numeric vector
#' @param error_type The type of error to use for the prediction intervals. Can be 'raw' or 'absolute'. If 'raw', bootstrapping will be done on the raw prediction errors. If 'absolute', bootstrapping will be done on the absolute prediction errors with random signs. Default is 'raw'
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param n_bootstraps The number of bootstraps to perform. Default is 1000
#' @param distance_weighted_bootstrap Logical. If TRUE, the function will use distance-weighted bootstrapping. Default is FALSE. If TRUE, the probability of selecting a prediction error is weighted by the distance to the predicted value using the specified distance function and weight function. If FALSE, standard bootstrapping is performed.
#' @param distance_features_calib A matrix, data frame, or numeric vector of features from which to compute distances when \code{distance_weighted_cp = TRUE}. This should contain the feature values for the calibration set. Must have the same number of rows as the calibration set. Can be the predicted values themselves, or any other features which give a meaningful distance measure.
#' @param distance_features_pred A matrix, data frame, or numeric vector of feature values for the prediction set. Must be the same features as specified in \code{distance_features_calib}. Required if \code{distance_weighted_cp = TRUE}.
#' @param distance_type The type of distance metric to use when computing distances between calibration and prediction points. Options are 'mahalanobis' (default) and 'euclidean'.
#'
#' @param normalize_distance Either 'minmax', 'sd', or 'none'. Indicates if and how to normalize the distances when distance_weighted_bootstrap is TRUE. Normalization helps ensure that distances are on a comparable scale across features. Default is 'minmax'.
#'
#' @param weight_function A character string specifying the weighting kernel to use for distance-weighted conformal prediction. Options are:
#' \itemize{
#'   \item \code{"gaussian_kernel"}: \eqn{ w(d) = e^{-d^2} }
#'   \item \code{"caucy_kernel"}: \eqn{ w(d) = 1/(1 + d^2) }
#'   \item \code{"logistic"}: \eqn{ w(d) = 1//(1 + e^{d}) }
#'   \item \code{"reciprocal_linear"}: \eqn{ w(d) = 1/(1 + d) }
#' }
#' The default is \code{"gaussian_kernel"}. Distances are computed as the Euclidean distance between the calibration and prediction feature vectors.
#' @param calibrate = FALSE Logical. If TRUE, the function will calibrate the predictions and intervals using the calibration set. Default is FALSE. See details for more information.
#' @param calibration_method The method to use for calibration. Can be "glm" or "isotonic". Default is "glm". Only used if calibrate = TRUE.
#' @param calibration_family The family used for the calibration model. Default is "gaussian". Only used if calibrate = TRUE and calibration_method = "glm".
#' @param calibration_transform Optional transformation to apply to the predictions before calibration. Default is NULL. Only used if calibrate = TRUE and calibration_method = "glm".
#'
#' #' @details
#' This function estimates prediction intervals using bootstrapped prediction errors derived from a calibration set. It supports both standard and distance-weighted bootstrapping. The calibration set must consist of predicted values and corresponding true values, either provided as separate vectors or as a two-column tibble or matrix. Alternatively, users may provide a vector of precomputed prediction errors if model predictions and truths are already processed.
#'
#' Two types of error can be used for bootstrapping:
#' - `"raw"`: bootstrapping is performed on the raw signed prediction errors (truth - prediction), allowing for asymmetric prediction intervals.
#' - `"absolute"`: bootstrapping is done on the absolute errors, and random signs are applied when constructing intervals. This results in (approximately) symmetric intervals around the prediction.
#'
#' Distance-weighted bootstrapping (`dw_bootstrap = TRUE`) can be used to give more weight to calibration errors closer to each test prediction. Distances are calculated using the specified `distance_function`, which defaults to the absolute difference. Distance weighted bootstrapping assigns a probability of selection for each calibration error based on its distance to the predicted value, computed as \eqn{1/(1 + distance)}. The distance function can be a custom function or a character string that matches a function, which takes two arguments: a single predicted value and a vector of calibration values, returning a numeric vector of distances.
#'
#' The number of bootstrap samples is controlled via the `n_bootstraps` parameter. For computational efficiency, this can be reduced at the cost of interval precision.
#'
#' Optionally, the predicted values can be calibrated in conjuncture with interval construction by setting `calibrate = TRUE`. In this case, the predictions are passed through `calibrate_predictions()` to adjust the predictions based on the calibration set. The calibration method can be specified using `calibration_method` and `calibration_family`, with "glm" being the default method. See \link[pintervals]{calibrate_predictions} for more information on calibration.
#'
#' @return A tibble with the predicted values, lower bounds, and upper bounds of the prediction intervals
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#'
#' # Simulate some data
#' set.seed(42)
#' x1 <- runif(1000)
#' x2 <- runif(1000)
#' y <- rlnorm(1000, meanlog = x1 + x2, sdlog = 0.4)
#' df <- tibble(x1, x2, y)
#'
#' # Split into train/calibration/test
#' df_train <- df[1:500, ]
#' df_cal <- df[501:750, ]
#' df_test <- df[751:1000, ]
#'
#' # Fit a log-linear model
#' model <- lm(log(y) ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' pred_cal <- exp(predict(model, newdata = df_cal))
#' pred_test <- exp(predict(model, newdata = df_test))
#'
#' # Compute bootstrap prediction intervals
#' intervals <- pinterval_bootstrap(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   error_type = "raw",
#'   alpha = 0.1,
#'   n_bootstraps = 1000
#' )
#'
pinterval_bootstrap <- function(
	pred,
	calib,
	calib_truth = NULL,
	error_type = c('raw', 'absolute'),
	alpha = 0.1,
	n_bootstraps = 1000,
	distance_weighted_bootstrap = FALSE,
	distance_features_calib = NULL,
	distance_features_pred = NULL,
	distance_type = c('mahalanobis', 'euclidean'),
	normalize_distance = TRUE,
	weight_function = c(
		'gaussian_kernel',
		'caucy_kernel',
		'logistic',
		'reciprocal_linear'
	),
	calibrate = FALSE,
	calibration_method = c('glm', 'isotonic'),
	calibration_family = 'gaussian',
	calibration_transform = NULL
) {
	i <- NA
	if (!is.numeric(pred)) {
		stop('pred must be a single number or a numeric vector')
	}

	if (is.numeric(calib) & (is.null(calib_truth))) {
		stop('If calib is numeric, calib_truth must be provided')
	}

	if (!is.numeric(calib) && ncol(calib) != 2) {
		stop(
			'calib must be a numeric vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values'
		)
	}

	if (calibrate && distance_weighted_bootstrap) {
		warning(
			'calibrate is TRUE, but dw_bootstrap is also TRUE. Calibration will be performed after bootstrapping, and the corresponding prediction intervals will be adjusted accordingly.'
		)
	}

	if (!is.numeric(calib)) {
		calib_org <- calib
		if (is.matrix(calib)) {
			calib <- as.numeric(calib_org[, 1])
			calib_truth <- as.numeric(calib_org[, 2])
		} else {
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
		}
	}

	if (distance_weighted_bootstrap) {
		if (is.null(distance_features_calib) || is.null(distance_features_pred)) {
			stop(
				'If distance_weighted_cp is TRUE, distance_features_calib and distance_features_pred must be provided'
			)
		}
		if (
			!is.matrix(distance_features_calib) &&
				!is.data.frame(distance_features_calib) &&
				!is.numeric(distance_features_calib)
		) {
			stop(
				'distance_features_calib must be a matrix, data frame, or numeric vector'
			)
		}
		if (
			!is.matrix(distance_features_pred) &&
				!is.data.frame(distance_features_pred) &&
				!is.numeric(distance_features_pred)
		) {
			stop(
				'distance_features_pred must be a matrix, data frame, or numeric vector'
			)
		}
		if (
			is.numeric(distance_features_calib) && is.numeric(distance_features_pred)
		) {
			if (
				length(distance_features_calib) != length(calib) ||
					length(distance_features_pred) != length(pred)
			) {
				stop(
					'If distance_features_calib and distance_features_pred are numeric vectors, they must have the same length as calib and pred, respectively'
				)
			}
		} else if (
			is.matrix(distance_features_calib) ||
				is.data.frame(distance_features_calib)
		) {
			if (nrow(distance_features_calib) != length(calib)) {
				stop(
					'If distance_features_calib is a matrix or data frame, it must have the same number of rows as calib'
				)
			}
			if (ncol(distance_features_calib) != ncol(distance_features_pred)) {
				stop(
					'distance_features_calib and distance_features_pred must have the same number of columns'
				)
			}
			if (nrow(distance_features_pred) != length(pred)) {
				stop(
					'If distance_features_pred is a matrix or data frame, it must have the same number of rows as pred'
				)
			}
		}

		distance_type <- match.arg(distance_type, c('mahalanobis', 'euclidean'))

		distance_features_calib <- as.matrix(distance_features_calib)
		distance_features_pred <- as.matrix(distance_features_pred)

		if (!is.function(weight_function)) {
			weight_function <- match.arg(
				weight_function,
				c('gaussian_kernel', 'caucy_kernel', 'logistic', 'reciprocal_linear')
			)
			weight_function <- switch(
				weight_function,
				'gaussian_kernel' = function(d) exp(-d^2),
				'caucy_kernel' = function(d) 1 / (1 + d^2),
				'logistic' = function(d) 1 / (1 + exp(d)),
				'reciprocal_linear' = function(d) 1 / (1 + d)
			)
		}
	}

	error_type <- match.arg(error_type, c('raw', 'absolute'))

	if (error_type == 'raw') {
		error <- calib - calib_truth
	} else if (error_type == 'absolute') {
		error <- abs(calib - calib_truth)
	}

	boot_set <- foreach::foreach(i = 1:length(pred)) %do%
		bootstrap_inner(
			pred = pred[i],
			calib = calib,
			error = error,
			nboot = n_bootstraps,
			alpha = alpha,
			error_type = error_type,
			distance_weighted_bootstrap = distance_weighted_bootstrap,
			distance_features_calib = distance_features_calib,
			distance_features_pred = distance_features_pred,
			distance_type = distance_type,
			normalize_distance = normalize_distance,
			weight_function = weight_function
		)

	boot_set <- dplyr::bind_rows(boot_set)

	if (calibrate) {
		calibration_method <- match.arg(calibration_method, c('glm', 'isotonic'))
		if (calibration_method == 'glm') {
			calibrated_predictions <- calibrate_predictions(
				pred,
				calib,
				calib_truth,
				method = calibration_method,
				family = calibration_family,
				transform = calibration_transform
			)
		} else if (calibration_method == 'isotonic') {
			calibrated_predictions <- calibrate_predictions(
				pred,
				calib,
				calib_truth,
				method = calibration_method
			)
		}

		calibrated_diffs <- pred - calibrated_predictions
		boot_set$lower_bound <- boot_set$lower_bound + calibrated_diffs
		boot_set$upper_bound <- boot_set$upper_bound + calibrated_diffs
		boot_set$pred <- calibrated_predictions
	}

	return(boot_set)
}
