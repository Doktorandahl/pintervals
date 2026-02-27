#' Conformal Prediction Intervals of Continuous Values
#'
#' @description
#' This function calculates conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction. The intervals are computed using either a calibration set with predicted and true values or a set of pre-computed non-conformity scores from the calibration set. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals.
#'
#' @param pred Vector of predicted values
#' @param calib  A numeric vector of predicted values in the calibration partition, or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values. If calib is a numeric vector, calib_truth must be provided.
#' @param calib_truth A numeric vector of true values in the calibration partition. Only required if calib is a numeric vector
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param ncs_type A string specifying the type of nonconformity score to use. Available options are:
#' \itemize{
#'   \item \code{"absolute_error"}: \eqn{|y - \hat{y}|}
#'   \item \code{"relative_error"}: \eqn{|y - \hat{y}| / \hat{y}}
#'   \item \code{"zero_adjusted_relative_error"}: \eqn{|y - \hat{y}| / (\hat{y} + 1)}
#'   \item \code{"heterogeneous_error"}: \eqn{|y - \hat{y}| / \sigma_{\hat{y}}} absolute error divided by a measure of heteroskedasticity, computed as the predicted value from a linear model of the absolute error on the predicted values
#'   \item \code{"raw_error"}: the signed error \eqn{y - \hat{y}}
#' }
#' The default is \code{"absolute_error"}.
#' @param grid_size The number of points to use in the grid search between the lower and upper bound. Default is 10,000. A larger grid size increases the resolution of the prediction intervals but also increases computation time.
#' @param resolution Alternatively to grid_size. The minimum step size between grid points. Useful if the a specific resolution is desired. Default is NULL.
#'
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used. Primarily useful when the possible outcome values are outside the range of values observed in the calibration set. If not provided, the minimum (true) value of the calibration partition will be used.
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used. Primarily useful when the possible outcome values are outside the range of values observed in the calibration set. If not provided, the maximum (true) value of the calibration partition will be used.
#' @param distance_weighted_cp Logical. If \code{TRUE}, weighted conformal prediction is performed where the non-conformity scores are weighted based on the distance between calibration and prediction points in feature space. Default is \code{FALSE}. See details for more information.
#'
#' @param distance_features_calib A matrix, data frame, or numeric vector of features from which to compute distances when \code{distance_weighted_cp = TRUE}. This should contain the feature values for the calibration set. Must have the same number of rows as the calibration set. Can be the predicted values themselves, or any other features which give a meaningful distance measure.
#' @param distance_features_pred A matrix, data frame, or numeric vector of feature values for the prediction set. Must be the same features as specified in \code{distance_features_calib}. Required if \code{distance_weighted_cp = TRUE}.
#'
#' @param distance_type The type of distance metric to use when computing distances between calibration and prediction points. Options are 'mahalanobis' (default) and 'euclidean'.
#'
#' @param normalize_distance Either 'minmax', 'sd', or 'none'. Indicates if and how to normalize the distances when distance_weighted_cp is TRUE. Normalization helps ensure that distances are on a comparable scale across features. Default is 'none'.
#'
#' @param weight_function A character string specifying the weighting kernel to use for distance-weighted conformal prediction. Options are:
#' \itemize{
#'   \item \code{"gaussian_kernel"}: \eqn{ w(d) = e^{-d^2} }
#'   \item \code{"caucy_kernel"}: \eqn{ w(d) = 1/(1 + d^2) }
#'   \item \code{"logistic"}: \eqn{ w(d) = 1//(1 + e^{d}) }
#'   \item \code{"reciprocal_linear"}: \eqn{ w(d) = 1/(1 + d) }
#' }
#' The default is \code{"gaussian_kernel"}. Distances are computed as the Euclidean distance between the calibration and prediction feature vectors.

#' @details
#' This function computes prediction intervals using inductive conformal prediction. The calibration set must include predicted values and true values. These can be provided either as separate vectors (`calib`and `calib_truth`) or as a two-column tibble or matrix where the first column contains the predicted values and the second column contains the true values. If `calib` is a numeric vector, `calib_truth` must also be provided.
#'
#' Non-conformity scores are calculated using the specified `ncs_type`, which determines how the prediction error is measured. Available options include:
#'
#' - `"absolute_error"`: the absolute difference between predicted and true values.
#' - `"relative_error"`: the absolute error divided by the true value.
#' - `"za_relative_error"`: zero-adjusted relative error, which replaces small or zero true values with a small constant to avoid division by zero.
#' - `"heterogeneous_error"`: absolute error scaled by a linear model of prediction error magnitude as a function of the predicted value.
#' - `"raw_error"`: the signed difference between predicted and true values.
#'
#' These options provide flexibility to adapt to different patterns of prediction error across the outcome space.
#'
#' To determine the prediction intervals, the function performs a grid search over a specified range of possible outcome values, identifying intervals that satisfy the desired confidence level of \eqn{1 - \alpha}. The user can define the range via the `lower_bound` and `upper_bound` parameters. If these are not supplied, the function defaults to using the minimum and maximum of the true values in the calibration data.
#'
#' The resolution of the grid search can be controlled by either the `resolution` argument, which sets the minimum step size, or the `grid_size` argument, which sets the number of grid points. For wide prediction spaces, the grid search may be computationally intensive. In such cases, increasing the `resolution` or reducing the `grid_size` may improve performance.
#'
#' When `distance_weighted_cp = TRUE`, the function applies distance-weighted conformal prediction, which adjusts the influence of calibration non-conformity scores based on how similar each calibration point is to the target prediction.  This approach preserves the distribution-free nature of conformal prediction while allowing intervals to adapt to local patterns, often yielding tighter and more responsive prediction sets in heterogeneous data environments.
#'
#' Distances are computed between the feature matrices or vectors supplied via `distance_features_calib` and `distance_features_pred`. These distances are then transformed into weights using the selected kernel in `weight_function`, with rapidly decaying kernels (e.g., Gaussian) emphasizing strong locality and slower decays (e.g., reciprocal or Cauchy) providing smoother influence. Distances can be geographic coordinates, predicted values, or any other relevant features that capture similarity in the context of the prediction task. The distance metric is specified via `distance_type`, with options for Mahalanobis or Euclidean distance. The default is Mahalanobis distance, which accounts for correlations between features. Normalization of distances can be applied using the `normalize_distance` parameter. Normalization is primarily useful for euclidean distances to ensure that features on different scales do not disproportionately influence the distance calculations.
#'
#'
#'
#'
#' @return A tibble with the predicted values and the lower and upper bounds of the prediction intervals.
#'
#' @export
#'
#' @examples
#'
#' # Generate example data
#' library(dplyr)
#' library(tibble)
#' x1 <- runif(1000)
#' x2 <- runif(1000)
#' y <- rlnorm(1000, meanlog = x1 + x2, sdlog = 0.5)
#' df <- tibble(x1, x2, y)
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model to the training data
#' mod <- lm(log(y) ~ x1 + x2, data=df_train)
#'
#' # Generate predictions on the original y scale for the calibration data
#' calib_pred  <- exp(predict(mod, newdata=df_cal))
#' calib_truth <- df_cal$y
#'
#' # Generate predictions for the test data
#' pred_test <- exp(predict(mod, newdata=df_test))
#'
#' # Calculate prediction intervals using conformal prediction.
#' pinterval_conformal(pred_test,
#' calib = calib_pred,
#' calib_truth = calib_truth,
#' alpha = 0.1,
#' lower_bound = 0)
#'
pinterval_conformal <- function(
	pred,
	calib = NULL,
	calib_truth = NULL,
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
	# --- Input validation ---

	if (!is.numeric(pred)) {
		stop(
			"pinterval_conformal: 'pred' must be a numeric vector.",
			call. = FALSE
		)
	}

	if (any(is.na(pred))) {
		warning(
			"pinterval_conformal: 'pred' contains NA values. Corresponding intervals will be NA.",
			call. = FALSE
		)
	}

	if (is.null(calib)) {
		stop(
			"pinterval_conformal: 'calib' must be provided.",
			call. = FALSE
		)
	}

	if (is.numeric(calib) && is.null(calib_truth)) {
		stop(
			"pinterval_conformal: 'calib_truth' must be provided when 'calib' is a numeric vector.",
			call. = FALSE
		)
	}

	if (!is.numeric(calib) && !is.matrix(calib) && !is.data.frame(calib)) {
		stop(
			"pinterval_conformal: 'calib' must be a numeric vector, matrix, or data frame.",
			call. = FALSE
		)
	}

	if (!is.numeric(calib) && (is.matrix(calib) || is.data.frame(calib)) && ncol(calib) != 2) {
		stop(
			"pinterval_conformal: 'calib' must be a numeric vector or a 2-column matrix/data frame with predicted values in column 1 and true values in column 2.",
			call. = FALSE
		)
	}

	if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
		stop(
			"pinterval_conformal: 'alpha' must be a single numeric value between 0 and 1 (exclusive).",
			call. = FALSE
		)
	}

	if (!is.null(grid_size) && (!is.numeric(grid_size) || length(grid_size) != 1 || grid_size < 1)) {
		stop(
			"pinterval_conformal: 'grid_size' must be a single positive integer.",
			call. = FALSE
		)
	}

	if (
		is.null(grid_size) &&
			(!is.numeric(resolution) ||
				resolution <= 0 ||
				length(resolution) != 1)
	) {
		stop(
			"pinterval_conformal: when 'grid_size' is NULL, 'resolution' must be a single positive numeric value.",
			call. = FALSE
		)
	}

	if (!is.null(grid_size) && !is.null(resolution)) {
		warning(
			"pinterval_conformal: both 'grid_size' and 'resolution' are provided; 'resolution' will be ignored.",
			call. = FALSE
		)
	}

	if (!is.null(lower_bound) && !is.null(upper_bound) && lower_bound >= upper_bound) {
		stop(
			"pinterval_conformal: 'lower_bound' must be strictly less than 'upper_bound'.",
			call. = FALSE
		)
	}

	normalize_distance <- match.arg(normalize_distance, c('none', 'minmax', 'sd'))

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

	distance_type <- match.arg(
		distance_type,
		c('mahalanobis', 'euclidean')
	)

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

	if (length(calib) != length(calib_truth)) {
		stop(
			"pinterval_conformal: 'calib' and 'calib_truth' must have the same length (got ",
			length(calib), " and ", length(calib_truth), ").",
			call. = FALSE
		)
	}

	if (any(is.na(calib)) || any(is.na(calib_truth))) {
		warning(
			"pinterval_conformal: 'calib' or 'calib_truth' contains NA values. This may affect interval estimation.",
			call. = FALSE
		)
	}

	if (ncs_type == 'heterogeneous_error') {
		coefs <- stats::coef(stats::lm(abs(calib - calib_truth) ~ calib))
	} else {
		coefs <- NULL
	}

	if (distance_weighted_cp) {
		validate_distance_inputs(
			distance_features_calib, distance_features_pred,
			length(calib), length(pred),
			fn_name = "pinterval_conformal"
		)
		distance_features_calib <- as.matrix(distance_features_calib)
		distance_features_pred <- as.matrix(distance_features_pred)
		weight_function <- resolve_weight_function(weight_function)
	}

	ncs <- ncs_compute(type = ncs_type, calib, calib_truth, coefs = coefs)

	if (is.null(lower_bound)) {
		lower_bound <- min(calib_truth)
	}
	if (is.null(upper_bound)) {
		upper_bound <- max(calib_truth)
	}

	cp_set <- grid_finder(
		y_min = lower_bound,
		y_max = upper_bound,
		ncs = ncs,
		y_hat = pred,
		min_step = resolution,
		alpha = alpha,
		grid_size = grid_size,
		ncs_type = ncs_type,
		coefs = coefs,
		calib = calib,
		distance_weighted_cp = distance_weighted_cp,
		distance_features_calib = distance_features_calib,
		distance_features_pred = distance_features_pred,
		distance_type = distance_type,
		normalize_distance = normalize_distance,
		weight_function = weight_function
	)

	return(cp_set)
}
