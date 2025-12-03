#' Mondrian conformal prediction intervals for continuous predictions
#'
#'@description
#'This function calculates Mondrian conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction on a Mondrian class-by-class basis. The intervals are computed using a calibration set with predicted and true values and their associated classes. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals. Mondrian conformal prediction intervals are useful when the prediction error is not constant across groups or classes, as they allow for locally valid coverage by ensuring that the coverage level \eqn{1 - \eqn{\alpha}} holds within each class—assuming exchangeability of non-conformity scores within classes.
#'
#' @param pred Vector of predicted values or a 2 column tibble or matrix with the first column being the predicted values and the second column being the Mondrian class labels. If pred is a numeric vector, pred_class must be provided.
#' @param pred_class A vector of class identifiers for the predicted values. This is used to group the predictions by class for Mondrian conformal prediction.
#' @param calib A numeric vector of predicted values in the calibration partition or a 3 column tibble or matrix with the first column being the predicted values and the second column being the truth values and the third column being the Mondrian class labels. If calib is a numeric vector, calib_truth and calib_class must be provided.
#' @param calib_truth A numeric vector of true values in the calibration partition
#' @param calib_class A vector of class identifiers for the calibration set.
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
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#' @param distance_weighted_cp Logical. If \code{TRUE}, weighted conformal prediction is performed where the non-conformity scores are weighted based on the distance between calibration and prediction points in feature space. Default is \code{FALSE}.
#'
#' @param distance_features_calib A matrix, data frame, or numeric vector of features from which to compute distances when \code{distance_weighted_cp = TRUE}. This should contain the feature values for the calibration set. Must have the same number of rows as the calibration set. Can be the predicted values themselves, or any other features which give a meaningful distance measure.
#' @param distance_features_pred A matrix, data frame, or numeric vector of feature values for the prediction set. Must be the same features as specified in \code{distance_features_calib}. Required if \code{distance_weighted_cp = TRUE}.
#' @param distance_type The type of distance metric to use when computing distances between calibration and prediction points. Options are 'mahalanobis' (default) and 'euclidean'.
#'
#' @param normalize_distance Either 'minmax', 'sd', or 'none'. Indicates if and how to normalize the distances when distance_weighted_cp is TRUE. Normalization helps ensure that distances are on a comparable scale across features. Default is 'minmax'.
#'
#' @param weight_function A character string specifying the weighting kernel to use for distance-weighted conformal prediction. Options are:
#' \itemize{
#'   \item \code{"gaussian_kernel"}: \eqn{ w(d) = e^{-d^2} }
#'   \item \code{"caucy_kernel"}: \eqn{ w(d) = 1/(1 + d^2) }
#'   \item \code{"logistic"}: \eqn{ w(d) = 1//(1 + e^{d}) }
#'   \item \code{"reciprocal_linear"}: \eqn{ w(d) = 1/(1 + d) }
#' }
#' The default is \code{"gaussian_kernel"}. Distances are computed as the Euclidean distance between the calibration and prediction feature vectors.
#' @param resolution The minimum step size for the grid search. Default is 0.01. See details for more information.
#' @param grid_size Alternative to `resolution`, the number of points to use in the grid search between the lower and upper bound. If provided, resolution will be ignored.
#'
#' @return A tibble with the predicted values, the lower and upper bounds of the prediction intervals. If treat_noncontiguous is 'non_contiguous', the lower and upper bounds are set in a list variable called 'intervals' where all non-contiguous intervals are stored.
#'
#' @details
#' This function computes Mondrian conformal prediction intervals using inductive conformal prediction applied separately within each class (also called strata or groups) of the calibration data. It is especially useful when prediction error varies systematically across known categories, allowing for class-conditional validity by ensuring that the prediction intervals attain the desired coverage level \eqn{1 - \eqn{\alpha}} within each class—under the assumption of exchangeability within classes.
#'
#' The calibration set must include predicted values, true values, and corresponding class labels. These can be supplied as separate vectors (`calib`, `calib_truth`, and `calib_class`) or as a single three-column matrix or tibble.
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
#' To determine the prediction intervals, the function performs a grid search over a specified range of possible outcome values, identifying intervals that satisfy the desired confidence level of \eqn{1 - \eqn{\alpha}}. The user can define the range via the `lower_bound` and `upper_bound` parameters. If these are not supplied, the function defaults to using the minimum and maximum of the true values in the calibration data.
#'
#' The resolution of the grid search can be controlled by either the `resolution` argument, which sets the minimum step size, or the `grid_size` argument, which sets the number of grid points. For wide prediction spaces, the grid search may be computationally intensive. In such cases, increasing the `resolution` or reducing the `grid_size` may improve performance.
#'
#'
#' The function returns a tibble with the original predictions and their corresponding lower and upper prediction interval bounds.
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
	lower_bound = NULL,
	upper_bound = NULL,
	alpha = 0.1,
	ncs_type = c(
		'absolute_error',
		'relative_error',
		'za_relative_error',
		'heterogeneous_error',
		'raw_error'
	),
	distance_weighted_cp = FALSE,
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
	calibration_method = 'glm',
	calibration_family = NULL,
	calibration_transform = NULL,
	resolution = 0.01,
	grid_size = NULL
) {
	i <- NA

	if (setdiff(unique(pred_class), unique(calib_class)) %>% length() > 0) {
		warning(
			'Some classes in pred_class are not present in calib_class. These will result in NA prediction intervals for those classes.'
		)
	}

	if (!is.numeric(pred) && ncol(pred) != 2) {
		stop(
			'pred must be a numeric scalar or vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the class labels'
		)
	}

	if (is.numeric(pred) && is.null(pred_class)) {
		stop('If pred is numeric, pred_class must be provided')
	}

	if (!is.numeric(pred)) {
		pred_class <- as.numeric(pred[[2]])
		pred <- as.numeric(pred[[1]])
	}

	if (is.numeric(calib) & is.null(calib_truth)) {
		stop('If calib is numeric, calib_truth must be provided')
	}

	if (is.numeric(calib) & is.null(calib_class)) {
		stop('If calib is numeric, calib_class must be provided')
	}
	if (!is.numeric(calib) && ncol(calib) < 3) {
		stop(
			'calib must be a numeric vector or a 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and the third column being the class labels'
		)
	}

	if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha) != 1) {
		stop('alpha must be a single numeric value between 0 and 1')
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

	if (!is.numeric(calib)) {
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

	if (ncs_type == 'heterogeneous_error') {
		coefs <- stats::coef(stats::lm(abs(calib - calib_truth) ~ calib))
	} else {
		coefs <- NULL
	}

	if (distance_weighted_cp) {
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

		distance_features_calib <- as.matrix(distance_features_calib)
		distance_features_pred <- as.matrix(distance_features_pred)
		distance_type <- match.arg(distance_type, c('mahalanobis', 'euclidean'))

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

	nobs_class <- as.numeric(table(calib_class))

	if (any(nobs_class * alpha / 2 < 1)) {
		warning(
			'Some classes have too few observations to calculate prediction intervals at the specified alpha level. Consider using a larger calibration set or a higher alpha level'
		)
	}

	class_labels <- sort(unique(pred_class))
	if (length(class_labels) < 2) {
		stop(
			'Calibration set must have at least two classes For continuous prediction intervals without classes, use pinterval_conformal() instead of pinterval_mondrian()'
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
					calibrate = calibrate,
					calibration_method = calibration_method,
					calibration_family = calibration_family,
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
