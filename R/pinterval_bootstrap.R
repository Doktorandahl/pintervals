#' Bootstrap prediction intervals
#'
#' @description
#' This function computes bootstrapped prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using bootstrapped prediction errors. The prediction errors to bootstrap from are computed using either a calibration set with predicted and true values or a set of pre-computed prediction errors from a calibration dataset or other data which the model was not trained on (e.g. OOB errors from a model using bagging). The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals.
#'
#' @inheritParams pinterval_conformal
#' @param error_type The type of error to use for the prediction intervals. Can be 'raw' or 'absolute'. If 'raw', bootstrapping will be done on the raw prediction errors. If 'absolute', bootstrapping will be done on the absolute prediction errors with random signs. Default is 'raw'
#' @param n_bootstraps The number of bootstraps to perform. Default is 1000
#' @param distance_weighted_bootstrap Logical. If TRUE, the function will use distance-weighted bootstrapping. Default is FALSE. If TRUE, the probability of selecting a prediction error is weighted by the distance to the predicted value using the specified distance function and weight function. If FALSE, standard bootstrapping is performed.
#' @details
#' This function estimates prediction intervals using bootstrapped prediction errors derived from a calibration set. It supports both standard and distance-weighted bootstrapping. The calibration set must consist of predicted values and corresponding true values, either provided as separate vectors or as a two-column tibble or matrix. Alternatively, users may provide a vector of precomputed prediction errors if model predictions and truths are already processed.
#'
#' Two types of error can be used for bootstrapping:
#' - `"raw"`: bootstrapping is performed on the raw signed prediction errors (truth - prediction), allowing for asymmetric prediction intervals.
#' - `"absolute"`: bootstrapping is done on the absolute errors, and random signs are applied when constructing intervals. This results in (approximately) symmetric intervals around the prediction.
#'
#' Distance-weighted bootstrapping (`distance_weighted_bootstrap = TRUE`) can be used to give more weight to calibration errors closer to each test prediction. Distances are computed between the feature matrices or vectors supplied via `distance_features_calib` and `distance_features_pred`. These distances are then transformed into weights using the selected kernel in `weight_function`, with rapidly decaying kernels (e.g., Gaussian) emphasizing strong locality and slower decays (e.g., reciprocal or Cauchy) providing smoother influence. Distances can be geographic coordinates, predicted values, or any other relevant features that capture similarity in the context of the prediction task. The distance metric is specified via `distance_type`, with options for Mahalanobis or Euclidean distance. The default is Mahalanobis distance, which accounts for correlations between features. Normalization of distances can be applied using the `normalize_distance` parameter. Normalization is primarily useful for euclidean distances to ensure that features on different scales do not disproportionately influence the distance calculations.
#'
#' The number of bootstrap samples is controlled via the `n_bootstraps` parameter. For computational efficiency, this can be reduced at the cost of interval precision.
#'
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
	)
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

	return(boot_set)
}
