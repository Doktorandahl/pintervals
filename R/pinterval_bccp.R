#' Bin-conditional conformal prediction intervals for continuous predictions
#'
#'@description
#'This function calculates bin-conditional conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction on a bin-by-bin basis. The intervals are computed using a calibration set with predicted and true values and their associated bins. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals. Bin-conditional conformal prediction intervals are useful when the prediction error is not constant across the range of predicted values and ensures that the coverage is (approximately) correct for each bin under the assumption that the non-conformity scores are exchangeable within each bin.
#'
#' @inheritParams pinterval_conformal
#' @param calib_bins A vector of bin identifiers for the calibration set. Not used if breaks are provided.
#' @param breaks A vector of break points for the bins to manually define the bins. If NULL, lower and upper bounds of the bins are calculated as the minimum and maximum values of each bin in the calibration set. Must be provided if calib_bins is not provided, either as a vector or as the last column of a calib tibble.
#' @param right Logical, if TRUE the bins are right-closed (a,b] and if FALSE the bins are left-closed `[ a,b)`. Only used if breaks are provided.
#' @param contiguize Logical indicating whether to contiguize the intervals. TRUE will consider all bins for each prediction using the lower and upper endpoints as interval limits to avoid non-contiguous intervals. FALSE will allows for non-contiguous intervals. TRUE guarantees at least appropriate coverage in each bin, but may suffer from over-coverage in certain bins. FALSE will have appropriate coverage in each bin but may have non-contiguous intervals. Default is FALSE.
#'
#' @details
#' `pinterval_bccp()` extends [pinterval_conformal()] to the
#' bin-conditional setting, where prediction intervals are
#' calibrated separately within user-specified bins. It is particularly useful when prediction error varies across the range of predicted values, as it enables locally valid coverage by ensuring that the coverage level \eqn{1 - \alpha} holds within each bin—assuming exchangeability of non-conformity scores within bins.
#'
#' For a detailed description of non-conformity scores, distance weighting and the general inductive conformal framework, see [pinterval_conformal()].
#'
#' For `pinterval_bccp()`, the calibration set must include predicted values, true values, and corresponding bin identifiers or breaks for the bins. These can be provided either as separate vectors (`calib`, `calib_truth`, and `calib_bins` or `breaks`).
#'
#' Bins endpoints can be defined manually via the `breaks` argument or inferred from the calibration data. If `contiguize = TRUE`, the function ensures the resulting prediction intervals are contiguous across bins, potentially increasing coverage beyond the nominal level in some bins. If `contiguize = FALSE`, the function may produce non-contiguous intervals, which are more efficient but may be harder to interpret.
#'
#' @seealso \code{\link[pintervals]{pinterval_conformal}}
#'
#' @return A tibble with the predicted values and the lower and upper bounds of the prediction intervals. If contiguize = FALSE, the intervals may consist of multiple disjoint segments; in this case, the tibble will contain a list-column with all segments for each prediction.
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
#'
#' # Create bins based on quantiles
#' 	bin <- cut(y, breaks = quantile(y, probs = seq(0, 1, 1/4)),
#' 	include.lowest = TRUE, labels =FALSE)
#' df <- tibble(x1, x2, y, bin)
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model to the training data
#' mod <- lm(log(y) ~ x1 + x2, data=df_train)
#'
#' # Generate predictions on the original y scale for the calibration data
#' calib <- exp(predict(mod, newdata=df_cal))
#' calib_truth <- df_cal$y
#' calib_bins <- df_cal$bin
#'
#' # Generate predictions for the test data
#'
#' pred_test <- exp(predict(mod, newdata=df_test))
#'
#' # Calculate bin-conditional conformal prediction intervals
#' pinterval_bccp(pred = pred_test,
#' calib = calib,
#' calib_truth = calib_truth,
#' calib_bins = calib_bins,
#' alpha = 0.1)
#'
pinterval_bccp = function(
	pred,
	calib = NULL,
	calib_truth = NULL,
	calib_bins = NULL,
	breaks = NULL,
	right = TRUE,
	contiguize = FALSE,
	alpha = 0.1,
	ncs_type = c(
		'absolute_error',
		'relative_error',
		'za_relative_error',
		'heterogeneous_error',
		'raw_error'
	),
	grid_size = 10000,
	resolution = NULL,
	distance_weighted_cp = FALSE,
	distance_features_calib = NULL,
	distance_features_pred = NULL,
	normalize_distance = TRUE,
	distance_type = c('mahalanobis', 'euclidean'),
	weight_function = c(
		'gaussian_kernel',
		'caucy_kernel',
		'logistic',
		'reciprocal_linear'
	)
) {
	i <- NA

	if (!is.numeric(pred)) {
		stop('pred must be a numeric scalar or vector')
	}

	if (is.numeric(calib) & is.null(calib_truth)) {
		stop('If calib is numeric, calib_truth must be provided')
	}
	if (!is.numeric(calib) && ncol(calib) < 2) {
		stop(
			'calib must be a numeric vector or a 2 or 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and (optionally) the third column being the bin values if bin structure is not provided in argument bins'
		)
	}

	if (
		(is.null(breaks)) &&
			(is.null(calib_bins) || (!is.numeric(calib) && ncol(calib) != 3))
	) {
		stop(
			'If breaks are not provided, bins for the calibration set must be provided as a vector or a as the last column of the calib if calib is a tibble or matrix'
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

	if ((is.null(breaks))) {
		warning(
			'No explicit bin structure provided, breaks are calculated as the minimum and maximum values for each bin in the calibration set'
		)
	}

	if (!is.null(breaks) & !is.null(calib_bins)) {
		warning("If breaks are provided, calib_bins will be ignored")
	}

	distance_type <- match.arg(distance_type, c('mahalanobis', 'euclidean'))

	if (!is.numeric(calib)) {
		calib_org <- calib
		if (is.matrix(calib)) {
			calib <- as.numeric(calib_org[, 1])
			calib_truth <- as.numeric(calib_org[, 2])
			if (is.null(calib_bins) && ncol(calib_org) == 3 && !is.null(breaks)) {
				calib_bins <- as.numeric(calib_org[, 3])
			}
		} else {
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
			if (is.null(calib_bins) && ncol(calib_org) == 3 && !is.null(breaks)) {
				calib_bins <- as.numeric(calib_org[[3]])
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

	if (is.null(calib_bins)) {
		calib_bins <- cut(
			calib_truth,
			breaks = breaks,
			labels = FALSE,
			right = right
		)
	}

	nobs_bins <- as.numeric(table(calib_bins))

	if (any(nobs_bins * alpha / 2 < 1)) {
		warning(
			'Some bins have too few observations to calculate prediction intervals at the specified alpha level. Consider using a larger calibration set or a higher alpha level'
		)
	}

	bin_labels <- sort(unique(calib_bins))
	if (length(bin_labels) < 2) {
		stop(
			'Calibration set must have at least two bins. For continuous prediction intervals without bins, use pinterval_conformal() instead of pinterval_bccp()'
		)
	}

	lower_bounds <- foreach::foreach(i = bin_labels, .final = unlist) %do%
		min(calib_truth[calib_bins == i])
	upper_bounds <- foreach::foreach(i = bin_labels, .final = unlist) %do%
		max(calib_truth[calib_bins == i])

	cp_intervals <- foreach::foreach(i = 1:length(bin_labels)) %do%
		suppressWarnings(pinterval_conformal(
			pred = pred,
			lower_bound = lower_bounds[i],
			upper_bound = upper_bounds[i],
			ncs_type = ncs_type,
			calib = calib[calib_bins == bin_labels[i]],
			calib_truth = calib_truth[calib_bins == bin_labels[i]],
			distance_weighted_cp = distance_weighted_cp,
			distance_features_calib = distance_features_calib[
				calib_bins == bin_labels[i],
			],
			distance_features_pred = distance_features_pred,
			distance_type = distance_type,
			normalize_distance = normalize_distance,
			weight_function = weight_function,
			alpha = alpha,
			resolution = resolution,
			grid_size = grid_size
		))

	cp_intervals2 <- flatten_cp_bin_intervals(
		cp_intervals,
		contiguize = contiguize
	)

	return(cp_intervals2)
}
