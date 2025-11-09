#' Bin-conditional conformal prediction intervals for continuous predictions
#'
#'@description
#'This function calculates bin-conditional conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction on a bin-by-bin basis. The intervals are computed using a calibration set with predicted and true values and their associated bins. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals. Bin-conditional conformal prediction intervals are useful when the prediction error is not constant across the range of predicted values and ensures that the coverage is (approximately) correct for each bin under the assumption that the non-conformity scores are exchangeable within each bin.
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 or 3 column tibble or matrix with the first column being the predicted values and the second column being the truth values and the third, optional, column being the bin labels. If calib is a numeric vector, calib_truth and either calib_bins or breaks must be provided.
#' @param calib_truth A numeric vector of true values in the calibration partition
#' @param calib_bins A vector of bin identifiers for the calibration set. Not used if breaks are provided.
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
#' @param breaks A vector of break points for the bins to manually define the bins. If NULL, lower and upper bounds of the bins are calculated as the minimum and maximum values of each bin in the calibration set. Must be provided if calib_bins is not provided, either as a vector or as the last column of a calib tibble.
#' @param right Logical, if TRUE the bins are right-closed (a,b] and if FALSE the bins are left-closed `[ a,b)`. Only used if breaks are provided.
#' @param contiguize logical indicating whether to contiguize the intervals. TRUE will consider all bins for each prediction using the lower and upper endpoints as interval limits to avoid non-contiguous intervals. FALSE will allows for non-contiguous intervals. TRUE guarantees at least appropriate coverage in each bin, but may suffer from over-coverage in certain bins. FALSE will have appropriate coverage in each bin but may have non-contiguous intervals. Default is FALSE.
#' @param distance_weighted_cp Logical. If \code{TRUE}, weighted conformal prediction is performed where the non-conformity scores are weighted based on the distance between calibration and prediction points in feature space. Default is \code{FALSE}.
#'
#' @param distance_features_calib A matrix, data frame, or numeric vector of features from which to compute distances when \code{distance_weighted_cp = TRUE}. This should contain the feature values for the calibration set. Must have the same number of rows as the calibration set. Can be the predicted values themselves, or any other features which give a meaningful distance measure.
#' @param distance_features_pred A matrix, data frame, or numeric vector of feature values for the prediction set. Must be the same features as specified in \code{distance_features_calib}. Required if \code{distance_weighted_cp = TRUE}.
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
#' @param calibrate = FALSE Logical. If TRUE, the function will calibrate the predictions and intervals using the calibration set. Default is FALSE. See details for more information.
#' @param calibration_method The method to use for calibration. Can be "glm" or "isotonic". Default is "glm". Only used if calibrate = TRUE.
#' @param calibration_family The family used for the calibration model. Default is "gaussian". Only used if calibrate = TRUE and calibration_method = "glm".
#' @param calibration_transform Optional transformation to apply to the predictions before calibration. Default is NULL. Only used if calibrate = TRUE and calibration_method = "glm".
#' @param resolution The minimum step size for the grid search. Default is 0.01. See details for more information.
#' @param grid_size Alternative to `resolution`, the number of points to use in the grid search between the lower and upper bound. If provided, resolution will be ignored.
#'
#' @return A tibble with the predicted values, the lower and upper bounds of the prediction intervals. If treat_noncontiguous is 'non_contiguous', the lower and upper bounds are set in a list variable called 'intervals' where all non-contiguous intervals are stored.
#'
#' @details
#' This function computes bin-conditional conformal prediction intervals using inductive conformal prediction applied separately within each bin of the calibration data. It is particularly useful when prediction error varies across the range of predicted values, as it enables locally valid coverage by ensuring that the coverage level \eqn{1 - \eqn{\alpha}} holds within each bin—assuming exchangeability of non-conformity scores within bins.
#'
#' The calibration set must include predicted values, true values, and corresponding bin identifiers or breaks for the bins. These can be provided either as separate vectors (`calib`, `calib_truth`, and `calib_bins` or `breaks`).
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
#' Bins endpoints can be defined manually via the `breaks` argument or inferred from the calibration data. If `contiguize = TRUE`, the function ensures the resulting prediction intervals are contiguous across bins, potentially increasing coverage beyond the nominal level in some bins. If `contiguize = FALSE`, the function may produce non-contiguous intervals, which are more efficient but may be harder to interpret.
#'
#' The prediction intervals are constructed using a grid search over a user-defined range of outcome values. The resolution of the grid search can be controlled by either the `resolution` argument, which sets the minimum step size, or the `grid_size` argument, which sets the number of grid points. For wide prediction spaces, the grid search may be computationally intensive. In such cases, increasing the `resolution` or reducing the `grid_size` may improve performance.
#'
#' Optionally, the predicted values can be calibrated before interval construction by setting `calibrate = TRUE`. In this case, the predictions are passed through `calibrate_predictions()` to adjust the predictions based on the calibration set. The calibration method can be specified using `calibration_method` and `calibration_family`, with "glm" being the default method. See \link[pintervals]{calibrate_predictions} for more information on calibration.
#'
#' The function returns a tibble containing the predicted values and their corresponding lower and upper bounds. If `contiguize = FALSE` and the resulting intervals are non-contiguous, the lower and upper bounds are in list-columns.
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
pinterval_bccp = function(pred,
									 calib = NULL,
									 calib_truth = NULL,
									 calib_bins = NULL,
									 breaks = NULL,
									 alpha = 0.1,
									 ncs_type = c('absolute_error',
									 						 'relative_error',
									 						 'za_relative_error',
									 						 'heterogeneous_error',
									 						 'raw_error'),
									 distance_weighted_cp = FALSE,
									 distance_features_calib = NULL,
									 distance_features_pred = NULL,
									 normalize_distance = TRUE,
									 weight_function = c('gaussian_kernel', 'caucy_kernel','logistic','reciprocal_linear'),
									 calibrate = FALSE,
									 calibration_method = 'glm',
									 calibration_family = NULL,
									 calibration_transform = NULL,
									 resolution = 0.01,
									 grid_size = NULL,
									 right = TRUE,
									 contiguize = FALSE){

	i <- NA


	if(!is.numeric(pred)){
		stop('pred must be a numeric scalar or vector')
	}

		if(is.numeric(calib) & is.null(calib_truth)){
			stop('If calib is numeric, calib_truth must be provided')
		}
		if(!is.numeric(calib) && ncol(calib)<2){
			stop('calib must be a numeric vector or a 2 or 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and (optionally) the third column being the bin values if bin structure is not provided in argument bins')
		}

		if((is.null(breaks)) && (is.null(calib_bins) || (!is.numeric(calib) && ncol(calib)!=3))){
			stop('If breaks are not provided, bins for the calibration set must be provided as a vector or a as the last column of the calib if calib is a tibble or matrix')
		}

	if(!is.numeric(alpha) || alpha<=0 || alpha>=1 || length(alpha)!=1){
		stop('alpha must be a single numeric value between 0 and 1')
	}

	ncs_type <- match.arg(ncs_type, c('absolute_error',
																		'relative_error',
																		'za_relative_error',
																		'heterogeneous_error',
																		'raw_error'))


	if((is.null(breaks))){
		warning('No explicit bin structure provided, breaks are calculated as the minimum and maximum values for each bin in the calibration set')
	}

	if(!is.null(breaks) & !is.null(calib_bins)){
		warning("If breaks are provided, calib_bins will be ignored")
	}

		if(!is.numeric(calib)){
			calib_org <- calib
			if(is.matrix(calib)){
				calib <- as.numeric(calib_org[,1])
				calib_truth <- as.numeric(calib_org[,2])
				if(is.null(calib_bins) && ncol(calib_org) == 3 && !is.null(breaks)){
					calib_bins <- as.numeric(calib_org[,3])
				}
			}else{
				calib_truth <- as.numeric(calib_org[[2]])
				calib <- as.numeric(calib_org[[1]])
				if(is.null(calib_bins) && ncol(calib_org) == 3 && !is.null(breaks)){
					calib_bins <- as.numeric(calib_org[[3]])
				}
			}
		}

	if(ncs_type == 'heterogeneous_error'){
		coefs <- stats::coef(stats::lm(abs(calib - calib_truth) ~ calib))
	}else{
		coefs <- NULL
	}

	if(distance_weighted_cp){
		if(is.null(distance_features_calib) || is.null(distance_features_pred)){
			stop('If distance_weighted_cp is TRUE, distance_features_calib and distance_features_pred must be provided')
		}
		if(!is.matrix(distance_features_calib) && !is.data.frame(distance_features_calib) && !is.numeric(distance_features_calib)){
			stop('distance_features_calib must be a matrix, data frame, or numeric vector')
		}
		if(!is.matrix(distance_features_pred) && !is.data.frame(distance_features_pred) && !is.numeric(distance_features_pred)){
			stop('distance_features_pred must be a matrix, data frame, or numeric vector')
		}
		if(is.numeric(distance_features_calib) && is.numeric(distance_features_pred)){
			if(length(distance_features_calib) != length(calib) || length(distance_features_pred) != length(pred)){
				stop('If distance_features_calib and distance_features_pred are numeric vectors, they must have the same length as calib and pred, respectively')
			}
		}else if(is.matrix(distance_features_calib) || is.data.frame(distance_features_calib)){
			if(nrow(distance_features_calib) != length(calib)){
				stop('If distance_features_calib is a matrix or data frame, it must have the same number of rows as calib')
			}
			if(ncol(distance_features_calib) != ncol(distance_features_pred)){
				stop('distance_features_calib and distance_features_pred must have the same number of columns')
			}
			if(nrow(distance_features_pred) != length(pred)){
				stop('If distance_features_pred is a matrix or data frame, it must have the same number of rows as pred')
			}

		}

		distance_features_calib <- as.matrix(distance_features_calib)
		distance_features_pred <- as.matrix(distance_features_pred)

		if(!is.function(weight_function)){
			weight_function <- match.arg(weight_function, c('gaussian_kernel', 'caucy_kernel','logistic','reciprocal_linear'))
			weight_function <- switch(weight_function,
																'gaussian_kernel' = function(d) exp(-d^2),
																'caucy_kernel' = function(d) 1 / (1 + d^2),
																'logistic' = function(d) 1 / (1 + exp(d)),
																'reciprocal_linear' = function(d) 1 / (1 + d))
		}
	}



	if(is.null(calib_bins)){
			calib_bins <- cut(calib_truth,breaks = breaks,labels = FALSE,right = right)
	}

	nobs_bins <- as.numeric(table(calib_bins))

	if(any(nobs_bins*alpha/2<1)){
		warning('Some bins have too few observations to calculate prediction intervals at the specified alpha level. Consider using a larger calibration set or a higher alpha level')
	}

		bin_labels <- sort(unique(calib_bins))
		if(length(bin_labels)<2){
			stop('Calibration set must have at least two bins. For continuous prediction intervals without bins, use pinterval_conformal() instead of pinterval_bccp()')
		}

		lower_bounds <- foreach::foreach(i = bin_labels,.final = unlist) %do% min(calib_truth[calib_bins==i])
		upper_bounds <- foreach::foreach(i = bin_labels,.final = unlist) %do% max(calib_truth[calib_bins==i])

	cp_intervals <- foreach::foreach(i = 1:length(bin_labels)) %do%
		suppressWarnings(pinterval_conformal(pred = pred,
																			 lower_bound = lower_bounds[i],
																			 upper_bound = upper_bounds[i],
																			 ncs_type = ncs_type,
																			 calib = calib[calib_bins==bin_labels[i]],
																			 calib_truth = calib_truth[calib_bins==bin_labels[i]],
																			 calibrate = calibrate,
																			 calibration_method = calibration_method,
																			 calibration_family = calibration_family,
																			 alpha = alpha,
																			 resolution = resolution,
																			 grid_size = grid_size))




	cp_intervals2 <- flatten_cp_bin_intervals(cp_intervals, contiguize = contiguize)

	if(calibrate){
		calibrated_preds <- calibrate_predictions(pred = pred,
																							calib = calib,
																							calib_truth = calib_truth,
																							method = calibration_method,
																							family = calibration_family,
																							transform = calibration_transform)
		cp_intervals2 <- cp_intervals2 %>%
			dplyr::mutate(pred = calibrated_preds)
	}

	return(cp_intervals2)
}
