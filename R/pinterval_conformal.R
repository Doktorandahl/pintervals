#' Conformal Prediction Intervals of Continuous Values
#'
#' @description
#' This function calculates conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction. The intervals are computed using either a calibration set with predicted and true values or a set of pre-computed non-conformity scores from the calibration set. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals.
#'
#' @param pred Vector of predicted values
#' @param calib  A numeric vector of predicted values in the calibration partition, or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values. If calib is a numeric vector, calib_truth must be provided.
#' @param calib_truth A numeric vector of true values in the calibration partition. Only required if calib is a numeric vector
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param ncs_function The function to compute nonconformity scores. Default is 'absolute_error'. The user can also provide a custom function, or a string that matches a function, which computes the nonconformity scores. This function should take two arguments, a vector of predicted values and a vector of true values, in that order, and should return a numeric vector of nonconformity scores.
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#' @param calibrate = FALSE Logical. If TRUE, the function will calibrate the predictions and intervals using the calibration set. Default is FALSE. See details for more information on calibration.
#' @param calibration_method The method to use for calibration. Can be "glm" or "isotonic". Default is "glm". Only used if calibrate = TRUE.
#' @param calibration_family The family used for the calibration model. Default is "gaussian". Only used if calibrate = TRUE and calibration_method = "glm".
#' @param calibration_transform Optional transformation to apply to the predictions before calibration. Default is NULL. Only used if calibrate = TRUE and calibration_method = "glm".
#' @param resolution The minimum step size for the grid search. Default is 0.01. See details for more information.
#' @param grid_size Alternative to `resolution`, the number of points to use in the grid search between the lower and upper bound. If provided, resolution will be ignored.
#'
#' @details
#' This function computes prediction intervals using inductive conformal prediction. The calibration set must include predicted values and true values. These can be provided either as separate vectors (`calib`and `calib_truth`) or as a two-column tibble or matrix where the first column contains the predicted values and the second column contains the true values. If `calib` is a numeric vector, `calib_truth` must also be provided.
#'
#' Non-conformity scores are calculated using the specified `ncs_function`, which can be `"absolute_error"` or a user-defined custom function which computes the non-conformity scores. A custom function should take two arguments, a vector of predicted values and a vector of true values, in that order, and return a numeric vector of non-conformity scores.
#'
#' To determine the prediction intervals, the function performs a grid search over a specified range of possible outcome values, identifying intervals that satisfy the desired confidence level of \(1 - \eqn{\alpha}\). The user can define the range via the `lower_bound` and `upper_bound` parameters. If these are not supplied, the function defaults to using the minimum and maximum of the true values in the calibration data.
#'
#' The resolution of the grid search can be controlled by either the `resolution` argument, which sets the minimum step size, or the `grid_size` argument, which sets the number of grid points. For wide prediction spaces, the grid search may be computationally intensive. In such cases, increasing the `resolution` or reducing the `grid_size` may improve performance.
#'
#' Optionally, the predicted values can be calibrated before interval construction by setting `calibrate = TRUE`. In this case, the predictions are passed through `calibrate_predictions()` to adjust the predictions based on the calibration set. The calibration method can be specified using `calibration_method` and `calibration_family`, with "glm" being the default method. See \link[pintervals]{calibrate_predictions} for more information on calibration.
#'
#' The function returns a tibble with the predicted values and their corresponding lower and upper prediction interval bounds.
#'
#'
#' @return A tibble with the predicted values and the lower and upper bounds of the prediction intervals.
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
#'' # Calculate prediction intervals using conformal prediction. Setting lower_bound to 0 as this is the the minimum value of the true values.
#' pinterval_conformal(pred_test,
#' calib = calib_pred
#' calib_truth = calib_truth,
#' alpha = 0.1,
#' lower_bound = 0)
#'
pinterval_conformal <- function(pred,
													calib = NULL,
													alpha = 0.1,
													ncs_function = 'absolute_error',
													lower_bound = NULL,
													upper_bound = NULL,
													calibrate = FALSE,
													calibration_method = 'glm',
													calibration_family = 'gaussian',
													calibration_transform = NULL,
													calib_truth = NULL,
													resolution = 0.01,
													grid_size = NULL){

		if(is.numeric(calib) && is.null(calib_truth)){
			stop('If calib is numeric, calib_truth must be provided')
		}
		if(!is.numeric(calib) && ncol(calib)!=2){
			stop('calib must be a numeric vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values')
		}

	if(!is.numeric(alpha) || alpha<=0 || alpha>=1 || length(alpha)!=1){
		stop('alpha must be a single numeric value between 0 and 1')
	}
	if(is.character(ncs_function)){

	if(ncs_function == 'absolute_error'){
		ncs_function <- abs_error
	}else{
		ncs_function <- match.fun(ncs_function)
	}
		}else if(!is.function(ncs_function)){
		stop('ncs_function must be a function or a character string matching a function. Note, the ncs_function must take two arguments, a vector of predicted values and a vector of true values, in that order')
	}

	if(!is.numeric(pred)){
		stop('pred must be a numeric vector')
	}

	if(!is.numeric(resolution) || resolution<=0 || length(resolution)!=1){
		stop('resolution must be a single numeric value greater than 0')
	}



		if(!is.numeric(calib)){
			calib_org <- calib
			if(is.matrix(calib)){
				calib <- as.numeric(calib_org[,1])
				calib_truth <- as.numeric(calib_org[,2])
			}else{
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
			}
		}
		ncs <- ncs_function(calib,calib_truth)

	if(is.null(lower_bound)){
		lower_bound <- min(calib_truth)
	}
	if(is.null(upper_bound)){
		upper_bound <- max(calib_truth)
	}

	cp_set <- grid_finder(y_min = lower_bound,
												y_max = upper_bound,
												ncs = ncs,
												y_hat = pred,
												min_step=resolution,
												alpha=alpha,
												grid_size=grid_size,
												ncs_function = ncs_function,
												calib = calib)

	if(calibrate){
		calibrated_preds <- calibrate_predictions(pred = pred,
																							calib = calib,
																							calib_truth = calib_truth,
																							method = calibration_method,
																							family = calibration_family,
																							transform = calibration_transform)
		cp_set <- cp_set %>%
			dplyr::mutate(pred = calibrated_preds)
	}

	return(cp_set)
}
