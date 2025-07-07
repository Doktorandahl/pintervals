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
#' @param dw_bootstrap Logical. If TRUE, the function will use distance-weighted bootstrapping. Default is FALSE. If TRUE, the probability of selecting a prediction error is weighted by the distance to the predicted value. The probability is computed as 1/(1 + distance), where distance is the distance between the predicted value and the calibration values.
#' @param distance_function A function to compute the distance between the predicted value and the calibration values. Default is 'absolute' which computes the absolute distance between the predicted value and the calibration values. Can also be a function, or a character string matching a function, that takes two arguments, a single predicted value and a vector of calibration values, and returns a numeric vector of distances. See details for more information.
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
#' Distance-weighted bootstrapping (`dw_bootstrap = TRUE`) can be used to give more weight to calibration errors closer to each test prediction. Distances are calculated using the specified `distance_function`, which defaults to the absolute difference. Distance weighted bootstrapping assigns a probability of selection for each calibration error based on its distance to the predicted value, computed as \(1/(1 + distance)\). The distance function can be a custom function or a character string that matches a function, which takes two arguments: a single predicted value and a vector of calibration values, returning a numeric vector of distances.
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
pinterval_bootstrap <- function(pred,
															 calib,
															 calib_truth = NULL,
															 error_type = c('raw','absolute'),
															 alpha = 0.1,
															 n_bootstraps=1000,
																dw_bootstrap = FALSE,
																distance_function = 'absolute',
															 calibrate = FALSE,
															 calibration_method = c('glm', 'isotonic'),
															 calibration_family = 'gaussian',
															 calibration_transform = NULL){

	i <- NA
	if(!is.numeric(pred)){
		stop('pred must be a single number or a numeric vector')
	}

	if(is.numeric(calib) & (is.null(calib_truth))){
		stop('If calib is numeric, calib_truth must be provided')
	}

	if(!is.numeric(calib) && ncol(calib)!=2){
		stop('calib must be a numeric vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values')
	}

	if(calibrate && dw_bootstrap){
		warning('calibrate is TRUE, but dw_bootstrap is also TRUE. Calibration will be performed after bootstrapping, and the corresponding prediction intervals will be adjusted accordingly.')
	}

	if(is.character(distance_function)){
		if(distance_function == 'absolute'){
			distance_function <- function(x, y) abs(x - y)}
		else{
			distance_function <- match.fun(distance_function)
		}
	}else if(!is.function(distance_function)){
		stop('distance_function must be a character string or a function')
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

	error_type <- match.arg(error_type, c('raw','absolute'))

		if(error_type == 'raw'){
			error <- calib - calib_truth
		}else if(error_type == 'absolute'){
			error <- abs(calib - calib_truth)
		}

	boot_set <- foreach::foreach(i = 1:length(pred)) %do%
		bootstrap_inner(pred = pred[i], calib = calib, error = error, nboot = n_bootstraps,
										alpha = alpha,error_type = error_type,
										dw_bootstrap = dw_bootstrap,
										distance_function = distance_function)

	boot_set <- dplyr::bind_rows(boot_set)

	if(calibrate){
		calibration_method <- match.arg(calibration_method, c('glm', 'isotonic'))
		if(calibration_method == 'glm'){
			calibrated_predictions <- calibrate_predictions(pred, calib, calib_truth, method = calibration_method,
																										family = calibration_family,
																										transform = calibration_transform)
		}else if(calibration_method == 'isotonic'){
			calibrated_predictions <- calibrate_predictions(pred, calib, calib_truth, method = calibration_method)
		}

		calibrated_diffs <- pred - calibrated_predictions
		boot_set <- boot_set %>%
			dplyr::mutate(lower_bound = rlang::.data$lower_bound + calibrated_diffs,
									 upper_bound = rlang::.data$upper_bound + calibrated_diffs,
																			 pred = calibrated_predictions)

	}

	return(boot_set)
}
