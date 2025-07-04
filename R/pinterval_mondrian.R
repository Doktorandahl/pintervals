#' Mondrian conformal prediction intervals for continuous predictions
#'
#'@description
#'This function calculates Mondrian conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction on a Mondrian class-by-class basis. The intervals are computed using a calibration set with predicted and true values and their associated classes. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals. Mondrian conformal prediction intervals are useful when the prediction error is not constant across groups or classes, as they allow for locally valid coverage by ensuring that the coverage level \(1 - \alpha\) holds within each class—assuming exchangeability of non-conformity scores within classes.
#'
#' @param pred Vector of predicted values or a 2 column tibble or matrix with the first column being the predicted values and the second column being the Mondrian class labels. If pred is a numeric vector, pred_class must be provided.
#' @param pred_class A vector of class identifiers for the predicted values. This is used to group the predictions by class for Mondrian conformal prediction.
#' @param calib A numeric vector of predicted values in the calibration partition or a 3 column tibble or matrix with the first column being the predicted values and the second column being the truth values and the third column being the Mondrian class labels. If calib is a numeric vector, calib_truth and calib_class must be provided.
#' @param calib_truth A numeric vector of true values in the calibration partition
#' @param calib_class A vector of class identifiers for the calibration set.
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param ncs_function The function to compute nonconformity scores. Default is 'absolute_error'. The user can also provide a custom function, or a string that matches a function, which computes the nonconformity scores. This function should take two arguments, a vector of predicted values and a vector of true values, in that order, and should return a numeric vector of nonconformity scores.
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#' @param calibrate = FALSE Logical. If TRUE, the function will calibrate the predictions and intervals using the calibration set. Default is FALSE.
#' @param calibration_method The method to use for calibration. Can be "glm" or "isotonic". Default is "glm". Only used if calibrate = TRUE.
#' @param calibration_family The family used for the calibration model. Default is "gaussian". Only used if calibrate = TRUE and calibration_method = "glm". See `calibrate_predictions()` for more information.
#' @param calibration_transform Optional transformation to apply to the predictions before calibration. Default is NULL. Only used if calibrate = TRUE and calibration_method = "glm". See `calibrate_predictions()` for more information.
#' @param resolution The minimum step size for the grid search. Default is 0.01. See details for more information.
#' @param grid_size Alternative to `resolution`, the number of points to use in the grid search between the lower and upper bound. If provided, resolution will be ignored.
#'
#' @return A tibble with the predicted values, the lower and upper bounds of the prediction intervals. If treat_noncontiguous is 'non_contiguous', the lower and upper bounds are set in a list variable called 'intervals' where all non-contiguous intervals are stored.
#'
#' @details
#' This function computes Mondrian conformal prediction intervals using inductive conformal prediction applied separately within each class (also called strata or groups) of the calibration data. It is especially useful when prediction error varies systematically across known categories, allowing for class-conditional validity by ensuring that the prediction intervals attain the desired coverage level \(1 - \alpha\) within each class—under the assumption of exchangeability within classes.
#'
#' The calibration set must include predicted values, true values, and corresponding class labels. These can be supplied as separate vectors (`calib`, `calib_truth`, and `calib_class`) or as a single three-column matrix or tibble.
#'
#' Non-conformity scores are calculated using the specified `ncs_function`, which can be `"absolute_error"` or a user-defined custom function which computes the non-conformity scores. A custom function should take two arguments, a vector of predicted values and a vector of true values, in that order, and return a numeric vector of non-conformity scores.
#'
#' To determine the prediction intervals, the function performs a grid search over a specified range of possible outcome values, identifying intervals that satisfy the desired confidence level of \(1 - \alpha\). The user can define the range via the `lower_bound` and `upper_bound` parameters. If these are not supplied, the function defaults to using the minimum and maximum of the true values in the calibration data.
#'
#' The resolution of the grid search can be controlled by either the `resolution` argument, which sets the minimum step size, or the `grid_size` argument, which sets the number of grid points. For wide prediction spaces, the grid search may be computationally intensive. In such cases, increasing the `resolution` or reducing the `grid_size` may improve performance.
#'
#' Optionally, the predicted values can be calibrated before interval construction by setting `calibrate = TRUE`. In this case, the predictions are passed through `calibrate_predictions()` to adjust the predictions based on the calibration set. The calibration method can be specified using `calibration_method` and `calibration_family`, with "glm" being the default method.
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
#'
#' # Apply Mondrian conformal prediction
#' pinterval_mondrian(pred = pred_test,
#'                    calib = calib,
#'                    calib_truth = calib_truth,
#'                    calib_class = calib_class,
#'                    alpha = 0.1)
pinterval_bccp = function(pred,
													pred_class = NULL,
													calib = NULL,
													calib_truth = NULL,
													calib_class = NULL,
													lower_bound = NULL,
													upper_bound = NULL,
													alpha = 0.1,
													ncs_function = 'absolute_error',
													calibrate = FALSE,
													calibration_method = 'glm',
													calibration_family = NULL,
													resolution = 0.01,
													grid_size = NULL){

	i <- NA


	if(!is.numeric(pred) && ncol(pred) != 2){
		stop('pred must be a numeric scalar or vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the class labels')
	}

	if(is.numeric(pred) && is.null(pred_class)){
		stop('If pred is numeric, pred_class must be provided')
	}

	if(!is.numeric(pred)){
		pred_class <- as.numeric(pred[[2]])
		pred <- as.numeric(pred[[1]])
	}

	if(is.numeric(calib) & is.null(calib_truth)){
		stop('If calib is numeric, calib_truth must be provided')
	}

	if(is.numeric(calib) & is.null(calib_class)){
		stop('If calib is numeric, calib_class must be provided')
	}
	if(!is.numeric(calib) && ncol(calib)<3){
		stop('calib must be a numeric vector or a 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and the third column being the class labels')
	}

	if(!is.numeric(alpha) || alpha<=0 || alpha>=1 || length(alpha)!=1){
		stop('alpha must be a single numeric value between 0 and 1')
	}

	if(is.character(ncs_function)){
		ncs_function <- match.arg(ncs_function, c('absolute_error'))
	}

	if(ncs_function == 'absolute_error'){
		ncs_function <- abs_error
	}else if(is.character(ncs_function)){
		ncs_function <- match.fun(ncs_function)
	}else if(!is.function(ncs_function) & is.null(ncs)){
		stop('ncs_function must be a function or a character string matching a function if ncs is not provided. The ncs_function must take two arguments, a vector of predicted values and a vector of true values, in that order')
	}

	if(!is.numeric(calib)){
		calib_org <- calib
		if(is.matrix(calib)){
			calib <- as.numeric(calib_org[,1])
			calib_truth <- as.numeric(calib_org[,2])
			if(is.null(calib_class) && ncol(calib_org) == 3 && !is.null(breaks)){
				calib_class <- as.numeric(calib_org[,3])
			}
		}else{
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
			if(is.null(calib_class) && ncol(calib_org) == 3 && !is.null(breaks)){
				calib_class <- as.numeric(calib_org[[3]])
			}
		}
	}


	nobs_class <- as.numeric(table(calib_class))

	if(any(nobs_class*alpha/2<1)){
		warning('Some classes have too few observations to calculate prediction intervals at the specified alpha level. Consider using a larger calibration set or a higher alpha level')
	}

	class_labels <- sort(unique(calib_class))
	if(length(calib_class)<2){
		stop('Calibration set must have at least two classes For continuous prediction intervals without classes, use pinterval_conformal() instead of pinterval_mondrian()')
	}

	cp_intervals <- foreach::foreach(i = 1:length(class_labels),.final=bind_rows) %do%{
		indices <- which(calib_class == class_labels[i])
		res <- suppressWarnings(pinterval_conformal(pred = pred[pred_class==class_labels[i]],
																				 lower_bound = lower_bound,
																				 upper_bound = upper_bound,
																				 ncs_function = ncs_function,
																				 calib = calib[calib_class==class_labels[i]],
																				 calib_truth = calib_truth[calib_class==class_labels[i]],
																				 calibrate = calibrate,
																				 calibration_method = calibration_method,
																				 calibration_family = calibration_family,
																				 alpha = alpha,
																				 resolution = resolution,
																				 grid_size = grid_size))
		res$indices <- indices
		res
	}

cp_intervals2 <- cp_intervals %>%
	arrange(indices) %>%
	select(-indices) %>%
	mutate(class = class_labels)

	return(cp_intervals2)
}
