#' Continuous Conformal Prediction Intervals
#'
#' @description
#' This function calculates conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction. The intervals are computed using either a calibration set with predicted and true values or a set of pre-computed non-conformity scores from the calibration set. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals.
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values
#' @param calib_truth A numeric vector of true values in the calibration partition. Only required if calib is a numeric vector
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param ncs_function A function or a character string matching a function that takes two arguments, a vector of predicted values and a vector of true values, in that order. The function should return a numeric vector of nonconformity scores. Default is 'absolute_error' which returns the absolute difference between the predicted and true values.
#' @param ncs A numeric vector of pre-computed nonconformity scores from a calibration partition. If provided, calib will be ignored
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#' @param min_step The minimum step size for the grid search. Default is 0.01. Useful to change if predictions are made on a discrete grid or if the resolution of the interval is too coarse or too fine.
#' @param grid_size Alternative to min_step, the number of points to use in the grid search between the lower and upper bound. If provided, min_step will be ignored.
#'
#' @return A tibble with the predicted values and the lower and upper bounds of the prediction intervals.
#' @export
#'
#' @examples
#'
pinterval_cp_cont <- function(pred,
													calib = NULL,
													calib_truth = NULL,
													alpha = 0.1,
													ncs_function = 'absolute_error',
													ncs = NULL,
													lower_bound = NULL,
													upper_bound = NULL,
													min_step = 0.01,
													grid_size = NULL){

	if(is.null(calib) & is.null(ncs)){
		stop('Either calib or ncs must be provided')
	}
	if(is.null(ncs)){
		if(is.numeric(calib) & is.null(calib_truth)){
			stop('If calib is numeric, calib_truth must be provided')
		}
		if(!is.numeric(calib) && ncol(calib)!=2){
			stop('calib must be a numeric vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values')
		}
	}else{
		if(!is.numeric(ncs)){
			stop('ncs must be a numeric vector')
		}
		if(!is.null(calib)){
			warning('ncs provided, calib will be ignored')
		}
	}
	if(!is.numeric(alpha) || alpha<=0 || alpha>=1 || length(alpha)!=1){
		stop('alpha must be a single numeric value between 0 and 1')
	}

	if(ncs_function == 'absolute_error'){
		ncs_function <- abs_error
	}else if(is.character(ncs_function)){
		ncs_function <- match.fun(ncs_function)
	}else if(!is.function(ncs_function)){
		stop('ncs_function must be a function or a character string matching a function. Note, the ncs_function must take two arguments, a vector of predicted values and a vector of true values, in that order')
	}

	if(!is.numeric(pred)){
		stop('pred must be a numeric vector')
	}

	if(!is.numeric(min_step) || min_step<=0 || length(min_step)!=1){
		stop('min_step must be a single numeric value greater than 0')
	}

	# if(is.null(trans)){
	# 	trans <- dummy_fun
	# }else{
	# 	trans <- match.fun(trans)
	# }

	if(is.null(ncs)){
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
	}

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
												min_step=min_step,
												alpha=alpha,
												grid_size=grid_size,
												ncs_function = ncs_function)

	return(cp_set)
}
