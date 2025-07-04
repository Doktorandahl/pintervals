#' Bootstrap prediction intervals
#'
#' @description
#' This function computes bootstrapped prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using bootstrapped prediction errors. The prediction errors to bootstrap from are computed using either a calibration set with predicted and true values or a set of pre-computed prediction errors from a calibration dataset or other data which the model was not trained on (e.g. OOB errors from a model using bagging). The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals.
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values
#' @param calib_truth A numeric vector of true values in the calibration partition. Only required if calib is a numeric vector
#' @param error An optional numeric vector of pre-computed prediction errors from a calibration partition or other test data. If provided, calib will be ignored
#' @param error_type The type of error to use for the prediction intervals. Can be 'raw' or 'absolute'. If 'raw', bootstrapping will be done on the raw prediction errors. If 'absolute', bootstrapping will be done on the absolute prediction errors with random signs. Default is 'raw'
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param n_bootstraps The number of bootstraps to perform. Default is 1000
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#'
#' @return A tibble with the predicted values, lower bounds, and upper bounds of the prediction intervals
#' @export
#'
#' @examples
#'
#' library(dplyr)
#' library(tibble)
#' x1 <- runif(1000)
#' x2 <- runif(1000)
#' y <- rlnorm(1000, meanlog = x1 + x2, sdlog = 0.5)
#' df <- tibble(x1, x2, y)
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#' mod <- lm(log(y) ~ x1 + x2, data=df_train)
#' calib <- exp(predict(mod, newdata=df_cal))
#' calib_truth <- df_cal$y
#' pred_test <- exp(predict(mod, newdata=df_test))
#'
#' pinterval_bootstrap(pred = pred_test,
#' calib = calib,
#' calib_truth = calib_truth,
#' error_type = 'raw',
#' alpha = 0.1,
#' lower_bound = 0)
pinterval_bootstrap <- function(pred,
															 calib = NULL,
															 calib_truth = NULL,
															 error = NULL,
															 error_type = c('raw','absolute'),
															 alpha = 0.1,
															 n_bootstraps=1000,
															 lower_bound = NULL,
															 upper_bound = NULL){

	i <- NA
	if(!is.numeric(pred)){
		stop('pred must be a single number or a numeric vector')
	}

	if(is.numeric(calib) & (is.null(calib_truth))){
		stop('If calib is numeric, calib_truth must be provided')
	}

	if(is.null(error)){
		if(is.null(calib)){
			stop('Either calib or error must be provided')
		}

	if(!is.numeric(calib) && ncol(calib)!=2){
		stop('calib must be a numeric vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values, unless error is provided')
	}
	}

	if(!is.null(error) & !is.null(calib)){
		warning("Both error and calib provided, error will be used")
	}

	if(!is.numeric(calib) & !is.null(error)){
		calib_org <- calib
		if(is.matrix(calib)){
			calib <- as.numeric(calib_org[,1])
			calib_truth <- as.numeric(calib_org[,2])
		}else{
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
		}
	}

if(is.null(lower_bound)){
		lower_bound <- -Inf
	}
	if(is.null(upper_bound)){
		upper_bound <- Inf
	}

	error_type <- match.arg(error_type, c('raw','absolute'))

	if(is.null(error)){
		if(error_type == 'raw'){
			error <- calib - calib_truth
		}
		else if(error_type == 'absolute'){
			error <- abs(calib - calib_truth)
			error <- c(error, -error)
		}
	}

	boot_set <- foreach::foreach(i = 1:length(pred)) %do%
		bootstrap_inner(pred = pred[i], error = error, nboot = n_bootstraps,
										alpha = alpha, lower_bound = lower_bound, upper_bound = upper_bound)

	return(dplyr::bind_rows(boot_set))
}
