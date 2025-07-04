

#' Parametric prediction intervals for continuous predictions
#'
#' @description
#' This function computes parametric prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using a user specified parametric distribution and parameters. The distribution can be any distribution available in R or a user defined distribution as long as a quantile function is available. The parameters should be estimated on calibration data. The prediction intervals are calculated as the quantiles of the distribution at the specified confidence level.
#'
#'
#' @param pred Vector of predicted values
#' @param dist Distribution to use for the prediction intervals. Can be a character string matching any available distribution in R or a function representing a distribution. If a function is provided, it must be a quantile function (e.g. qnorm, qgamma, etc.)
#' @param pars List of named parameters for the distribution for each prediction. See details for more information.
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#'
#' @details
#' The distributions are not limited to the standard distributions available in R. Any distribution can be used as long as a quantile function is available. Users may create their own distribution functions and plug in the resulting quantile function or create compositie or mixture distributions using for instance the package `mistr` and plug in the resulting quantile function.
#'
#' The list of parameters should be constructed such that when the distribution function is called with the parameters, it returns a vector of the same length as the predictions. In most cases the parameters should ensure that the predicted value corresponds to the mean, median, or mode of the resulting distribution. Parameters relating to the prediction error should be estimated on calibration data. For example, if normal prediction intervals are desired, the mean parameter should be the predicted value and the standard deviation parameter should be the estimated standard deviation of the prediction errors in the calibration set. If the distribution is a negative binomial distribution with a fixed size parameter, the size parameter should be estimated on the calibration data and the mu parameter should be the predicted value.
#' @return A tibble with the predicted values and the lower and upper bounds of the prediction intervals
#' @export
#'
#' @examples
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
#' # Normal prediction intervals
#' pinterval_parametric(pred = pred_test,
#' dist = 'norm',
#' pars = list(mean = pred_test,
#'							sd = sqrt(mean((calib - calib_truth)^2))))
#'
#' # Log-normal prediction intervals
#' pinterval_parametric(pred = pred_test,
#' dist = 'lnorm',
#' pars = list(meanlog = pred_test,
#'						sdlog = sqrt(mean((log(calib) - log(calib_truth))^2))))
#'
pinterval_parametric <- function(pred,
																 dist = c('norm',
																 				 'lnorm',
																 				 'pois',
																 				 'nbinom',
																 				 'gamma',
																 				 'logis',
																 				 'beta'),
																 pars = list(),
																 alpha = 0.1,
																 lower_bound = NULL,
																 upper_bound = NULL){


	if(!is.numeric(pred)){
		stop('pred must be a single number or a numeric vector')
	}

	if(length(dist) > 1){
		stop('dist must be a single distribution')
	}

	if(!is.character(dist) & !is.function(dist)){
		stop('dist must be a character string matching a distribution or a function representing a distribution')
	}

	if(!is.list(pars) | length(pars) == 0){
		stop('pars must be a list of parameters for the distribution for each prediction')
	}

	if(!is.function(dist)){
		if(substring(dist, 1, 1) == 'q'){
			dist <- match.fun(dist)
		}else{
			dist <- match.fun(paste0('q', dist))
		}
	}

	lower_bounds <- do.call(dist, args = c(list(p = alpha/2), pars))
	upper_bounds <- do.call(dist, args = c(list(p = 1-alpha/2), pars))

	lower_bounds[lower_bounds < lower_bound] <- lower_bound
	upper_bounds[upper_bounds > upper_bound] <- upper_bound

	if(length(lower_bounds) == 1){
		warning('The parameters provided to the distribution produced a vector of length 1 for the lower bound. This will be recycled to the length of the predictions.')
		lower_bounds <- rep(lower_bounds, length(pred))
	}
	if(length(upper_bounds) == 1){
		warning('The parameters provided to the distribution produced a vector of length 1 for the upper bound. This will be recycled to the length of the predictions.')
		upper_bounds <- rep(upper_bounds, length(pred))
	}

	if(length(pred) != length(lower_bounds) | length(pred) != length(upper_bounds)){
		stop("The parameters provided to the distribution do not produce the same length of lower and upper bounds as the length of the predictions. Make sure that at least one argument in pars is a vector of the same length as pred.")
	}

	return(dplyr::tibble(pred = pred,
												lower_bound = lower_bounds,
												upper_bound = upper_bounds))

}
