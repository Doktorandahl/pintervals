#' #' Parametric prediction intervals for continuous predictions
#'
#' This function computes parametric prediction intervals at a confidence level of \eqn{1 - \eqn{\alpha}} for a vector of continuous predictions. The intervals are based on a user-specified probability distribution and associated parameters, either estimated from calibration data or supplied directly. Supported distributions include common options like the normal, log-normal, gamma, beta, and negative binomial, as well as any user-defined distribution with a quantile function. Prediction intervals are calculated by evaluating the appropriate quantiles for each predicted value.
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values
#' @param calib_truth A numeric vector of true values in the calibration partition. Only required if calib is a numeric vector
#' @param dist Distribution to use for the prediction intervals. Can be a character string matching any available distribution in R or a function representing a distribution, e.g. `qnorm`, `qgamma`, or a user defined quantile function. Default options are 'norm', 'lnorm','exp, 'pois', 'nbinom', 'chisq', 'gamma', 'logis', and 'beta' for which parameters can be computed from the calibration set. If a custom function is provided, parameters need to be provided in `pars`.
#' @param pars List of named parameters for the distribution for each prediction. Not needed if calib is provided and the distribution is one of the default options. If a custom distribution function is provided, this list should contain the parameters needed for the quantile function, with names matching the corresponding arguments for the parameter names of the distribution function. See details for more information.
#'
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#'
#' @details
#' This function supports a wide range of distributions for constructing prediction intervals. Built-in support is provided for the following distributions: `"norm"`, `"lnorm"`, `"exp"`, `"pois"`, `"nbinom"`, `"chisq"`, `"gamma"`, `"logis"`, and `"beta"`. For each of these, parameters can be automatically estimated from a calibration set if not supplied directly via the `pars` argument.
#'
#' The calibration set (`calib` and `calib_truth`) is used to estimate error dispersion or shape parameters. For example:
#' - **Normal**: standard deviation of errors
#' - **Log-normal**: standard deviation of log-errors
#' - **Gamma**: dispersion via `glm`
#' - **Negative binomial**: dispersion via `glm.nb()`
#' - **Beta**: precision estimated from error variance
#'
#' If `pars` is supplied, it should be a list of named arguments corresponding to the distribution’s quantile function. Parameters may be scalars or vectors (one per prediction). When both `pars` and `calib` are provided, the values in `pars` are used.
#'
#' Users may also specify a custom distribution by passing a quantile function directly (e.g., a function with the signature `function(p, ...)`) as the `dist` argument, in which case `pars` must be provided explicitly.
#'
#' @return A tibble with the predicted values and the lower and upper bounds of the prediction intervals
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#'
#' # Simulate example data
#' set.seed(123)
#' x1 <- runif(1000)
#' x2 <- runif(1000)
#' y <- rlnorm(1000, meanlog = x1 + x2, sdlog = 0.5)
#' df <- tibble(x1, x2, y)
#'
#' # Split into training, calibration, and test sets
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model on the log-scale
#' mod <- lm(log(y) ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' pred_cal <- exp(predict(mod, newdata = df_cal))
#' pred_test <- exp(predict(mod, newdata = df_test))
#'
#' # Estimate log-normal prediction intervals from calibration data
#' log_resid_sd <- sqrt(mean((log(pred_cal) - log(df_cal$y))^2))
#' pinterval_parametric(
#'   pred = pred_test,
#'   dist = "lnorm",
#'   pars = list(meanlog = log(pred_test), sdlog = log_resid_sd)
#' )
#'
#' # Alternatively, use calibration data directly to estimate parameters
#' pinterval_parametric(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   dist = "lnorm"
#' )
#'
#' # Use the normal distribution with direct parameter input
#' norm_sd <- sqrt(mean((pred_cal - df_cal$y)^2))
#' pinterval_parametric(
#'   pred = pred_test,
#'   dist = "norm",
#'   pars = list(mean = pred_test, sd = norm_sd)
#' )
#'
#' # Use the gamma distribution with parameters estimated from calibration data
#' pinterval_parametric(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   dist = "gamma"
#' )
pinterval_parametric <- function(
	pred,
	calib = NULL,
	calib_truth = NULL,
	dist = c(
		'norm',
		'lnorm',
		'exp',
		'pois',
		'nbinom',
		'gamma',
		'chisq',
		'logis',
		'beta'
	),
	pars = list(),
	alpha = 0.1
) {
	if (!is.numeric(pred)) {
		stop('pred must be a single number or a numeric vector')
	}

	if (length(dist) > 1) {
		stop('dist must be a single distribution')
	}

	if (!is.character(dist) & !is.function(dist)) {
		stop(
			'dist must be a character string matching a distribution or a function representing a distribution'
		)
	}

	if (is.null(calib) && is.null(pars) && length(pars) == 0) {
		stop('Either calib or pars must be provided.')
	}

	if (!is.null(calib) && is.numeric(calib) && is.null(calib_truth)) {
		stop('If calib is numeric, calib_truth must be provided')
	}

	if (!is.null(calib) && !is.numeric(calib) && ncol(calib) != 2) {
		stop(
			'calib must be a numeric vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values'
		)
	}

	if (length(pars) != 0 && !is.list(pars)) {
		stop(
			'pars must be a list of parameters for the distribution for each prediction'
		)
	}

	if (length(pars) != 0 && !is.null(calib)) {
		warning(
			'pars is provided, but calib is also provided. The provided parameters will be used for the prediction intervals and calib for the calibration of the predictions.'
		)
	}

	if (length(pars) == 0) {
		if (dist %in% c('norm', 'qnorm')) {
			message(
				'The distribution is a normal distribution. The standard deviation parameters will be estimated from the calibration data.'
			)
			pars$mean <- pred
			pars$sd <- sqrt(var((calib - calib_truth)))
		} else if (dist %in% c('lnorm', 'qlnorm')) {
			message(
				'The distribution is a log-normal distribution. The standard deviation parameters will be estimated from the calibration data.'
			)
			pars$meanlog <- log(pred)
			pars$sdlog <- sqrt(var((log(calib) - log(calib_truth))))
		} else if (dist %in% c('pois', 'qpois')) {
			pars$lambda <- pred
		} else if (dist %in% c('nbinom', 'qnbinom')) {
			message(
				'The distribution is a negative binomial distribution. The size (dispersion) parameter will be estimated from the calibration data. The dispersion parameter is assumed to be constant across predictions.'
			)
			pars$mu <- pred
			mle_theta <- MASS::glm.nb(calib_truth ~ offset(log(calib)))
			pars$size <- mle_theta$theta
		} else if (dist %in% c('gamma', 'qgamma')) {
			message(
				'The distribution is a gamma distribution. The rate (dispersion) parameter will be estimated from the calibration data. The dispersion parameter is assumed to be constant across predictions.'
			)

			dispersion <- summary(stats::glm(
				calib_truth ~ offset(log(calib)),
				family = stats::Gamma(link = 'log')
			))$dispersion

			pars$shape <- pred / dispersion
			pars$rate <- dispersion
		} else if (dist %in% c('logis', 'qlogis')) {
			message(
				'The distribution is a logistic distribution. The scale parameter will be estimated from the calibration data.'
			)
			pars$location <- pred
			var <- var((calib - calib_truth))
			pars$scale <- sqrt(3 * var / pi^2)
		} else if (dist %in% c('beta', 'qbeta')) {
			if (
				any(calib > 1) ||
					any(calib < 0) ||
					any(calib_truth > 1) ||
					any(calib_truth < 0)
			) {
				stop(
					'The beta distribution requires values between 0 and 1. Please ensure that calib is between 0 and 1.'
				)
			}
			message(
				'The distribution is a beta distribution. The precision parameter will be estimated from the calibration data. The precision parameter is assumed to be constant across predictions.'
			)

			var_hat <- var((calib - calib_truth))
			phi_hat <- (mean(calib) * (1 - mean(calib))) / var_hat - 1
			if (phi_hat <= 0) {
				stop(
					'The estimated precision parameter is non-positive. Please check the calibration data.'
				)
			}
			pars$shape1 <- pred * phi_hat
			pars$shape2 <- (1 - pred) * phi_hat
		} else {
			stop(
				'The distribution is not supported or not implemented yet. Please provide the parameters in pars.'
			)
		}
	}

	if (!is.function(dist)) {
		if (substring(dist, 1, 1) == 'q') {
			dist <- match.fun(dist)
		} else {
			dist <- match.fun(paste0('q', dist))
		}
	}

	lower_bounds <- do.call(dist, args = c(list(p = alpha / 2), pars))
	upper_bounds <- do.call(dist, args = c(list(p = 1 - alpha / 2), pars))

	out <- dplyr::tibble(
		pred = pred,
		lower_bound = lower_bounds,
		upper_bound = upper_bounds
	)

	attr(out, 'dist') <- dist

	return(out)
}
