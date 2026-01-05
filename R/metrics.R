#' Empirical coverage of prediction intervals
#'
#' @description Calculates the mean empirical coverage rate of prediction intervals, i.e., the proportion of true values that fall within their corresponding prediction intervals.
#'
#' @param truth A numeric vector of true outcome values.
#' @param lower_bound A numeric vector of lower bounds of the prediction intervals.
#' @param upper_bound A numeric vector of upper bounds of the prediction intervals.
#' @param intervals Alternative input for prediction intervals as a list-column, where each element is a list with components 'lower_bound' and 'upper_bound'. Useful with non-contigous intervals, for instance constructed using the bin conditional conformal method wich can yield multiple intervals per prediction. See details.
#' @param return_vector Logical, whether to return the coverage vector (TRUE) or the mean coverage (FALSE). Default is FALSE.
#' @param na.rm Logical, whether to remove NA values before calculation. Default is FALSE.
#'
#' @details
#' If the `intervals` argument is provided, it should be a list-column where each element is a list containing 'lower_bound' and 'upper_bound' vectors. This allows for the calculation of coverage for non-contiguous intervals, such as those produced by certain conformal prediction methods such as the bin conditional conformal method. In this case, coverage is determined by checking if the true value falls within any of the specified intervals for each observation. If the user has some observations with contiguous intervals and others with non-contiguous intervals, they can provide both `lower_bound` and `upper_bound` vectors along with the `intervals` list-column. The function will compute coverage accordingly for each observation based on the available information.
#'
#' @return A single numeric value between 0 and 1 representing the proportion of covered values.
#'
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
#' y <- rnorm(1000, mean = x1 + x2, sd = 1)
#' df <- tibble(x1, x2, y)
#'
#' # Split into training, calibration, and test sets
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model on the log-scale
#' mod <- lm(y ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' pred_cal <- predict(mod, newdata = df_cal)
#' pred_test <- predict(mod, newdata = df_test)
#'
#' # Estimate normal prediction intervals from calibration data
#' intervals <- pinterval_parametric(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   dist = "norm",
#'   alpha = 0.1
#' )
#'
#' # Calculate empirical coverage
#' interval_coverage(truth = df_test$y,
#'          lower_bound = intervals$lower_bound,
#'          upper_bound = intervals$upper_bound)
#'
interval_coverage <- function(
	truth,
	lower_bound = NULL,
	upper_bound = NULL,
	intervals = NULL,
	return_vector = FALSE,
	na.rm = FALSE
) {
	if (!is.null(intervals)) {
		intervals_idx <- which(purrr::map_lgl(intervals, ~ !is.null(.x)))
		covered_interval <- foreach::foreach(i = intervals_idx) %do%
			{
				cover(truth[i], intervals[[i]])
			}
	}

	if (!is.null(lower_bound) & !is.null(upper_bound)) {
		# Calculate coverage
		covered <- (truth >= lower_bound) & (truth <= upper_bound)
	}

	if (is.null(intervals)) {
		out <- covered
	} else if (is.null(lower_bound) | is.null(upper_bound)) {
		out <- unlist(covered_interval)
	} else {
		out <- rep(NA, length(truth))
		out[intervals_idx] <- unlist(covered_interval)
		out[setdiff(1:length(truth), intervals_idx)] <- covered[setdiff(
			1:length(truth),
			intervals_idx
		)]
	}

	if (return_vector) {
		return(out)
	} else {
		return(mean(out, na.rm = na.rm))
	}
}


cover <- function(truth, intervals) {
	lower_bound <- intervals$lower_bound
	upper_bound <- intervals$upper_bound
	covered <- (truth >= lower_bound) & (truth <= upper_bound)
	return(max(covered))
}

#' Empirical miscoverage of prediction intervals
#'
#' @description Calculates the empirical miscoverage rate of prediction intervals, i.e., the difference between proportion of true values that fall within their corresponding prediction intervals and the nominal coverage rate (1 - alpha).
#'
#' @param truth A numeric vector of true outcome values.
#' @param lower_bound A numeric vector of lower bounds of the prediction intervals.
#' @param upper_bound A numeric vector of upper bounds of the prediction intervals.
#' @param alpha The nominal miscoverage rate (e.g., 0.1 for 90\% prediction intervals).
#' @param na.rm Logical, whether to remove NA values before calculation. Default is FALSE.
#'
#' @return A single numeric value between -1 and 1 representing the empirical miscoverage rate. A value close to 0 indicates that the prediction intervals are well-calibrated.
#'
#'
#'
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
#' y <- rnorm(1000, mean = x1 + x2, sd = 1)
#' df <- tibble(x1, x2, y)
#'
#' # Split into training, calibration, and test sets
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model on the log-scale
#' mod <- lm(y ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' pred_cal <- predict(mod, newdata = df_cal)
#' pred_test <- predict(mod, newdata = df_test)
#'
#' # Estimate normal prediction intervals from calibration data
#' intervals <- pinterval_parametric(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   dist = "norm",
#'   alpha = 0.1
#' )
#'
#' # Calculate empirical coverage
#' interval_miscoverage(truth = df_test$y,
#'          lower_bound = intervals$lower_bound,
#'          upper_bound = intervals$upper_bound,
#'          alpha = 0.1)
#'
interval_miscoverage <- function(
	truth,
	lower_bound,
	upper_bound,
	alpha,
	na.rm = FALSE
) {
	# Check if the lengths of the vectors are equal
	if (
		length(truth) != length(lower_bound) || length(truth) != length(upper_bound)
	) {
		stop("All input vectors must have the same length.")
	}

	# Calculate empirical coverage
	covered <- (truth >= lower_bound) & (truth <= upper_bound)

	# Return the proportion of covered values
	return(mean(covered, na.rm = na.rm) - (1 - alpha))
}

#' Mean interval score (MIS) for prediction intervals
#'
#' @description Computes the mean interval score, a proper scoring rule that penalizes both the width of prediction intervals and any lack of coverage. Lower values indicate better interval quality.
#' @param truth A numeric vector of true outcome values.
#' @param lower_bound A numeric vector of lower bounds of the prediction intervals.
#' @param upper_bound A numeric vector of upper bounds of the prediction intervals.
#' @param alpha The nominal miscoverage rate (e.g., 0.1 for 90\% prediction intervals).
#' @param intervals Alternative input for prediction intervals as a list-column, where each element is a list with components 'lower_bound' and 'upper_bound'. Useful with non-contigous intervals, for instance constructed using the bin conditional conformal method wich can yield multiple intervals per prediction. See details.
#' @param return_vector Logical, whether to return the interval score vector (TRUE) or the mean interval score (FALSE). Default is FALSE.
#' @param na.rm Logical, whether to remove NA values before calculation. Default is FALSE.
#'
#' @details
#'
#' The mean interval score (MIS) is defined as:
#' \deqn{
#' MIS = (ub - lb) + \frac{2}{\alpha}(lb - y) \cdot 1_{y < lb} + \frac{2}{\alpha}(y - ub) \cdot 1_{y > ub}
#' }
#' where \( y \) is the true value, and \( [lb, ub] \) is the prediction interval.
#'
#' If the `intervals` argument is provided, it should be a list-column where each element is a list containing 'lower_bound' and 'upper_bound' vectors. This allows for the calculation of coverage for non-contiguous intervals, such as those produced by certain conformal prediction methods such as the bin conditional conformal method. In this case, coverage is determined by checking if the true value falls within any of the specified intervals for each observation. If the user has some observations with contiguous intervals and others with non-contiguous intervals, they can provide both `lower_bound` and `upper_bound` vectors along with the `intervals` list-column. The function will compute coverage accordingly for each observation based on the available information.
#'
#' @return A single numeric value representing the mean interval score across all observations.
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
#' y <- rnorm(1000, mean = x1 + x2, sd = 1)
#' df <- tibble(x1, x2, y)
#'
#' # Split into training, calibration, and test sets
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model on the log-scale
#' mod <- lm(y ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' pred_cal <- predict(mod, newdata = df_cal)
#' pred_test <- predict(mod, newdata = df_test)
#'
#' # Estimate normal prediction intervals from calibration data
#' intervals <- pinterval_parametric(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   dist = "norm",
#'   alpha = 0.1
#' )
#'
#' # Calculate empirical coverage
#' interval_score(truth = df_test$y,
#'          lower_bound = intervals$lower_bound,
#'          upper_bound = intervals$upper_bound,
#'          alpha = 0.1)
interval_score <- function(
	truth,
	lower_bound = NULL,
	upper_bound = NULL,
	intervals = NULL,
	return_vector = FALSE,
	alpha,
	na.rm = FALSE
) {
	if (!is.null(intervals)) {
		intervals_idx <- which(purrr::map_lgl(intervals, ~ !is.null(.x)))
		mis_intervals <- foreach::foreach(i = intervals_idx) %do%
			{
				interval_score_interval(intervals[[i]], truth[i], alpha)
			}
	}

	if (!is.null(lower_bound) & !is.null(upper_bound)) {
		# Calculate mean interval score (MIS)
		mis_values <- (upper_bound - lower_bound) +
			2 / alpha * abs(lower_bound - truth) * (truth < lower_bound) +
			2 / alpha * abs(truth - upper_bound) * (truth > upper_bound)
	}

	if (is.null(intervals)) {
		out <- mis_values
	} else if (is.null(lower_bound) | is.null(upper_bound)) {
		out <- unlist(mis_intervals)
	} else {
		out <- rep(NA, length(truth))
		out[intervals_idx] <- unlist(mis_intervals)
		out[setdiff(1:length(truth), intervals_idx)] <- mis_values[setdiff(
			1:length(truth),
			intervals_idx
		)]
	}

	if (return_vector) {
		return(out)
	} else {
		return(mean(out, na.rm = na.rm))
	}
}


interval_score_interval <- function(intervals, truth, alpha, na.rm = FALSE) {
	lower_bound <- intervals$lower_bound
	upper_bound <- intervals$upper_bound
	if (any((truth >= lower_bound) & (truth <= upper_bound))) {
		return(sum(upper_bound - lower_bound))
	} else {
		min_penalty <- min(c(
			2 / alpha * abs(lower_bound - truth)[truth < lower_bound],
			2 / alpha * abs(truth - upper_bound)[truth > upper_bound]
		))
		return(sum(upper_bound - lower_bound) + min_penalty)
	}
}


#' Mean width of prediction intervals
#' @description Computes the mean width of prediction intervals, defined as the average difference between upper and lower bounds.
#' @param lower_bound A numeric vector of lower bounds of the prediction intervals.
#' @param upper_bound A numeric vector of upper bounds of the prediction intervals.
#' @param intervals Alternative input for prediction intervals as a list-column, where each element is a list with components 'lower_bound' and 'upper_bound'. Useful with non-contigous intervals, for instance constructed using the bin conditional conformal method wich can yield multiple intervals per prediction. See details.
#' @param return_vector Logical, whether to return the width vector (TRUE) or the mean width (FALSE). Default is FALSE.
#' @param na.rm Logical, whether to remove NA values before calculation. Default is FALSE.
#' @details
#' The mean width is calculated as:
#' \deqn{
#' \text{Mean Width} = \frac{1}{n} \sum_{i=1}^{n} (ub_i - lb_i)
#' }
#'
#' where \( ub_i \) and \( lb_i \) are the upper and lower bounds of the prediction interval for observation \( i \), and \( n \) is the total number of observations.
#'
#' If the `intervals` argument is provided, it should be a list-column where each element is a list containing 'lower_bound' and 'upper_bound' vectors. This allows for the calculation of coverage for non-contiguous intervals, such as those produced by certain conformal prediction methods such as the bin conditional conformal method. In this case, coverage is determined by checking if the true value falls within any of the specified intervals for each observation. If the user has some observations with contiguous intervals and others with non-contiguous intervals, they can provide both `lower_bound` and `upper_bound` vectors along with the `intervals` list-column. The function will compute coverage accordingly for each observation based on the available information.
#'
#' @return A single numeric value representing the mean width of the prediction intervals.
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
#' y <- rnorm(1000, mean = x1 + x2, sd = 1)
#' df <- tibble(x1, x2, y)
#'
#' # Split into training, calibration, and test sets
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit a model on the log-scale
#' mod <- lm(y ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' pred_cal <- predict(mod, newdata = df_cal)
#' pred_test <- predict(mod, newdata = df_test)
#'
#' # Estimate normal prediction intervals from calibration data
#' intervals <- pinterval_parametric(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   dist = "norm",
#'   alpha = 0.1
#' )
#'
#' # Calculate empirical coverage
#' interval_width(lower_bound = intervals$lower_bound,
#'          upper_bound = intervals$upper_bound)
interval_width <- function(
	lower_bound = NULL,
	upper_bound = NULL,
	intervals = NULL,
	return_vector = FALSE,
	na.rm = FALSE
) {
	if (!is.null(intervals)) {
		intervals_idx <- which(purrr::map_lgl(intervals, ~ !is.null(.x)))
		width_intervals <- foreach::foreach(i = intervals_idx) %do%
			{
				width_interval(intervals[[i]])
			}
	}

	if (!is.null(lower_bound) & !is.null(upper_bound)) {
		# Calculate the width of the intervals
		widths <- upper_bound - lower_bound
	}
	if (is.null(intervals)) {
		out <- widths
	} else if (is.null(lower_bound) | is.null(upper_bound)) {
		out <- unlist(width_intervals)
	} else {
		out <- rep(NA, length(lower_bound))
		out[intervals_idx] <- unlist(width_intervals)
		out[setdiff(1:length(lower_bound), intervals_idx)] <- widths[setdiff(
			1:length(lower_bound),
			intervals_idx
		)]
	}

	if (return_vector) {
		return(out)
	} else {
		return(mean(out, na.rm = na.rm))
	}
}


width_interval <- function(intervals) {
	lower_bound <- intervals$lower_bound
	upper_bound <- intervals$upper_bound
	return(sum(upper_bound - lower_bound))
}
