#' Empirical coverage of prediction intervals
#'
#' @description Calculates the mean empirical coverage rate of prediction intervals, i.e., the proportion of true values that fall within their corresponding prediction intervals.
#'
#' @param truth A numeric vector of true outcome values.
#' @param lower_bound A numeric vector of lower bounds of the prediction intervals.
#' @param upper_bound A numeric vector of upper bounds of the prediction intervals.
#' @param na.rm Logical, whether to remove NA values before calculation. Default is FALSE.
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
interval_coverage <- function(truth, lower_bound, upper_bound,na.rm=FALSE) {
  # Check if the lengths of the vectors are equal
  if (length(truth) != length(lower_bound) || length(truth) != length(upper_bound)) {
    stop("All input vectors must have the same length.")
  }

  # Calculate coverage
  covered <- (truth >= lower_bound) & (truth <= upper_bound)

  # Return the proportion of covered values
  return(mean(covered,na.rm=na.rm))
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
interval_miscoverage <- function(truth, lower_bound, upper_bound, alpha,na.rm=FALSE) {
  # Check if the lengths of the vectors are equal
  if (length(truth) != length(lower_bound) || length(truth) != length(upper_bound)) {
    stop("All input vectors must have the same length.")
  }

  # Calculate empirical coverage
  covered <- (truth >= lower_bound) & (truth <= upper_bound)

  # Return the proportion of covered values
  return(mean(covered,na.rm=na.rm) - (1 - alpha))
}

#' Mean interval score (MIS) for prediction intervals
#'
#' @description Computes the mean interval score, a proper scoring rule that penalizes both the width of prediction intervals and any lack of coverage. Lower values indicate better interval quality.
#' @param truth A numeric vector of true outcome values.
#' @param lower_bound A numeric vector of lower bounds of the prediction intervals.
#' @param upper_bound A numeric vector of upper bounds of the prediction intervals.
#' @param alpha The nominal miscoverage rate (e.g., 0.1 for 90\% prediction intervals).
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
interval_score <- function(truth, lower_bound, upper_bound,alpha,na.rm=FALSE) {
  # Check if the lengths of the vectors are equal
  if (length(truth) != length(lower_bound) || length(truth) != length(upper_bound)) {
    stop("All input vectors must have the same length.")
  }

  # Calculate mean interval score (MIS)
	mis_values <- (upper_bound - lower_bound) +
    2/alpha * abs(lower_bound-truth) * (truth < lower_bound) +
		    2/alpha * abs(truth-upper_bound) * (truth > upper_bound)

  # Return the mean of the MIS values
  return(mean(mis_values, na.rm = na.rm))
}



#' Mean width of prediction intervals
#' @description Computes the mean width of prediction intervals, defined as the average difference between upper and lower bounds.
#' @param lower_bound A numeric vector of lower bounds of the prediction intervals.
#' @param upper_bound A numeric vector of upper bounds of the prediction intervals.
#' @param na.rm Logical, whether to remove NA values before calculation. Default is FALSE.
#' @details
#' The mean width is calculated as:
#' \deqn{
#' \text{Mean Width} = \frac{1}{n} \sum_{i=1}^{n} (ub_i - lb_i)
#' }
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
interval_width <- function(lower_bound, upper_bound,na.rm=FALSE) {
  # Check if the lengths of the vectors are equal
  if (length(lower_bound) != length(upper_bound)) {
    stop("lower_bound and upper_bound must have the same length.")
  }

  # Calculate the width of the intervals
  widths <- upper_bound - lower_bound

  # Return the mean width
  return(mean(widths,na.rm=na.rm))
}
