% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{interval_miscoverage}
\alias{interval_miscoverage}
\title{Empirical miscoverage of prediction intervals}
\usage{
interval_miscoverage(truth, lower_bound, upper_bound, alpha)
}
\arguments{
\item{truth}{A numeric vector of true outcome values.}

\item{lower_bound}{A numeric vector of lower bounds of the prediction intervals.}

\item{upper_bound}{A numeric vector of upper bounds of the prediction intervals.}

\item{alpha}{The nominal miscoverage rate (e.g., 0.1 for 90\% prediction intervals).}
}
\value{
A single numeric value between -1 and 1 representing the empirical miscoverage rate. A value close to 0 indicates that the prediction intervals are well-calibrated.
}
\description{
Calculates the empirical miscoverage rate of prediction intervals, i.e., the difference between proportion of true values that fall within their corresponding prediction intervals and the nominal coverage rate (1 - alpha).
}
\examples{
library(dplyr)
library(tibble)

# Simulate example data
set.seed(123)
x1 <- runif(1000)
x2 <- runif(1000)
y <- rnorm(1000, mean = x1 + x2, sd = 1)
df <- tibble(x1, x2, y)

# Split into training, calibration, and test sets
df_train <- df \%>\% slice(1:500)
df_cal <- df \%>\% slice(501:750)
df_test <- df \%>\% slice(751:1000)

# Fit a model on the log-scale
mod <- lm(y ~ x1 + x2, data = df_train)

# Generate predictions
pred_cal <- predict(mod, newdata = df_cal)
pred_test <- predict(mod, newdata = df_test)

# Estimate normal prediction intervals from calibration data
intervals <- pinterval_parametric(
  pred = pred_test,
  calib = pred_cal,
  calib_truth = df_cal$y,
  dist = "norm",
  alpha = 0.1
)

# Calculate empirical coverage
coverage(truth = df_test$y,
         lower_bound = intervals$lower_bound,
         upper_bound = intervals$upper_bound)

}
