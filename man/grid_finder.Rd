% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{grid_finder}
\alias{grid_finder}
\title{Grid search for lower and upper bounds of continuous conformal prediction intervals}
\usage{
grid_finder(
  y_min,
  y_max,
  ncs,
  ncs_type,
  y_hat,
  alpha,
  min_step = NULL,
  grid_size = NULL,
  calib = NULL,
  coefs = NULL,
  distance_weighted_cp = FALSE,
  distance_features_calib = NULL,
  distance_features_pred = NULL,
  normalize_distance = TRUE,
  weight_function = gauss_kern
)
}
\arguments{
\item{y_min}{minimum value to search}

\item{y_max}{maximum value to search}

\item{ncs}{vector of non-conformity scores}

\item{ncs_type}{String indicating the non-conformity score function to use}

\item{y_hat}{vector of predicted values}

\item{alpha}{confidence level}

\item{min_step}{The minimum step size for the grid search}

\item{grid_size}{Alternative to min_step, the number of points to use in the grid search between the lower and upper bound}

\item{calib}{a tibble with the predicted values and the true values of the calibration partition. Used when weighted_cp is TRUE. Default is NULL}

\item{coefs}{a numeric vector of coefficients for the heterogeneous error model. Must be of length 2, where the first element is the intercept and the second element is the slope. Used when ncs_type is 'heterogeneous_error'. Default is NULL}

\item{distance_weighted_cp}{logical. If TRUE, the non-conformity scores will be weighted according to the distance function}

\item{distance_features_calib}{a matrix of features for the calibration partition. Used when distance_weighted_cp is TRUE}

\item{distance_features_pred}{a matrix of features for the prediction partition. Used when distance_weighted_cp is TRUE}

\item{normalize_distance}{logical. If TRUE, the distances will be normalized to the range [0,1]}

\item{weight_function}{a function to use for weighting the distances. Can be 'gaussian_kernel', 'caucy_kernel', 'logistic', or 'reciprocal_linear'. Default is 'gaussian_kernel'}
}
\value{
a tibble with the predicted values and the lower and upper bounds of the prediction intervals
}
\description{
Grid search for lower and upper bounds of continuous conformal prediction intervals
}
