% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pinterval_bccp.R
\name{pinterval_bccp}
\alias{pinterval_bccp}
\title{Bin-conditional conformal prediction intervals for continuous predictions}
\usage{
pinterval_bccp(
  pred,
  calib = NULL,
  calib_truth = NULL,
  calib_bins = NULL,
  breaks = NULL,
  alpha = 0.1,
  ncs_type = c("absolute_error", "relative_error", "za_relative_error",
    "heterogeneous_error", "raw_error"),
  distance_weighted_cp = FALSE,
  distance_features_calib = NULL,
  distance_features_pred = NULL,
  normalize_distance = TRUE,
  weight_function = c("gaussian_kernel", "caucy_kernel", "logistic", "reciprocal_linear"),
  calibrate = FALSE,
  calibration_method = "glm",
  calibration_family = NULL,
  calibration_transform = NULL,
  resolution = 0.01,
  grid_size = NULL,
  right = TRUE,
  contiguize = FALSE
)
}
\arguments{
\item{pred}{Vector of predicted values}

\item{calib}{A numeric vector of predicted values in the calibration partition or a 2 or 3 column tibble or matrix with the first column being the predicted values and the second column being the truth values and the third, optional, column being the bin labels. If calib is a numeric vector, calib_truth and either calib_bins or breaks must be provided.}

\item{calib_truth}{A numeric vector of true values in the calibration partition}

\item{calib_bins}{A vector of bin identifiers for the calibration set. Not used if breaks are provided.}

\item{breaks}{A vector of break points for the bins to manually define the bins. If NULL, lower and upper bounds of the bins are calculated as the minimum and maximum values of each bin in the calibration set. Must be provided if calib_bins is not provided, either as a vector or as the last column of a calib tibble.}

\item{alpha}{The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1}

\item{ncs_type}{A string specifying the type of nonconformity score to use. Available options are:
\itemize{
  \item \code{"absolute_error"}: \eqn{|y - \hat{y}|}
  \item \code{"relative_error"}: \eqn{|y - \hat{y}| / \hat{y}}
  \item \code{"zero_adjusted_relative_error"}: \eqn{|y - \hat{y}| / (\hat{y} + 1)}
  \item \code{"heterogeneous_error"}: \eqn{|y - \hat{y}| / \sigma_{\hat{y}}} absolute error divided by a measure of heteroskedasticity, computed as the predicted value from a linear model of the absolute error on the predicted values
  \item \code{"raw_error"}: the signed error \eqn{y - \hat{y}}
}
The default is \code{"absolute_error"}.}

\item{distance_weighted_cp}{Logical. If \code{TRUE}, weighted conformal prediction is performed where the non-conformity scores are weighted based on the distance between calibration and prediction points in feature space. Default is \code{FALSE}.}

\item{distance_features_calib}{A matrix, data frame, or numeric vector of features from which to compute distances when \code{distance_weighted_cp = TRUE}. This should contain the feature values for the calibration set. Must have the same number of rows as the calibration set. Can be the predicted values themselves, or any other features which give a meaningful distance measure.}

\item{distance_features_pred}{A matrix, data frame, or numeric vector of feature values for the prediction set. Must be the same features as specified in \code{distance_features_calib}. Required if \code{distance_weighted_cp = TRUE}.}

\item{normalize_distance}{Logical. If \code{TRUE}, distances are normalized to the [0, 1] interval before applying the weight function. This is typically recommended to ensure consistent scaling across features. Default is \code{TRUE}.}

\item{weight_function}{A character string specifying the weighting kernel to use for distance-weighted conformal prediction. Options are:
\itemize{
  \item \code{"gaussian_kernel"}: \eqn{ w(d) = e^{-d^2} }
  \item \code{"caucy_kernel"}: \eqn{ w(d) = 1/(1 + d^2) }
  \item \code{"logistic"}: \eqn{ w(d) = 1//(1 + e^{d}) }
  \item \code{"reciprocal_linear"}: \eqn{ w(d) = 1/(1 + d) }
}
The default is \code{"gaussian_kernel"}. Distances are computed as the Euclidean distance between the calibration and prediction feature vectors.}

\item{calibrate}{= FALSE Logical. If TRUE, the function will calibrate the predictions and intervals using the calibration set. Default is FALSE. See details for more information.}

\item{calibration_method}{The method to use for calibration. Can be "glm" or "isotonic". Default is "glm". Only used if calibrate = TRUE.}

\item{calibration_family}{The family used for the calibration model. Default is "gaussian". Only used if calibrate = TRUE and calibration_method = "glm".}

\item{calibration_transform}{Optional transformation to apply to the predictions before calibration. Default is NULL. Only used if calibrate = TRUE and calibration_method = "glm".}

\item{resolution}{The minimum step size for the grid search. Default is 0.01. See details for more information.}

\item{grid_size}{Alternative to `resolution`, the number of points to use in the grid search between the lower and upper bound. If provided, resolution will be ignored.}

\item{right}{Logical, if TRUE the bins are right-closed (a,b] and if FALSE the bins are left-closed `[ a,b)`. Only used if breaks are provided.}

\item{contiguize}{logical indicating whether to contiguize the intervals. TRUE will consider all bins for each prediction using the lower and upper endpoints as interval limits to avoid non-contiguous intervals. FALSE will allows for non-contiguous intervals. TRUE guarantees at least appropriate coverage in each bin, but may suffer from over-coverage in certain bins. FALSE will have appropriate coverage in each bin but may have non-contiguous intervals. Default is FALSE.}
}
\value{
A tibble with the predicted values, the lower and upper bounds of the prediction intervals. If treat_noncontiguous is 'non_contiguous', the lower and upper bounds are set in a list variable called 'intervals' where all non-contiguous intervals are stored.
}
\description{
This function calculates bin-conditional conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction on a bin-by-bin basis. The intervals are computed using a calibration set with predicted and true values and their associated bins. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals. Bin-conditional conformal prediction intervals are useful when the prediction error is not constant across the range of predicted values and ensures that the coverage is (approximately) correct for each bin under the assumption that the non-conformity scores are exchangeable within each bin.
}
\details{
This function computes bin-conditional conformal prediction intervals using inductive conformal prediction applied separately within each bin of the calibration data. It is particularly useful when prediction error varies across the range of predicted values, as it enables locally valid coverage by ensuring that the coverage level \eqn{1 - \eqn{\alpha}} holds within each bin—assuming exchangeability of non-conformity scores within bins.

The calibration set must include predicted values, true values, and corresponding bin identifiers or breaks for the bins. These can be provided either as separate vectors (`calib`, `calib_truth`, and `calib_bins` or `breaks`).

Non-conformity scores are calculated using the specified `ncs_type`, which determines how the prediction error is measured. Available options include:

- `"absolute_error"`: the absolute difference between predicted and true values.
- `"relative_error"`: the absolute error divided by the true value.
- `"za_relative_error"`: zero-adjusted relative error, which replaces small or zero true values with a small constant to avoid division by zero.
- `"heterogeneous_error"`: absolute error scaled by a linear model of prediction error magnitude as a function of the predicted value.
- `"raw_error"`: the signed difference between predicted and true values.

These options provide flexibility to adapt to different patterns of prediction error across the outcome space.

Bins endpoints can be defined manually via the `breaks` argument or inferred from the calibration data. If `contiguize = TRUE`, the function ensures the resulting prediction intervals are contiguous across bins, potentially increasing coverage beyond the nominal level in some bins. If `contiguize = FALSE`, the function may produce non-contiguous intervals, which are more efficient but may be harder to interpret.

The prediction intervals are constructed using a grid search over a user-defined range of outcome values. The resolution of the grid search can be controlled by either the `resolution` argument, which sets the minimum step size, or the `grid_size` argument, which sets the number of grid points. For wide prediction spaces, the grid search may be computationally intensive. In such cases, increasing the `resolution` or reducing the `grid_size` may improve performance.

Optionally, the predicted values can be calibrated before interval construction by setting `calibrate = TRUE`. In this case, the predictions are passed through `calibrate_predictions()` to adjust the predictions based on the calibration set. The calibration method can be specified using `calibration_method` and `calibration_family`, with "glm" being the default method. See \link[pintervals]{calibrate_predictions} for more information on calibration.

The function returns a tibble containing the predicted values and their corresponding lower and upper bounds. If `contiguize = FALSE` and the resulting intervals are non-contiguous, the lower and upper bounds are in list-columns.
}
\examples{

# Generate example data
library(dplyr)
library(tibble)
x1 <- runif(1000)
x2 <- runif(1000)
y <- rlnorm(1000, meanlog = x1 + x2, sdlog = 0.5)

# Create bins based on quantiles
	bin <- cut(y, breaks = quantile(y, probs = seq(0, 1, 1/4)),
	include.lowest = TRUE, labels =FALSE)
df <- tibble(x1, x2, y, bin)
df_train <- df \%>\% slice(1:500)
df_cal <- df \%>\% slice(501:750)
df_test <- df \%>\% slice(751:1000)

# Fit a model to the training data
mod <- lm(log(y) ~ x1 + x2, data=df_train)

# Generate predictions on the original y scale for the calibration data
calib <- exp(predict(mod, newdata=df_cal))
calib_truth <- df_cal$y
calib_bins <- df_cal$bin

# Generate predictions for the test data

pred_test <- exp(predict(mod, newdata=df_test))

# Calculate bin-conditional conformal prediction intervals
pinterval_bccp(pred = pred_test,
calib = calib,
calib_truth = calib_truth,
calib_bins = calib_bins,
alpha = 0.1)

}
