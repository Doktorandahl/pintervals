% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{ncs_compute}
\alias{ncs_compute}
\title{Non-Conformity Score Computation Function}
\usage{
ncs_compute(type, pred, truth, coefs = NULL)
}
\arguments{
\item{type}{Type of non-conformity score to compute. Options include 'absolute_error', 'raw_error', 'relative_error', 'relative_error2', and 'heterogeneous_error'.}

\item{pred}{a numeric vector of predicted values}

\item{truth}{a numeric vector of true values}

\item{coefs}{a numeric vector of coefficients for the heterogeneous error model. Must be of length 2, where the first element is the intercept and the second element is the slope.}
}
\description{
Non-Conformity Score Computation Function
}
