% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{bcss_compute}
\alias{bcss_compute}
\title{Function to compute the between-cluster sum of squares (BCSS) for a set of clusters}
\usage{
bcss_compute(ncs, class_vec, clusters, q = seq(0.1, 0.9, by = 0.1))
}
\arguments{
\item{ncs}{Vector of non-conformity scores}

\item{class_vec}{Vector of class labels}

\item{clusters}{List of clusters, where each element is a vector of class labels assigned to that cluster}

\item{q}{Quantiles to use for the qECDFs, default is a sequence from 0.1 to 0.9 in steps of 0.1}
}
\value{
A numeric value representing the BCSS for the clusters
}
\description{
Function to compute the between-cluster sum of squares (BCSS) for a set of clusters
}
