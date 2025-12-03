#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector row_euclidean_distance(NumericMatrix X, NumericVector v) {
	int n = X.nrow();
	int p = X.ncol();

	if (v.size() != p) {
		stop("Length of vector v (%d) must match number of columns in X (%d)", v.size(), p);
	}

	NumericVector dists(n);

	for (int i = 0; i < n; i++) {
		double sum_sq = 0.0;
		for (int j = 0; j < p; j++) {
			double diff = X(i, j) - v[j];
			sum_sq += diff * diff;
		}
		dists[i] = std::sqrt(sum_sq);
	}

	return dists;
}

// [[Rcpp::export]]
NumericVector row_mahalanobis_distance(NumericMatrix X,
                                       NumericVector v,
                                       NumericMatrix S_inv) {
  int n = X.nrow();
  int p = X.ncol();

  if (v.size() != p) {
    stop("Length of v must match number of columns in X.");
  }
  if (S_inv.nrow() != p || S_inv.ncol() != p) {
    stop("S_inv must be a square matrix with dimensions equal to ncol(X).");
  }

  NumericVector dists(n);

  // Temporary storage for differences and S_inv * diff
  std::vector<double> diff(p);
  std::vector<double> Sinv_diff(p);

  for (int i = 0; i < n; i++) {
    // diff = X(i, :) - v
    for (int j = 0; j < p; j++) {
      diff[j] = X(i, j) - v[j];
    }

    // Sinv_diff = S_inv %*% diff
    for (int j = 0; j < p; j++) {
      double acc = 0.0;
      // Row j of S_inv times diff
      for (int k = 0; k < p; k++) {
        acc += S_inv(j, k) * diff[k];
      }
      Sinv_diff[j] = acc;
    }

    // quad = diff^T * (S_inv * diff)
    double quad = 0.0;
    for (int j = 0; j < p; j++) {
      quad += diff[j] * Sinv_diff[j];
    }

    dists[i] = std::sqrt(quad);
  }

  return dists;
}

