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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//


