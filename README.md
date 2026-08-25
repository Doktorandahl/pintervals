# pintervals

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/pintervals)](https://CRAN.R-project.org/package=pintervals)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md)
<!-- badges: end -->

**pintervals** is an R package for computing model-agnostic prediction intervals. It provides a unified, easy-to-use interface for generating such intervals using inductive conformal prediction, bootstrapping, and parametric methods. The model-agnostic nature of the functions allows the user to quantify the uncertainty around predictions from *any* model that produces point predictions — linear models, machine learning algorithms, or anything else.

## Why pintervals?

Point predictions on their own don't tell you how much to trust them. **pintervals** lets you attach calibrated uncertainty to predictions without being tied to a specific modelling framework. Prediction intervals are generated with the `pinterval_*` functions, where each function maps to a specific framework for generating the predictions. **pintervals** generally require an out-of-sample calibration data set on which errors can be computed.

- **Conformal prediction** (`pinterval_conformal()`) — distribution-free intervals with finite-sample marginal coverage guarantees under the assumption of exchangeability. **pintervals** implements the split/inductive conformal prediction algorithm.
- **Mondrian conformal prediction** (`pinterval_mondrian()`) — group-conditional coverage for heterogeneous data with finite-sample marginal coverage guarantees within each group, under the assumption of exchangeability.
- **Clustered conformal prediction** (`pinterval_ccp()`) — pools similar groups together for more efficient group-conditional intervals when groups are numerous or small.
- **Bin-conditional conformal prediction** (`pinterval_bccp()`) — adaptive conformal prediction intervals conditional on the outcome, for skewed or imbalanced outcome distributions with finite-sample marginal coverage guarantees within binned ranges of the outcome.
- **Distance-weighted conformal prediction** — available for all conformal methods above, giving more weight to calibration points close to the new observation.
- **Bootstrapped prediction intervals** (`pinterval_bootstrap()`) — resampling-based intervals that make no distributional assumptions.
- **Parametric prediction intervals** (`pinterval_parametric()`) — fast, interpretable intervals under a specified (or custom) error distribution.

All `pinterval_*()` functions share a common interface: point predictions for a test set, plus predicted and true values from a calibration set, in — a tibble with the resulting lower and upper bounds out. Functions to evaluate interval performance, such as `interval_coverage()`, `interval_width()`, and `interval_score()`, are also included.

## Installation

Install the released version from CRAN:

```r
install.packages("pintervals")
```

Or install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("Doktorandahl/pintervals")
```

## Quick example

```r
library(pintervals)

# Split predicted/true outcome pairs into calibration and test sets
data("county_turnout", package = "pintervals")
set.seed(1)
calib_idx <- sample(nrow(county_turnout), size = 0.5 * nrow(county_turnout))
calib_data <- county_turnout[calib_idx, ]
test_data  <- county_turnout[-calib_idx, ]

# Conformal prediction intervals at 90% coverage
intervals <- pinterval_conformal(
  pred = test_data$predicted_turnout,
  calib = calib_data$predicted_turnout,
  calib_truth = calib_data$turnout,
  alpha = 0.1
)

# Check empirical coverage and average width
interval_coverage(
  truth = test_data$turnout,
  lower_bound = intervals$lower_bound,
  upper_bound = intervals$upper_bound
)
```

## Learn more

For a full introduction to the package, the theory behind each method, and a worked example predicting U.S. county-level voter turnout, see the [package vignette](https://cran.r-project.org/package=pintervals/vignettes/pintervals.html):

```r
vignette("pintervals", package = "pintervals")
```

## Citation

If you use **pintervals** in your research, please cite:

> Randahl, D., Hjort, A., & Williams, J. P. (2026). pintervals: an R package for model-agnostic prediction intervals. arXiv preprint arXiv:2601.03994.

## License

GPL (>= 3). See [LICENSE.md](LICENSE.md) for details.
