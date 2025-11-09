#' Clustered conformal prediction intervals for continuous predictions
#'
#' @description
#' This function computes conformal prediction intervals with a confidence level of \eqn{1 - \alpha} by first grouping Mondrian classes into data-driven clusters based on the distribution of their nonconformity scores. The resulting clusters are used as strata for computing class-conditional (Mondrian-style) conformal prediction intervals. This approach improves local validity and statistical efficiency when there are many small or similar classes with overlapping prediction behavior. The coverage level \eqn{1 - \alpha} is approximate within each cluster, assuming exchangeability of nonconformity scores within clusters.
#'
#' The method supports additional features such as prediction calibration, distance-weighted conformal scores, and clustering optimization via internal validity measures (e.g., Calinski-Harabasz index or minimum cluster size heuristics).
#'
#'
#' @param pred Vector of predicted values or a 2 column tibble or matrix with the first column being the predicted values and the second column being the Mondrian class labels. If pred is a numeric vector, pred_class must be provided.
#' @param pred_class A vector of class identifiers for the predicted values. This is used for sorting the predictions into their corresponding clusters.
#' @param calib A numeric vector of predicted values in the calibration partition or a 3 column tibble or matrix with the first column being the predicted values and the second column being the truth values and the third column being the Mondrian class labels. If calib is a numeric vector, calib_truth and calib_class must be provided.
#' @param calib_truth A numeric vector of true values in the calibration partition
#' @param calib_class A vector of class identifiers for the calibration set which is used for clustering the nonconformity scores. If calib is a tibble or matrix, this can be extracted from the third column.
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param ncs_type A string specifying the type of nonconformity score to use. Available options are:
#' \itemize{
#'   \item \code{"absolute_error"}: \eqn{|y - \hat{y}|}
#'   \item \code{"relative_error"}: \eqn{|y - \hat{y}| / \hat{y}}
#'   \item \code{"zero_adjusted_relative_error"}: \eqn{|y - \hat{y}| / (\hat{y} + 1)}
#'   \item \code{"heterogeneous_error"}: \eqn{|y - \hat{y}| / \sigma_{\hat{y}}} absolute error divided by a measure of heteroskedasticity, computed as the predicted value from a linear model of the absolute error on the predicted values
#'   \item \code{"raw_error"}: the signed error \eqn{y - \hat{y}}
#' }
#' The default is \code{"absolute_error"}.
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#' @param n_clusters Number of clusters to use when combining Mondrian classes. Required if \code{optimize_n_clusters = FALSE}.
#' @param cluster_method Clustering method used to group Mondrian classes. Options are \code{"kmeans"} or \code{"ks"} (Kolmogorov-Smirnov). Default is \code{"kmeans"}.
#' @param cluster_train_fraction Fraction of the calibration data used to estimate nonconformity scores and compute clustering. Default is 1 (use all).
#' @param optimize_n_clusters Logical. If \code{TRUE}, the number of clusters is chosen automatically based on internal clustering criteria.
#' @param optimize_n_clusters_method Method used for cluster optimization. One of \code{"calinhara"} (Calinski-Harabasz index) or \code{"min_cluster_size"}. Default is \code{"calinhara"}.
#' @param min_cluster_size Minimum number of calibration points per cluster. Used only when \code{optimize_n_clusters_method = "min_cluster_size"}.
#' @param min_n_clusters Minimum number of clusters to consider when optimizing.
#' @param max_n_clusters Maximum number of clusters to consider. If \code{NULL}, the upper limit is set to the number of unique Mondrian classes minus 1.
#'
#' @param distance_weighted_cp Logical. If \code{TRUE}, weighted conformal prediction is performed where the non-conformity scores are weighted based on the distance between calibration and prediction points in feature space. Default is \code{FALSE}.
#'
#' @param distance_features_calib A matrix, data frame, or numeric vector of features from which to compute distances when \code{distance_weighted_cp = TRUE}. This should contain the feature values for the calibration set. Must have the same number of rows as the calibration set. Can be the predicted values themselves, or any other features which give a meaningful distance measure.
#' @param distance_features_pred A matrix, data frame, or numeric vector of feature values for the prediction set. Must be the same features as specified in \code{distance_features_calib}. Required if \code{distance_weighted_cp = TRUE}.
#'
#' @param normalize_distance Either 'minmax', 'sd', or 'none'. Indicates if and how to normalize the distances when distance_weighted_cp is TRUE. Normalization helps ensure that distances are on a comparable scale across features. Default is 'minmax'.
#'
#' @param weight_function A character string specifying the weighting kernel to use for distance-weighted conformal prediction. Options are:
#' \itemize{
#'   \item \code{"gaussian_kernel"}: \eqn{ w(d) = e^{-d^2} }
#'   \item \code{"caucy_kernel"}: \eqn{ w(d) = 1/(1 + d^2) }
#'   \item \code{"logistic"}: \eqn{ w(d) = 1//(1 + e^{d}) }
#'   \item \code{"reciprocal_linear"}: \eqn{ w(d) = 1/(1 + d) }
#' }
#' The default is \code{"gaussian_kernel"}. Distances are computed as the Euclidean distance between the calibration and prediction feature vectors.
#'
#' @param calibrate = FALSE Logical. If TRUE, the function will calibrate the predictions and intervals using the calibration set. Default is FALSE. See details for more information on calibration.
#' @param calibration_method The method to use for calibration. Can be "glm" or "isotonic". Default is "glm". Only used if calibrate = TRUE.
#' @param calibration_family The family used for the calibration model. Default is "gaussian". Only used if calibrate = TRUE and calibration_method = "glm".
#' @param calibration_transform Optional transformation to apply to the predictions before calibration. Default is NULL. Only used if calibrate = TRUE and calibration_method = "glm".
#' @param resolution The minimum step size for the grid search. Default is 0.01. See details for more information.
#' @param grid_size Alternative to `resolution`, the number of points to use in the grid search between the lower and upper bound. If provided, resolution will be ignored.
#'
##' @return A tibble with predicted values, lower and upper prediction interval bounds, class labels, and assigned cluster labels. Attributes include clustering diagnostics (e.g., cluster assignments, coverage gaps, internal validity scores).
#'
#' @details
#' This function implements a clustered conformal prediction approach by aggregating similar Mondrian classes into clusters based on the distribution of their nonconformity scores. Each cluster is then treated as a stratum for standard inductive conformal prediction. This improves robustness and efficiency in settings with many small or noisy groups, while still enabling conditional validity.
#'
#' #' Non-conformity scores are calculated using the specified `ncs_type`, which determines how the prediction error is measured. Available options include:
#'
#' - `"absolute_error"`: the absolute difference between predicted and true values.
#' - `"relative_error"`: the absolute error divided by the true value.
#' - `"za_relative_error"`: zero-adjusted relative error, which replaces small or zero true values with a small constant to avoid division by zero.
#' - `"heterogeneous_error"`: absolute error scaled by a linear model of prediction error magnitude as a function of the predicted value.
#' - `"raw_error"`: the signed difference between predicted and true values.
#'
#' These options provide flexibility to adapt to different patterns of prediction error across the outcome space.
#'
#' Clustering is performed using the nonconformity scores of calibration data, optionally on a subsample defined by \code{cluster_train_fraction}. Users can specify the number of clusters directly, or let the function choose the optimal number based on internal criteria.
#'
#' If distance-weighted conformal prediction is enabled, calibration examples are weighted based on their similarity to test points, with several kernel functions available.
#'
#' Optionally, the predicted values can be calibrated before interval construction by setting `calibrate = TRUE`. In this case, the predictions are passed through `calibrate_predictions()` to adjust the predictions based on the calibration set. The calibration method can be specified using `calibration_method` and `calibration_family`, with "glm" being the default method. See \link[pintervals]{calibrate_predictions} for more information on calibration.
#'
#'
#' @seealso \code{\link[pintervals]{pinterval_conformal}}, \code{\link[pintervals]{pinterval_mondrian}}, \code{\link[pintervals]{calibrate_predictions}}
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#'
#' # Simulate data with 6 Mondrian classes forming 3 natural clusters
#' set.seed(123)
#' x1 <- runif(1000)
#' x2 <- runif(1000)
#' class_raw <- sample(1:6, size = 1000, replace = TRUE)
#'
#' # Construct 3 latent clusters: (1,2), (3,4), (5,6)
#' mu <- ifelse(class_raw %in% c(1, 2), 1 + x1 + x2,
#'       ifelse(class_raw %in% c(3, 4), 2 + x1 + x2,
#'                                3 + x1 + x2))
#'
#' sds <- ifelse(class_raw %in% c(1, 2), 0.5,
#'       ifelse(class_raw %in% c(3, 4), 0.3,
#'                         0.4))
#'
#' y <- rlnorm(1000, meanlog = mu, sdlog = sds)
#'
#' df <- tibble(x1, x2, class = factor(class_raw), y)
#'
#' # Split into training, calibration, and test sets
#' df_train <- df %>% slice(1:500)
#' df_cal <- df %>% slice(501:750)
#' df_test <- df %>% slice(751:1000)
#'
#' # Fit model (on log-scale)
#' mod <- lm(log(y) ~ x1 + x2, data = df_train)
#'
#' # Generate predictions
#' pred_cal <- exp(predict(mod, newdata = df_cal))
#' pred_test <- exp(predict(mod, newdata = df_test))
#'
#' # Apply clustered conformal prediction
#' intervals <- pinterval_ccp(
#'   pred = pred_test,
#'   pred_class = df_test$class,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   calib_class = df_cal$class,
#'   alpha = 0.1,
#'   ncs_type = "absolute_error",
#'   optimize_n_clusters = TRUE,
#'   optimize_n_clusters_method = "calinhara",
#'   min_n_clusters = 2,
#'   max_n_clusters = 4
#' )
#'
#' # View clustered prediction intervals
#' head(intervals)
#'
pinterval_ccp = function(pred,
															pred_class = NULL,
															calib = NULL,
															calib_truth = NULL,
															calib_class = NULL,
															lower_bound = NULL,
															upper_bound = NULL,
															alpha = 0.1,
															ncs_type = c('absolute_error',
																					 'relative_error',
																					 'za_relative_error',
																					 'heterogeneous_error',
																					 'raw_error'),
												 n_clusters = NULL,
												 cluster_method = c('kmeans','ks'),
												 cluster_train_fraction = 1,
												 optimize_n_clusters = TRUE,
												 optimize_n_clusters_method = c('calinhara','min_cluster_size'),
												 min_cluster_size = 150,
												 min_n_clusters = 2,
												 max_n_clusters = NULL,
															distance_weighted_cp = FALSE,
															distance_features_calib = NULL,
															distance_features_pred = NULL,
															normalize_distance = TRUE,
															weight_function = c('gaussian_kernel', 'caucy_kernel','logistic','reciprocal_linear'),
															calibrate = FALSE,
															calibration_method = 'glm',
															calibration_family = NULL,
															calibration_transform = NULL,
															resolution = 0.01,
															grid_size = NULL){

	i <- NA

	min_class_size <- max(10, ceiling(1/alpha))


	if(setdiff(unique(pred_class), unique(calib_class)) %>% length() > 0){
		warning('Some classes in pred_class are not present in calib_class. These will result in NA prediction intervals for those classes.')
	}


	if(!is.numeric(pred) && ncol(pred) != 2){
		stop('pred must be a numeric scalar or vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the class labels')
	}

	if(is.numeric(pred) && is.null(pred_class)){
		stop('If pred is numeric, pred_class must be provided')
	}

	if(!is.numeric(pred)){
		pred_class <- as.numeric(pred[[2]])
		pred <- as.numeric(pred[[1]])
	}

	if(is.numeric(calib) & is.null(calib_truth)){
		stop('If calib is numeric, calib_truth must be provided')
	}

	if(is.numeric(calib) & is.null(calib_class)){
		stop('If calib is numeric, calib_class must be provided')
	}
	if(!is.numeric(calib) && ncol(calib)<3){
		stop('calib must be a numeric vector or a 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and the third column being the class labels')
	}

	if(!is.numeric(alpha) || alpha<=0 || alpha>=1 || length(alpha)!=1){
		stop('alpha must be a single numeric value between 0 and 1')
	}

	ncs_type <- match.arg(ncs_type, c('absolute_error',
																		'relative_error',
																		'za_relative_error',
																		'heterogeneous_error',
																		'raw_error'))



	if(!is.numeric(calib)){
		calib_org <- calib
		if(is.matrix(calib)){
			calib <- as.numeric(calib_org[,1])
			calib_truth <- as.numeric(calib_org[,2])
			if(is.null(calib_class) && ncol(calib_org) == 3){
				calib_class <- as.numeric(calib_org[,3])
			}
		}else{
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
			if(is.null(calib_class) && ncol(calib_org) == 3){
				calib_class <- as.numeric(calib_org[[3]])
			}
		}
	}

	if(ncs_type == 'heterogeneous_error'){
		coefs <- stats::coef(stats::lm(abs(calib - calib_truth) ~ calib))
	}else{
		coefs <- NULL
	}

	if(distance_weighted_cp){
		if(is.null(distance_features_calib) || is.null(distance_features_pred)){
			stop('If distance_weighted_cp is TRUE, distance_features_calib and distance_features_pred must be provided')
		}
		if(!is.matrix(distance_features_calib) && !is.data.frame(distance_features_calib) && !is.numeric(distance_features_calib)){
			stop('distance_features_calib must be a matrix, data frame, or numeric vector')
		}
		if(!is.matrix(distance_features_pred) && !is.data.frame(distance_features_pred) && !is.numeric(distance_features_pred)){
			stop('distance_features_pred must be a matrix, data frame, or numeric vector')
		}
		if(is.numeric(distance_features_calib) && is.numeric(distance_features_pred)){
			if(length(distance_features_calib) != length(calib) || length(distance_features_pred) != length(pred)){
				stop('If distance_features_calib and distance_features_pred are numeric vectors, they must have the same length as calib and pred, respectively')
			}
		}else if(is.matrix(distance_features_calib) || is.data.frame(distance_features_calib)){
			if(nrow(distance_features_calib) != length(calib)){
				stop('If distance_features_calib is a matrix or data frame, it must have the same number of rows as calib')
			}
			if(ncol(distance_features_calib) != ncol(distance_features_pred)){
				stop('distance_features_calib and distance_features_pred must have the same number of columns')
			}
			if(nrow(distance_features_pred) != length(pred)){
				stop('If distance_features_pred is a matrix or data frame, it must have the same number of rows as pred')
			}

		}

		distance_features_calib <- as.matrix(distance_features_calib)
		distance_features_pred <- as.matrix(distance_features_pred)

		if(!is.function(weight_function)){
			weight_function <- match.arg(weight_function, c('gaussian_kernel', 'caucy_kernel','logistic','reciprocal_linear'))
			weight_function <- switch(weight_function,
																'gaussian_kernel' = function(d) exp(-d^2),
																'caucy_kernel' = function(d) 1 / (1 + d^2),
																'logistic' = function(d) 1 / (1 + exp(d)),
																'reciprocal_linear' = function(d) 1 / (1 + d))
		}
	}

class_labels <- sort(unique(calib_class))
	if(length(class_labels)<2){
		stop('Calibration set must have at least two classes For continuous prediction intervals without classes, use pinterval_conformal() instead of pinterval_mondrian()')
	}


if(!optimize_n_clusters && is.null(n_clusters)){
		stop('If optimize_n_clusters is FALSE, n_clusters must be provided')
	}

optimize_n_clusters_method <- match.arg(optimize_n_clusters_method, c('calinhara','min_cluster_size'))
cluster_method <- match.arg(cluster_method, c('kmeans','ks'))

if(cluster_train_fraction == 1){
	warning('cluster_train_fraction is set to 1, which means the entire calibration set will be used for clustering. This may lead to overfitting.')
	calib_cluster <- calib
	calib_cluster_class <- calib_class
	calib_cluster_truth <- calib_truth
}else{
	if(cluster_train_fraction <= 0 || cluster_train_fraction >= 1){
		stop('cluster_train_fraction must be a numeric value between 0 and 1')
	}
	calib_cluster_ids <- sample(1:length(calib), size = floor(length(calib) * cluster_train_fraction), replace = FALSE)
	calib_cluster <- calib[calib_cluster_ids]
	calib_cluster_class <- calib_class[calib_cluster_ids]
	calib_cluster_truth <- calib_truth[calib_cluster_ids]
	calib <- calib[-calib_cluster_ids]
	calib_truth <- calib_truth[-calib_cluster_ids]
	calib_class <- calib_class[-calib_cluster_ids]
}

ncs_calib_cluster <- ncs_compute(ncs_type, calib_cluster, calib_cluster_truth, coefs)


if(optimize_n_clusters && optimize_n_clusters_method == 'min_cluster_size'){
	if(is.null(min_cluster_size) || !is.numeric(min_cluster_size) || length(min_cluster_size) != 1 || min_cluster_size <= 0){
		stop('If optimize_n_clusters_method is "min_cluster_size", min_cluster_size must be a single positive numeric value')
	}
	ntilde <- max(1/alpha - 1, min(table(calib_class)))
	ktilde <- sum(table(calib_class) >= ntilde)
	gamma <- ktilde/(ktilde + min_cluster_size/2)

	n_clusters <- floor(ntilde/(gamma*2))



	calib_cluster_vec <- clusterer(ncs_calib_cluster, n_clusters, calib_cluster_class,
																 method = cluster_method, min_class_size = min_class_size)


}else if(optimize_n_clusters && optimize_n_clusters_method == 'calinhara'){

	# if(is.null(n_clusters) && is.null(min_n_clusters) && is.null(max_n_clusters)){
	# 	stop('If optimize_n_clusters_method is "calinhara", min_n_clusters, and max_n_clusters must be provided, or n_clusters must be provided as a vector of n_clusters to optimize over')
	# }


	if(is.null(max_n_clusters)){
		max_n_clusters <- length(unique(calib_cluster_class)) - 1
		warning('max_n_clusters is not provided, setting to number of unique Mondrian classes minus 1: ', max_n_clusters)
	}

if(!is.null(n_clusters) && length(n_clusters == 1) && is.numeric(min_n_clusters) && is.numeric(max_n_clusters)){
	warning('Optimize clusters is set to TRUE, but n_clusters is provided as a single value. This will be ignored and the number of clusters will be optimized using the Calinhara method with the min_n_clusters and max_n_clusters parameters.')
}

if(!((is.numeric(min_n_clusters) && length(min_n_clusters) == 1 && min_n_clusters > 0) && (is.numeric(max_n_clusters) && length(max_n_clusters) == 1 && max_n_clusters > min_n_clusters)) && is.null(n_clusters)){
		stop('If optimize_n_clusters_method is "calinhara", min_n_clusters and max_n_clusters must be single positive numeric values with max_n_clusters > min_n_clusters, or n_clusters must be provided as a vector of n_clusters to optimize over')
}
	if(is.numeric(n_clusters) && length(n_clusters) > 1){
		if(!is.null(min_n_clusters) || !is.null(max_n_clusters)){
			warning('n_clusters is provided as a vector, so min_n_clusters and max_n_clusters will be ignored')


	}
		ms <- n_clusters
	}else{
		ms <- seq(from = min_n_clusters, to = max_n_clusters, by = 1)
	}




	calib_cluster_vec <- optimize_clusters(ncs_calib_cluster,
																				 calib_cluster_class,
																				 method = cluster_method,
																				 ms = ms)


}else if(!optimize_n_clusters){
	if(is.null(n_clusters) || !is.numeric(n_clusters) || length(n_clusters) != 1 || n_clusters <= 0){
		stop('If optimize_n_clusters is FALSE, n_clusters must be a single positive numeric value')
	}
	calib_cluster_vec <- clusterer(ncs_calib_cluster, n_clusters, calib_cluster_class,
																 method = cluster_method, min_class_size = min_class_size)

}


if(any(table(calib_cluster_class) < min_class_size)){
	warning('Some classes in in the cluster calibration set have less than ', min_class_size, ' calibration points. These classes are assigned to the NULL cluster and will utilize the full calibration set for prediction intervals, rather than cluster-specific intervals.')
}

calib_clusters <- class_to_clusters(calib_class, calib_cluster_vec)
pred_clusters <- class_to_clusters(pred_class, calib_cluster_vec)

cluster_labels <- sort(unique(pred_clusters))
if(any(is.na(pred_clusters))){
	cluster_labels <- c(cluster_labels, NA)
}


	cp_intervals <- foreach::foreach(i = 1:length(cluster_labels),.final=dplyr::bind_rows) %do%{
		indices <- which(pred_clusters == cluster_labels[i])
		if(is.na(cluster_labels[i])){
			res <- suppressWarnings(pinterval_conformal(pred = pred[is.na(pred_clusters)],
																									lower_bound = lower_bound,
																									upper_bound = upper_bound,
																									ncs_type = ncs_type,
																									calib = calib,
																									calib_truth = calib_truth,
																									calibrate = calibrate,
																									calibration_method = calibration_method,
																									calibration_family = calibration_family,
																									alpha = alpha,
																									resolution = resolution,
																									grid_size = grid_size))
			res$indices <- which(is.na(pred_clusters))
		}else{
		res <- suppressWarnings(pinterval_conformal(pred = na.omit(pred[pred_clusters==cluster_labels[i]]),
																								lower_bound = lower_bound,
																								upper_bound = upper_bound,
																								ncs_type = ncs_type,
																								calib = na.omit(calib[calib_clusters==cluster_labels[i]]),
																								calib_truth = na.omit(calib_truth[calib_clusters==cluster_labels[i]]),
																								calibrate = calibrate,
																								calibration_method = calibration_method,
																								calibration_family = calibration_family,
																								alpha = alpha,
																								resolution = resolution,
																								grid_size = grid_size))
		res$indices <- indices
		}
		res
		}


	cp_intervals2 <- cp_intervals %>%
		dplyr::arrange(indices) %>%
		dplyr::select(-indices) %>%
		dplyr::mutate(class = pred_class,
									cluster = pred_clusters)

	attr(cp_intervals2, 'clusters') = attr(calib_cluster_vec, 'clusters')
	attr(cp_intervals2, 'method') <- attr(calib_cluster_vec, 'method')
	attr(cp_intervals2, 'n_clusters') <- attr(calib_cluster_vec, 'm')
	attr(cp_intervals2, 'calib_ch_index') <- attr(calib_cluster_vec, 'ch_index')
	attr(cp_intervals2, 'calib_coverage_gaps') <- attr(calib_cluster_vec, 'coverage_gaps')

	if(calibrate){
		calibrated_preds <- calibrate_predictions(pred = pred,
																							calib = calib,
																							calib_truth = calib_truth,
																							method = calibration_method,
																							family = calibration_family,
																							transform = calibration_transform)
		cp_intervals2 <- cp_intervals2 %>%
			dplyr::mutate(pred = calibrated_preds)
	}

	return(cp_intervals2)
}
