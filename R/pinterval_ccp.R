#' Clustered Conformal Prediction Intervals for Continuous Predictions
#'
#' @description
#' This function computes conformal prediction intervals with a confidence level of \eqn{1 - \alpha} by first grouping Mondrian classes into data-driven clusters based on the distribution of their nonconformity scores. The resulting clusters are used as strata for computing class-conditional (Mondrian-style) conformal prediction intervals. This approach improves local validity and statistical efficiency when there are many small or similar classes with overlapping prediction behavior. The coverage level \eqn{1 - \alpha} is approximate within each cluster, assuming exchangeability of nonconformity scores within clusters.
#'
#' The method supports additional features such as prediction calibration, distance-weighted conformal scores, and clustering optimization via internal validity measures (e.g., Calinski-Harabasz index or minimum cluster size heuristics).
#'
#'
#' @inheritParams pinterval_conformal
#' @inheritParams pinterval_mondrian
#' @param n_clusters Number of clusters to use when combining Mondrian classes. Required if \code{optimize_n_clusters = FALSE}.
#' @param cluster_method Clustering method used to group Mondrian classes. Options are \code{"kmeans"} or \code{"ks"} (Kolmogorov-Smirnov). Default is \code{"kmeans"}.
#' @param cluster_train_fraction Fraction of the calibration data used to estimate nonconformity scores and compute clustering. Default is 1, which uses the entire calibration set for both clustering and interval estimation. See details for more discussion.
#' @param optimize_n_clusters Logical. If \code{TRUE}, the number of clusters is chosen automatically based on internal clustering criteria.
#' @param optimize_n_clusters_method Method used for cluster optimization. One of \code{"calinhara"} (Calinski-Harabasz index) or \code{"min_cluster_size"}. Default is \code{"calinhara"}.
#' @param min_cluster_size Minimum number of calibration points per cluster. Used only when \code{optimize_n_clusters_method = "min_cluster_size"}.
#' @param min_n_clusters Minimum number of clusters to consider when optimizing.
#' @param max_n_clusters Maximum number of clusters to consider. If \code{NULL}, the upper limit is set to the number of unique Mondrian classes minus 1.
#'
#' @return A tibble with predicted values, lower and upper prediction interval bounds, class labels, and assigned cluster labels. Attributes include clustering diagnostics (e.g., cluster assignments, coverage gaps, internal validity scores).
#'
#' @details
#' `pinterval_ccp()` builds on [pinterval_mondrian()] by introducing a
#' clustered conformal prediction framework. Instead of requiring a separate calibration distribution for every Mondrian class, which may lead to unstable or noisy intervals when there are many small groups, the method groups similar Mondrian classes into clusters with similar nonconformity score distributions. Classes with similar prediction-error behavior are assigned to the same cluster. Each resulting cluster is then treated as a stratum for standard inductive conformal prediction.
#'
#' Users may specify the number of clusters directly using the `n_clusters` argument or optimize the number of clusters using the Calinski–Harabasz index or minimum cluster size heuristics.
#'
#' Clustering can be computed using all calibration data or a subsample defined by `cluster_train_fraction`. By default, the entire calibration set is used for both clustering and interval estimation, which may lead to overfitting. Setting `cluster_train_fraction` to a value less than 1 (e.g., 0.5) can help mitigate this risk by using separate data for clustering and interval estimation, at the cost of potentially less stable cluster assignments with smaller calibration subsets. If data is limited, using the full calibration set for clustering may still be preferable, but users should be aware of the potential for overfitting and optimistic coverage estimates in this case.
#'
#' Clustering is based on either k-means or Kolmogorov-Smirnov distance between nonconformity score distributions of the Mondrian classes, selected via the `cluster_method` argument.
#'
#' For a detailed description of non-conformity scores, distance-weighting, and the general conformal  prediction framework, see [pinterval_conformal()], and for a description of Mondrian conformal prediction, see [pinterval_mondrian()].
#'
#' @seealso \code{\link[pintervals]{pinterval_conformal}}, \code{\link[pintervals]{pinterval_mondrian}}
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
pinterval_ccp = function(
	pred,
	pred_class = NULL,
	calib = NULL,
	calib_truth = NULL,
	calib_class = NULL,
	lower_bound = NULL,
	upper_bound = NULL,
	alpha = 0.1,
	ncs_type = c(
		'absolute_error',
		'relative_error',
		'za_relative_error',
		'heterogeneous_error',
		'raw_error'
	),
	grid_size = 10000,
	resolution = NULL,
	n_clusters = NULL,
	cluster_method = c('kmeans', 'ks'),
	cluster_train_fraction = 0.5,
	optimize_n_clusters = TRUE,
	optimize_n_clusters_method = c('calinhara', 'min_cluster_size'),
	min_cluster_size = 150,
	min_n_clusters = 2,
	max_n_clusters = NULL,
	distance_weighted_cp = FALSE,
	distance_features_calib = NULL,
	distance_features_pred = NULL,
	distance_type = c('mahalanobis', 'euclidean'),
	normalize_distance = 'none',
	weight_function = c(
		'gaussian_kernel',
		'caucy_kernel',
		'logistic',
		'reciprocal_linear'
	)
) {
	i <- NA

	min_class_size <- max(10, ceiling(1 / alpha))

	if (any(is.na(n_clusters)) || any(n_clusters <= 0)) {
		stop(
			'pinterval_ccp: n_clusters be a single positive numeric value or a vector of positive numeric values to optimize over',
			call. = FALSE
		)
	}

	# Validate pred
	if (!is.numeric(pred) && (!is.matrix(pred) && !is.data.frame(pred))) {
		stop(
			'pinterval_ccp: pred must be a numeric scalar or vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the class labels',
			call. = FALSE
		)
	}

	if (!is.numeric(pred) && ncol(pred) != 2) {
		stop(
			'pinterval_ccp: pred must be a numeric scalar or vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the class labels',
			call. = FALSE
		)
	}

	if (is.numeric(pred) && is.null(pred_class)) {
		stop(
			'pinterval_ccp: If pred is numeric, pred_class must be provided',
			call. = FALSE
		)
	}

	if (!is.numeric(pred)) {
		pred_class <- as.numeric(pred[[2]])
		pred <- as.numeric(pred[[1]])
	}

	# Check for NAs in pred
	if (any(is.na(pred))) {
		warning('pinterval_ccp: pred contains NA values', call. = FALSE)
	}

	# Validate pred_class length
	if (length(pred_class) != length(pred)) {
		stop(
			'pinterval_ccp: pred_class must have the same length as pred',
			call. = FALSE
		)
	}

	# Validate calib - NULL check and type check
	if (is.null(calib)) {
		stop('pinterval_ccp: calib must be provided', call. = FALSE)
	}

	if (!is.numeric(calib) && !is.matrix(calib) && !is.data.frame(calib)) {
		stop(
			'pinterval_ccp: calib must be a numeric vector or a 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and the third column being the class labels',
			call. = FALSE
		)
	}

	# Check calib column count if not numeric
	if (!is.vector(calib) && ncol(calib) < 3) {
		stop(
			'pinterval_ccp: calib must be a numeric vector or a 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and the third column being the class labels',
			call. = FALSE
		)
	}

	if (is.numeric(calib) && is.null(calib_truth)) {
		stop(
			'pinterval_ccp: If calib is a vector, calib_truth must be provided',
			call. = FALSE
		)
	}

	if (is.vector(calib) && is.vector(calib_class)) {
		stop(
			'pinterval_ccp: If calib is a vector, calib_class must be provided',
			call. = FALSE
		)
	}

	if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha) != 1) {
		stop(
			'pinterval_ccp: alpha must be a single numeric value between 0 and 1',
			call. = FALSE
		)
	}

	ncs_type <- match.arg(
		ncs_type,
		c(
			'absolute_error',
			'relative_error',
			'za_relative_error',
			'heterogeneous_error',
			'raw_error'
		)
	)

	# Parse calib into components
	if (!is.vector(calib)) {
		calib_org <- calib
		if (is.matrix(calib)) {
			calib <- as.numeric(calib_org[, 1])
			calib_truth <- as.numeric(calib_org[, 2])
			if (is.null(calib_class) && ncol(calib_org) == 3) {
				calib_class <- as.numeric(calib_org[, 3])
			}
		} else {
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
			if (is.null(calib_class) && ncol(calib_org) == 3) {
				calib_class <- as.numeric(calib_org[[3]])
			}
		}
	}

	# Check for NAs in calib and calib_truth
	if (any(is.na(calib))) {
		warning('pinterval_ccp: calib contains NA values', call. = FALSE)
	}
	if (any(is.na(calib_truth))) {
		warning('pinterval_ccp: calib_truth contains NA values', call. = FALSE)
	}

	# Validate calib_class and calib_truth length
	if (length(calib_class) != length(calib)) {
		stop(
			'pinterval_ccp: calib_class must have the same length as calib',
			call. = FALSE
		)
	}
	if (length(calib_truth) != length(calib)) {
		stop(
			'pinterval_ccp: calib_truth must have the same length as calib',
			call. = FALSE
		)
	}

	# Check setdiff warning AFTER calib_class and pred_class have been validated/parsed
	if (setdiff(unique(pred_class), unique(calib_class)) %>% length() > 0) {
		warning(
			'pinterval_ccp: Some classes in pred_class are not present in calib_class. These will result in NA prediction intervals for those classes.',
			call. = FALSE
		)
	}

	if (ncs_type == 'heterogeneous_error') {
		coefs <- stats::coef(stats::lm(abs(calib - calib_truth) ~ calib))
	} else {
		coefs <- NULL
	}

	# Validate normalize_distance
	normalize_distance <- match.arg(normalize_distance, c('none', 'minmax', 'sd'))

	if (distance_weighted_cp) {
		validate_distance_inputs(
			distance_features_calib,
			distance_features_pred,
			length(calib),
			length(pred),
			fn_name = "pinterval_ccp"
		)
		distance_features_calib <- as.matrix(distance_features_calib)
		distance_features_pred <- as.matrix(distance_features_pred)
		distance_type <- match.arg(distance_type, c('mahalanobis', 'euclidean'))
		weight_function <- resolve_weight_function(weight_function)
	}

	class_labels <- sort(unique(calib_class))
	if (length(class_labels) < 2) {
		stop(
			"pinterval_ccp: calibration set must have at least two classes. For continuous prediction intervals without classes, use pinterval_conformal() instead.",
			call. = FALSE
		)
	}

	if (!optimize_n_clusters && is.null(n_clusters)) {
		stop(
			'pinterval_ccp: If optimize_n_clusters is FALSE, n_clusters must be provided',
			call. = FALSE
		)
	}

	optimize_n_clusters_method <- match.arg(
		optimize_n_clusters_method,
		c('calinhara', 'min_cluster_size')
	)
	cluster_method <- match.arg(cluster_method, c('kmeans', 'ks'))

	if (
		!is.numeric(cluster_train_fraction) ||
			length(cluster_train_fraction) != 1 ||
			cluster_train_fraction <= 0 ||
			cluster_train_fraction > 1
	) {
		stop(
			"pinterval_ccp: 'cluster_train_fraction' must be a single numeric value in (0, 1].",
			call. = FALSE
		)
	}

	if (cluster_train_fraction == 1) {
		warning(
			"pinterval_ccp: 'cluster_train_fraction' is set to 1, meaning the entire calibration set is used for both clustering and interval estimation. This may lead to overfitting. Consider setting a value < 1 (e.g., 0.5).",
			call. = FALSE
		)
		calib_cluster <- calib
		calib_cluster_class <- calib_class
		calib_cluster_truth <- calib_truth
	} else {
		calib_cluster_ids <- sample(
			1:length(calib),
			size = floor(length(calib) * cluster_train_fraction),
			replace = FALSE
		)
		calib_cluster <- calib[calib_cluster_ids]
		calib_cluster_class <- calib_class[calib_cluster_ids]
		calib_cluster_truth <- calib_truth[calib_cluster_ids]
		calib <- calib[-calib_cluster_ids]
		calib_truth <- calib_truth[-calib_cluster_ids]
		calib_class <- calib_class[-calib_cluster_ids]
	}

	ncs_calib_cluster <- ncs_compute(
		ncs_type,
		calib_cluster,
		calib_cluster_truth,
		coefs
	)

	if (optimize_n_clusters && optimize_n_clusters_method == 'min_cluster_size') {
		if (
			is.null(min_cluster_size) ||
				!is.numeric(min_cluster_size) ||
				length(min_cluster_size) != 1 ||
				min_cluster_size <= 0
		) {
			stop(
				'pinterval_ccp: If optimize_n_clusters_method is "min_cluster_size", min_cluster_size must be a single positive numeric value',
				call. = FALSE
			)
		}
		ntilde <- max(1 / alpha - 1, min(table(calib_class)))
		ktilde <- sum(table(calib_class) >= ntilde)
		gamma <- ktilde / (ktilde + min_cluster_size / 2)

		n_clusters <- floor(ntilde / (gamma * 2))

		calib_cluster_vec <- clusterer(
			ncs_calib_cluster,
			n_clusters,
			calib_cluster_class,
			method = cluster_method,
			min_class_size = min_class_size
		)
	} else if (optimize_n_clusters && optimize_n_clusters_method == 'calinhara') {
		# if(is.null(n_clusters) && is.null(min_n_clusters) && is.null(max_n_clusters)){
		# 	stop('pinterval_ccp: If optimize_n_clusters_method is "calinhara", min_n_clusters, and max_n_clusters must be provided, or n_clusters must be provided as a vector of n_clusters to optimize over', call. = FALSE)
		# }

		if (is.null(max_n_clusters)) {
			max_n_clusters <- length(unique(calib_cluster_class)) - 1
			warning(
				'pinterval_ccp: max_n_clusters is not provided, setting to number of unique Mondrian classes minus 1: ',
				max_n_clusters,
				call. = FALSE
			)
		}

		if (
			!is.null(n_clusters) &&
				length(n_clusters) == 1 &&
				is.numeric(min_n_clusters) &&
				is.numeric(max_n_clusters)
		) {
			warning(
				'pinterval_ccp: Optimize clusters is set to TRUE, but n_clusters is provided as a single value. This will be ignored and the number of clusters will be optimized using the Calinhara method with the min_n_clusters and max_n_clusters parameters.',
				call. = FALSE
			)
		}

		if (
			!((is.numeric(min_n_clusters) &&
				length(min_n_clusters) == 1 &&
				min_n_clusters > 0) &&
				(is.numeric(max_n_clusters) &&
					length(max_n_clusters) == 1 &&
					max_n_clusters > min_n_clusters)) &&
				is.null(n_clusters)
		) {
			stop(
				'pinterval_ccp: If optimize_n_clusters_method is "calinhara", min_n_clusters and max_n_clusters must be single positive numeric values with max_n_clusters > min_n_clusters, or n_clusters must be provided as a vector of n_clusters to optimize over',
				call. = FALSE
			)
		}
		if (is.numeric(n_clusters) && length(n_clusters) > 1) {
			if (!is.null(min_n_clusters) || !is.null(max_n_clusters)) {
				warning(
					'pinterval_ccp: n_clusters is provided as a vector, so min_n_clusters and max_n_clusters will be ignored',
					call. = FALSE
				)
			}
			ms <- n_clusters
		} else {
			ms <- seq(from = min_n_clusters, to = max_n_clusters, by = 1)
		}

		calib_cluster_vec <- optimize_clusters(
			ncs_calib_cluster,
			calib_cluster_class,
			method = cluster_method,
			ms = ms
		)
	} else if (!optimize_n_clusters) {
		if (
			is.null(n_clusters) ||
				!is.numeric(n_clusters) ||
				length(n_clusters) != 1 ||
				n_clusters <= 0
		) {
			stop(
				'pinterval_ccp: If optimize_n_clusters is FALSE, n_clusters must be a single positive numeric value',
				call. = FALSE
			)
		}

		# Warn if n_clusters is not an integer
		if (!is.integer(n_clusters) && (n_clusters %% 1 != 0)) {
			warning(
				'pinterval_ccp: n_clusters should be an integer value',
				call. = FALSE
			)
		}

		calib_cluster_vec <- clusterer(
			ncs_calib_cluster,
			n_clusters,
			calib_cluster_class,
			method = cluster_method,
			min_class_size = min_class_size
		)
	}

	if (any(table(calib_cluster_class) < min_class_size)) {
		warning(
			'pinterval_ccp: Some classes in in the cluster calibration set have less than ',
			min_class_size,
			' calibration points. These classes are assigned to the NULL cluster and will utilize the full calibration set for prediction intervals, rather than cluster-specific intervals.',
			call. = FALSE
		)
	}

	calib_clusters <- class_to_clusters(calib_class, calib_cluster_vec)
	pred_clusters <- class_to_clusters(pred_class, calib_cluster_vec)

	cluster_labels <- sort(unique(pred_clusters))
	if (any(is.na(pred_clusters))) {
		cluster_labels <- c(cluster_labels, NA)
	}

	cp_intervals <- foreach::foreach(
		i = 1:length(cluster_labels),
		.final = dplyr::bind_rows
	) %do%
		{
			indices <- which(pred_clusters == cluster_labels[i])
			if (is.na(cluster_labels[i])) {
				res <- suppressWarnings(pinterval_conformal(
					pred = pred[is.na(pred_clusters)],
					lower_bound = lower_bound,
					upper_bound = upper_bound,
					ncs_type = ncs_type,
					calib = calib,
					calib_truth = calib_truth,
					distance_weighted_cp = distance_weighted_cp,
					distance_features_calib = distance_features_calib,
					distance_features_pred = distance_features_pred,
					distance_type = distance_type,
					normalize_distance = normalize_distance,
					weight_function = weight_function,
					alpha = alpha,
					resolution = resolution,
					grid_size = grid_size
				))
				res$indices <- which(is.na(pred_clusters))
			} else {
				res <- suppressWarnings(pinterval_conformal(
					pred = na.omit(pred[pred_clusters == cluster_labels[i]]),
					lower_bound = lower_bound,
					upper_bound = upper_bound,
					ncs_type = ncs_type,
					calib = na.omit(calib[calib_clusters == cluster_labels[i]]),
					calib_truth = na.omit(calib_truth[
						calib_clusters == cluster_labels[i]
					]),
					distance_weighted_cp = distance_weighted_cp,
					distance_features_calib = distance_features_calib[
						calib_clusters == cluster_labels[i],
					],
					distance_features_pred = distance_features_pred[
						pred_clusters == cluster_labels[i],
					],
					distance_type = distance_type,
					normalize_distance = normalize_distance,
					weight_function = weight_function,
					alpha = alpha,
					resolution = resolution,
					grid_size = grid_size
				))
				res$indices <- indices
			}
			res
		}

	cp_intervals2 <- cp_intervals %>%
		dplyr::arrange(indices) %>%
		dplyr::select(-indices) %>%
		dplyr::mutate(class = pred_class, cluster = pred_clusters)

	attr(cp_intervals2, 'clusters') = attr(calib_cluster_vec, 'clusters')
	attr(cp_intervals2, 'method') <- attr(calib_cluster_vec, 'method')
	attr(cp_intervals2, 'n_clusters') <- attr(calib_cluster_vec, 'm')
	attr(cp_intervals2, 'calib_ch_index') <- attr(calib_cluster_vec, 'ch_index')
	attr(cp_intervals2, 'calib_coverage_gaps') <- attr(
		calib_cluster_vec,
		'coverage_gaps'
	)

	return(cp_intervals2)
}
