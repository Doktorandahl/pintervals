
#' Non-Conformity Score Computation Function
#' @param type Type of non-conformity score to compute. Options include 'absolute_error', 'raw_error', 'relative_error', 'relative_error2', and 'heterogeneous_error'.
#' @param pred a numeric vector of predicted values
#' @param truth a numeric vector of true values
#' @param coefs a numeric vector of coefficients for the heterogeneous error model. Must be of length 2, where the first element is the intercept and the second element is the slope.
ncs_compute <- function(type, pred, truth, coefs = NULL){
	if(type == 'absolute_error'){
		return(abs_error(pred, truth))
	}else if(type == 'raw_error'){
		return(raw_error(pred, truth))
	}else if(type == 'relative_error'){
		return(rel_error(pred, truth))
	}else if(type == 'za_relative_error'){
			return(za_rel_error(pred, truth))
	}else if(type == 'heterogeneous_error'){
		if(is.null(coefs)){
			stop('coefs must be provided for heterogeneous error')
		}
		return(heterogeneous_error(pred, truth, coefs))
	}else{
		stop('Unknown non-conformity score type')
	}
}

#' Gaussian Kernel Function
#' @param d a numeric vector of distances
#' @return a numeric vector of Gaussian kernel values
gauss_kern <- function(d){
	return(exp(-d^2/2))
}

#' Cauchy Kernel Function
#' @param d a numeric vector of distances
#' @return a numeric vector of Cauchy kernel values
cauchy_kern <- function(d){
	return(1/(1 + d^2))
}

#' Logistic Kernel Function
#' @param d a numeric vector of distances
#' @return a numeric vector of logistic kernel values
logistic_kern <- function(d){
	return(1/(1 + exp(d)))
}

#' Reciprocal Linear Kernel Function
#' @param d a numeric vector of distances
#' @return a numeric vector of reciprocal linear kernel values
reciprocal_linear_kern <- function(d){
	return(1/(1 + d))
}



#' Absolute Error Function for Non-Conformity Scores
#'
#' @param pred a numeric vector of predicted values
#' @param truth a numeric vector of true values
#'
#' @return a numeric vector of absolute errors
#'
abs_error <- function(pred, truth){
	return(abs(pred-truth))
}

#' Raw Error Function for Non-Conformity Scores
#'
#' @param pred a numeric vector of predicted values
#' @param truth a numeric vector of true values
#' @return a numeric vector of raw errors
raw_error <- function(pred, truth){
	return(truth-pred)
}

#' Relative Error Function for Non-Conformity Scores by predicted Values
#'
#' @param pred a numeric vector of predicted values
#' @param truth a numeric vector of true values
#'
#' @return a numeric vector of relative errors
#'
rel_error <- function(pred, truth){
	return(abs((pred-truth)/pred))
}

#' Zero-adjusted Relative Error Function for Non-Conformity Scores by predicted Values with a small adjustment
#' @param pred a numeric vector of predicted values
#' @param truth a numeric vector of true values
#' @return a numeric vector of zero-adjusted relative errors
za_rel_error <- function(pred, truth){
		return(abs((pred-truth)/(1+pred)))
}


#' Heterogeneous Error Function for Non-Conformity Scores
#' @param pred a numeric vector of predicted values
#' @param truth a numeric vector of true values
#' @param coefs a numeric vector of coefficients for the heterogeneous error model. Must be of length 2, where the first element is the intercept and the second element is the slope.
heterogeneous_error <- function(pred, truth, coefs){
	if(length(coefs) != 2){
		stop('coefs must be a vector of length 2')
	}

	est_heterogeneous_error <- coefs[1] + coefs[2] * pred
	return(abs(pred - truth) / est_heterogeneous_error)

}

#' Grid search for lower and upper bounds of continuous conformal prediction intervals
#'
#' @param y_min minimum value to search
#' @param y_max maximum value to search
#' @param ncs vector of non-conformity scores
#' @param y_hat vector of predicted values
#' @param alpha confidence level
#' @param min_step The minimum step size for the grid search
#' @param grid_size  Alternative to min_step, the number of points to use in the grid search between the lower and upper bound
#' @param ncs_type String indicating the non-conformity score function to use
#' @param calib a tibble with the predicted values and the true values of the calibration partition. Used when weighted_cp is TRUE. Default is NULL
#' @param coefs a numeric vector of coefficients for the heterogeneous error model. Must be of length 2, where the first element is the intercept and the second element is the slope. Used when ncs_type is 'heterogeneous_error'. Default is NULL
#' @param distance_weighted_cp logical. If TRUE, the non-conformity scores will be weighted according to the distance function
#' @param distance_features_calib a matrix of features for the calibration partition. Used when distance_weighted_cp is TRUE
#' @param distance_features_pred a matrix of features for the prediction partition. Used when distance_weighted_cp is TRUE
#' @param normalize_distance logical. If TRUE, the distances will be normalized to the range [0,1]
#' @param weight_function a function to use for weighting the distances. Can be 'gaussian_kernel', 'caucy_kernel', 'logistic', or 'reciprocal_linear'. Default is 'gaussian_kernel'
#'
#' @return a tibble with the predicted values and the lower and upper bounds of the prediction intervals
#'
grid_finder <- function(y_min,y_max,ncs,ncs_type,y_hat, alpha, min_step = NULL, grid_size = NULL, calib = NULL, coefs = NULL,distance_weighted_cp = FALSE,
												distance_features_calib = NULL,
												distance_features_pred = NULL,
												normalize_distance = TRUE,
												weight_function = gauss_kern){

	i <- NA
	if(is.null(grid_size)){
		pos_vals <- seq(from=y_min,to=y_max,by=min_step)
		if(length(pos_vals)>10000){
			warning('Grid size with set step size is large, consider adjusting min_step or using grid_size instead of min_step if the search is too slow')
		}
	}else{
		pos_vals <- seq(from=y_min,to=y_max,length.out=grid_size)
	}




		out <- foreach::foreach(i = 1:length(y_hat)) %do%
			grid_inner(ncs_compute(ncs_type,y_hat[i],pos_vals, coefs = coefs),y_hat[i],ncs,pos_vals,alpha,ncs_type,distance_weighted_cp,
								 distance_features_calib,
								 distance_features_pred[i,],
								 normalize_distance,
								 weight_function)


	return(dplyr::bind_rows(out))
}

#' Inner function for grid search
#'
#' @param hyp_ncs vector of hypothetical non-conformity scores
#' @param y_hat predicted value
#' @param ncs vector of non-conformity scores
#' @param pos_vals vector of possible values for the lower and upper bounds of the prediction interval
#' @param alpha confidence level
#' @param ncs_type type of non-conformity score
#' @param distance_weighted_cp logical. If TRUE, the non-conformity scores will be weighted according to the distance function
#' @param distance_features_calib a matrix of features for the calibration partition. Used when distance_weighted_cp is TRUE
#' @param distance_features_pred a matrix of features for the prediction partition. Used when distance_weighted_cp is TRUE
#' @param normalize_distance logical. If TRUE, the distances will be normalized to the range [0,1]
#' @param weight_function a function to use for weighting the distances. Can be 'gaussian_kernel', 'caucy_kernel', 'logistic', or 'reciprocal_linear'. Default is 'gaussian_kernel'
#'
#' @return a numeric vector with the predicted value and the lower and upper bounds of the prediction interval
grid_inner <- function(hyp_ncs,y_hat,ncs,pos_vals,alpha,ncs_type,
											 distance_weighted_cp,
											 distance_features_calib,
											 distance_features_pred,
											 normalize_distance,
											 weight_function){
i <- NULL
if(distance_weighted_cp){
	distances <- foreach::foreach(i = 1:nrow(distance_features_calib)) %do%
		stats::dist(rbind(distance_features_calib[i,]),distance_features_pred)[1]

	if(normalize_distance){
		distances <- distances/max(distances)
	}

	wt <- weight_function(distances)


}else{
	wt <- rep(1, length(ncs))
}


	if(ncs_type != 'raw_error'){

	if(sum(hyp_ncs<Hmisc::wtd.quantile(ncs, weights = wt, normwt=TRUE, probs = 1-alpha))==0){
			return(c(pred = as.numeric(y_hat), lower_bound = NA_real_, upper_bound = NA_real_))

	}else{
		lb <- min(pos_vals[hyp_ncs<=Hmisc::wtd.quantile(ncs, weights = wt, normwt=TRUE, probs =1-alpha)])
		ub <- max(pos_vals[hyp_ncs<=Hmisc::wtd.quantile(ncs, weights = wt, normwt=TRUE, probs =1-alpha)])

		return(c(pred = as.numeric(y_hat), lower_bound = lb, upper_bound = ub))
		}
	}else{
		if(sum(hyp_ncs<Hmisc::wtd.quantile(ncs,weights = wt, normwt=TRUE, probs =1-alpha/2) &
					 hyp_ncs>Hmisc::wtd.quantile(ncs,weights = wt, normwt=TRUE, probs =alpha/2))==0){
			return(c(pred = as.numeric(y_hat), lower_bound = NA_real_, upper_bound = NA_real_))
		}else{
			lb <- min(pos_vals[hyp_ncs>=Hmisc::wtd.quantile(ncs,weights = wt, normwt=TRUE, probs =alpha/2)])
			ub <- max(pos_vals[hyp_ncs<=Hmisc::wtd.quantile(ncs,weights = wt, normwt=TRUE, probs =1-alpha/2)])

			return(c(pred = as.numeric(y_hat), lower_bound = lb, upper_bound = ub))
		}
	}



	}

#' Bootstrap function for bootstrapping the prediction intervals
#'
#' @param pred predicted value
#' @param error vector of errors
#' @param nboot number of bootstrap samples
#' @param alpha confidence level
#' @param error vector of errors.
#' @param dw_bootstrap logical. If TRUE, the bootstrap samples will be weighted according to the distance function
#' @param calib a vector of true values of the calibration partition. Used when weighted_bootstrap is TRUE
#' @param error_type The type of error to use for the prediction intervals. Can be 'raw' or 'absolute'. If 'raw', bootstrapping will be done on the raw prediction errors. If 'absolute', bootstrapping will be done on the absolute prediction errors with random signs. Default is 'raw'
#' @param distance_weighted_bootstrap logical. If TRUE, the bootstrap samples will be weighted according to the distance function
#' @param distance_function a function that takes two numeric vectors and returns a numeric vector of distances. Default is NULL, in which case the absolute error will be used
#' @param distance_features_calib a matrix of features for the calibration partition. Used when distance_weighted_bootstrap is TRUE
#' @param distance_features_pred a matrix of features for the prediction partition. Used when distance_weighted_bootstrap is TRUE
#' @param normalize_distance logical. If TRUE, the distances will be normalized to the range [0,1]
#' @param weight_function a function to use for weighting the distances. Can be 'gaussian_kernel', 'caucy_kernel', 'logistic', or 'reciprocal_linear'. Default is 'gaussian_kernel'
#'
#' @return a numeric vector with the predicted value and the lower and upper bounds of the prediction interval
bootstrap_inner <- function(pred, calib, error, nboot, alpha, dw_bootstrap = FALSE, distance_function = NULL, error_type = c('raw','absolute'),
														distance_weighted_bootstrap = FALSE,
														distance_features_calib = NULL,
														distance_features_pred = NULL,
														normalize_distance = TRUE,
														weight_function = gauss_kern){
	i <- NA
if(!distance_weighted_bootstrap){
	boot_error <- sample(error, size = nboot, replace = TRUE)
}else{

	distances <- foreach::foreach(i = 1:nrow(distance_features_calib)) %do%
		stats::dist(rbind(distance_features_calib[i,]),distance_features_pred)[1]

	if(normalize_distance){
		distances <- distances/max(distances)
	}
	boot_error <- sample(error, size = nboot, replace = TRUE, prob = weight_function(distances))
}

	if(error_type == 'absolute'){
		boot_error <- sample(c(-1,1), size = nboot, replace = TRUE) * boot_error
	}

	boot_pred <- pred + boot_error
	lb <- as.numeric(stats::quantile(boot_pred, alpha/2))
	ub <- as.numeric(stats::quantile(boot_pred, 1-alpha/2))


	return(c(pred = as.numeric(pred), lower_bound = lb, upper_bound = ub))
}


#' Bin chopper function for binned bootstrapping
#'
#' @param x vector of values to be binned
#' @param nbins number of bins
#' @param return_breaks logical indicating whether to return the bin breaks
bin_chopper <- function(x, nbins, return_breaks = FALSE){
	if(nbins < 2){
		stop('nbins must be greater than 1')
	}
	if(nbins > length(x)){
		stop('nbins must be less than or equal to the length of x')
	}
	if(length(unique(x)) == 1){
		stop('x must have more than one unique value')
	}
	if(length(unique(x)) < nbins){
		stop('x must have more unique values than nbins')
	}

	target_num <- ceiling(length(x)/nbins)

	qtiles <- seq(from = 0, to = 1, length.out = nbins+1)
	qtiles <- qtiles[-c(1,length(qtiles))]
	cutpoints_qtiles <- as.numeric(stats::quantile(x,qtiles))
	init_cut <- cut(x, breaks = c(-Inf,cutpoints_qtiles,Inf), labels = F)
	if(max(table(init_cut)) <= target_num + 2 & min(table(init_cut)) >= target_num - 2){
		cutpoints <- cutpoints_qtiles
	}else{
	nobs_per_value <- table(x)
	binsizes <- rep(target_num, nbins)
	binsizes2 <- rep(0, nbins)
	cutpoints <- rep(0, nbins-1)
	k <- 0
	while(!identical(binsizes,binsizes2) & k<10){
	for(i in 1:(nbins-1)){
		ccs <- 0
		j <- 0
		while(ccs < sum(binsizes[1:i])){
			j <- j + 1
			ccs <- sum(nobs_per_value[1:j])
		}
		if(i>1){
			ccs <- ccs - sum(binsizes[1:(i-1)])
		}
		binsizes[i] <- ccs
		cutpoints[i] <- as.numeric(names(nobs_per_value)[j])
		if(ccs > target_num & i<nbins){
			binsizes[(i+1):nbins] <- (length(x) - sum(binsizes[1:i]))/(nbins-i)
			target_num <- (length(x) - sum(binsizes[1:i]))/(nbins-i)
		}
	}
	binsizes[length(binsizes)] <- length(x) - sum(binsizes[1:(nbins-1)])
	binsizes2 <- binsizes
	k <- k + 1
	}
	}

if(!return_breaks){
	return(cut(x, breaks = c(-Inf,cutpoints, Inf), labels = FALSE))
}else{
	return(c(-Inf,cutpoints, Inf))
}
	}


#' Bin-individual alpha function for conformal prediction
#'
#' @param minqs Minimum quantiles
#' @param alpha alpha level
bindividual_alpha <- function(minqs, alpha){
	a <- alpha
	rem_bins <- sum(minqs >= a, na.rm=T)
	minqs[which(minqs < a)] <- NA

	if(all(is.na(minqs))){
		return(list(power = 0, bins = !is.na(minqs)))
	}

	a_tot <- prod(minq_to_alpha(minqs, a),na.rm=T)
	rem_bins_old <- rem_bins + 1

	while(rem_bins != rem_bins_old){

		if(prod(minq_to_alpha(minqs[-which.min(minqs)],a),na.rm=T)<=alpha){
			minqs[which.min(minqs)] <- NA
			rem_bins <- rem_bins - 1
		}

	if(min(minqs,na.rm=T) > a^(1/rem_bins)){
		minqs[which(minqs < a^(1/rem_bins))] <- a^(1/rem_bins)
		return(list(power = rem_bins, bins = !is.na(minqs)))
	}



		rem_bins_old <- rem_bins
	}

return(list(power = rem_bins, bins = !is.na(minqs)))

}

#' Helper for minimum quantile to alpha function
#'
#' @param minq minimum quantile
#' @param alpha alpha level
minq_to_alpha <- function(minq, alpha){
	minq[which(minq > alpha)] <- alpha
	return(minq)
}

#' Flatten binned conformal prediction intervals to contiguous intervals
#'
#' @param lst list of binned conformal prediction intervals
#' @param contiguize logical indicating whether to contiguize the intervals
flatten_cp_bin_intervals <- function(lst,
																		 contiguize = FALSE){

	i <- 'tmp'

	pred <- lst[[1]]$pred
	lower_bound <- foreach::foreach(i = 1:length(lst), .final = unlist) %do% lst[[i]]$lower_bound
	lower_bound <- matrix(lower_bound, nrow = length(pred), ncol = length(lst), byrow = FALSE)
	upper_bound <- foreach::foreach(i = 1:length(lst), .final = unlist) %do% lst[[i]]$upper_bound
	upper_bound <- matrix(upper_bound, nrow = length(pred), ncol = length(lst), byrow = FALSE)

	if(contiguize){
	lower_bound <- apply(lower_bound, 1, min, na.rm = TRUE)
	lower_bound[which(is.infinite(lower_bound))] <- NA

	upper_bound <- apply(upper_bound, 1, max, na.rm = TRUE)
	upper_bound[which(is.infinite(upper_bound))] <- NA
	return(tibble::tibble(pred = pred, lower_bound = lower_bound, upper_bound = upper_bound))
	}else{
		empirical_lower_bounds <- apply(lower_bound, 2, min, na.rm = TRUE)
		empirical_upper_bounds <- apply(upper_bound, 2, max, na.rm = TRUE)

		contiguous_intervals <- foreach::foreach(i = 1:length(pred)) %do%
			contiguize_intervals(lower_bound[i,], upper_bound[i,], empirical_lower_bounds, empirical_upper_bounds,return_all = T)

		return(tibble::tibble(pred = pred, intervals = contiguous_intervals))


}

	}

#' Contiguize non-contiguous intervals
#'
#' @param pot_lower_bounds Potential non-contiguous lower bounds
#' @param pot_upper_bounds Potential non-contiguous upper bounds
#' @param empirical_lower_bounds Observed lower bounds
#' @param empirical_upper_bounds Observed upper bounds
#' @param return_all Return all intervals or just contiguous intervals
contiguize_intervals <- function(pot_lower_bounds,
																 pot_upper_bounds,
																 empirical_lower_bounds,
																 empirical_upper_bounds,
																 return_all = FALSE){

	if(all(is.na(pot_lower_bounds))){
		return(tibble::tibble(lower_bound = NA, upper_bound = NA))
	}

	intervals <- matrix(c(pot_lower_bounds,pot_upper_bounds, empirical_lower_bounds, empirical_upper_bounds), nrow = length(pot_lower_bounds))
	intervals <- stats::na.omit(intervals)

	i <- 1
	while(i < nrow(intervals) & nrow(intervals) > 1){
		if(intervals[i,2] == intervals[i, 4] & intervals[i+1, 3] == intervals[i+1, 1]){
			intervals[i,2] <- intervals[i+1,2]
			intervals[i,4] <- intervals[i+1,4]
			intervals <- matrix(intervals[-(i+1),], ncol = 4)
		}else{
			i <- i + 1
		}
	}

	widths <- intervals[,2] - intervals[,1]
	if(return_all){
		return(tibble::tibble(lower_bound = as.numeric(intervals[,1]), upper_bound = as.numeric(intervals[,2])))
	}else{
		colnames(intervals) <- c('lower_bound', 'upper_bound', 'empirical_lower_bound', 'empirical_upper_bound')

	return(intervals[which.min(widths)[1],1:2])
	}
}

#' Function to optimize clusters based on Calibrated Clustering
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param method Clustering method to use, either 'ks' for Kolmogorov-Smirnov or 'kmeans' for K-means clustering
#' @param min_m Minimum number of clusters to consider
#' @param max_m Maximum number of clusters to consider. If NULL, defaults to the number of unique classes minus one
#' @param ms Vector of specific numbers of clusters to consider. If NULL, defaults to a sequence from min_m to max_m
#' @param maxit Maximum number of iterations for the clustering algorithm
#' @param q Quantiles to use for K-means clustering, default is a sequence from 0.1 to 0.9 in steps of 0.1
#' @return A vector of cluster assignments, with attributes containing the clusters, coverage gaps, method used, number of clusters, and Calibrated Clustering index
optimize_clusters <- function(ncs, class_vec, method = c('ks', 'kmeans'), min_m = 2, max_m = NULL,
															ms = NULL, maxit = 100, q = seq(0.1, 0.9, by=0.1)){
	method <- match.arg(method, c('ks', 'kmeans'))
m <- NULL
	if(is.null(max_m)){
		max_m <- length(unique(class_vec))-1
	}

	if(is.null(ms)){
		ms <- min_m:max_m
	}

	if(method== 'ks'){
		clusters_pot <- foreach::foreach(m = ms) %do%
			clusterer(ncs, m, class_vec, maxit = maxit, method = 'ks')
	}else if(method == 'kmeans'){
		clusters_pot <- foreach::foreach(m = ms) %do%
			clusterer(ncs, m, class_vec, maxit = maxit, method = 'kmeans', q = q)
	}

	ch_indices <- purrr::map_dbl(clusters_pot, ~ attr(.x,'ch_index'))

clusters <- clusters_pot[[which.max(ch_indices)]]

		return(clusters)

}

#' Function to cluster non-conformity scores using either Kolmogorov-Smirnov or K-means clustering
#' @param ncs Vector of non-conformity scores
#' @param m Number of clusters to form
#' @param class_vec Vector of class labels
#' @param maxit Maximum number of iterations for the clustering algorithm
#' @param method Clustering method to use, either 'ks' for Kolmogorov-Smirnov or 'kmeans' for K-means clustering
#' @param q Quantiles to use for K-means clustering, default is a sequence from 0.1 to 0.9 in steps of 0.1
#' @return A vector of cluster assignments, with attributes containing the clusters, coverage gaps, method used, number of clusters, and Calibrated Clustering index
clusterer <- function(ncs, m, class_vec, maxit = 100, method = c('ks','kmeans'),q = seq(0.1, 0.9, by=0.1)){
	i <- NULL

	method <- match.arg(method, c('ks', 'kmeans'))

	if(method == 'ks'){
	clusters <- ks_cluster(ncs, class_vec, m, maxit)
	}else if(method == 'kmeans'){
		clusters <- kmeans_cluster_qecdf(ncs, class_vec, m = m, q = q)
	}else{
		stop('Unknown clustering method')
	}

	coverage_gaps <- foreach::foreach(i = 1:length(clusters), .final = unlist) %do%
		coverage_gap_finder(ncs, class_vec, clusters[[i]])

	cluster_vec <- rep(NA, length(class_vec))
	for(i in 1:length(clusters)){
		cluster_vec[class_vec %in% clusters[[i]]] <- i
	}

	ch_score <- ch_index(ncs, class_vec, clusters, q = q)

	attributes(cluster_vec) <- list(clusters = clusters,
																	coverage_gaps = coverage_gaps,
																	method = method, m = m,
																	ch_index = ch_score)

	return(cluster_vec)
}


#' Function to convert class vector to cluster vector based on calibrated clusters
#' @param class_vec Vector of class labels
#' @param cluster_vec_calib Vector of calibrated clusters
#' @return A vector of cluster assignments, with attributes containing the clusters, method used, number of clusters, Calibrated Clustering index, and coverage gaps
class_to_clusters <- function(class_vec, cluster_vec_calib){
	i <- NULL

	clusters <- attr(cluster_vec_calib, 'clusters')
	cluster_vec <- rep(NA, length(class_vec))
	for(i in 1:length(clusters)){
		cluster_vec[class_vec %in% clusters[[i]]] <- i
	}

	attributes(cluster_vec) <- list(clusters = clusters,
																	method = attr(cluster_vec_calib, 'method'),
																	m = attr(cluster_vec_calib, 'm'),
																	calib_ch_index = attr(cluster_vec_calib, 'ch_index'),
																	calib_coverage_gaps = attr(cluster_vec_calib, 'coverage_gaps'))
	return(cluster_vec)

}

#' Function to perform Kolmogorov-Smirnov clustering on non-conformity scores
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param m Number of clusters to form
#' @param maxit Maximum number of iterations for the clustering algorithm
#' @param nrep Number of repetitions for the clustering algorithm
#' @return A vector of cluster assignments, with attributes containing the clusters, coverage gaps, method used, number of clusters, and Calibrated Clustering index
ks_cluster <- function(ncs, class_vec, m, maxit = 100, nrep = 10){
	r <- NULL
	class_labels <- unique(class_vec)

	if(m < 2){
		stop('m must be greater than or equal to 2')
	}

	if(length(class_labels) < m){
		stop('Number of classes must be greater than or equal to m')
	}

	ch_score <- 0

	for(r in 1:nrep){

	tmp_clusters <- ks_cluster_init_step(ncs, class_vec, m)

	for(i in 1:maxit){
		clusters <- ks_cluster_assignment_step(ncs, class_vec, class_labels, tmp_clusters)

		if(identical(clusters, tmp_clusters)){
			break
		}else{
			tmp_clusters <- clusters
		}
		if(i == maxit){
			warning('Maximum number of iterations reached without convergence')
		}
	}


	score <- ch_index(ncs, class_vec, clusters)
	if(score > ch_score){
		ch_score <- score
		clusters_final <- clusters
	}
	}


return(clusters_final)
}



#' Function to initialize clusters for Kolmogorov-Smirnov clustering
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param m Number of clusters to form
ks_cluster_init_step <- function(ncs, class_vec, m){
	j <- i <- NULL

	# randomly select 1 start class
	clusters <- list()
	class_labels <- unique(class_vec)
	clusters[[1]] <- sample(class_labels, 1)

	for(j in 2:m){
		dms <- foreach::foreach(i = 1:length(class_labels), .final = unlist) %do%
			Dm_finder(ncs, class_vec, class_labels[i], clusters, return = 'min')

		probs <- foreach::foreach(i = 1:length(class_labels),.final=unlist) %do%
			dm_to_prob(dms[i], dms)

		clusters[[j]] <- sample(class_labels, 1, prob = probs)

	}

	return(clusters)


}

#' Function to assign classes to clusters based on Kolmogorov-Smirnov clustering
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param class_labels Vector of unique class labels
#' @param clusters List of clusters
#' @param m Number of clusters
#' @return A list of clusters, where each element is a vector of class labels assigned to that cluster
ks_cluster_assignment_step <- function(ncs, class_vec, class_labels, clusters, m){
	i <- NULL

	cluster_vec <- foreach::foreach(i = 1:length(class_labels), .final = unlist) %do% {
		Dm_finder(ncs, class_vec, class_labels[i], clusters, return = "which.min")
	}

	return(split(class_labels, cluster_vec))

	}

#' Function to find the minimum distance between a class and a set of clusters
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param class Class label to compare against the clusters
#' @param clusters List of clusters
#' @param return Character string indicating what to return. Options are 'min' for the minimum distance, 'which.min' for the index of the cluster with the minimum distance, or 'vec' for a vector of distances to each cluster.
#' @return A numeric value or vector depending on the value of the `return` parameter. If `return` is 'min', returns the minimum distance. If `return` is 'which.min', returns the index of the cluster with the minimum distance. If `return` is 'vec', returns a vector of distances to each cluster.
Dm_finder <- function(ncs, class_vec, class, clusters, return = c('min','which.min','vec')){
	i <- NULL

	dnorms <- foreach::foreach(c = 1:length(clusters),.final = unlist) %do% {
		d <- suppressWarnings(stats::ks.test(x = ncs[class_vec == class], y = ncs[class_vec %in% clusters[[c]]]))$statistic

d_norm <- sqrt((length(ncs[class_vec == class]) * length(ncs[class_vec %in% clusters[[c]]])) / (length(ncs[class_vec == class]) + length(ncs[class_vec %in% clusters[[c]]]))) * d
		d_norm
	}
	if(return == 'which.min'){
		return(which.min(dnorms))
	}else if(return == 'min'){
		return(min(dnorms))
	}else if(return == 'vec'){
		return(dnorms)
	}
}

#' Function to convert distance measure to probability
#' @param dm Distance measure
#' @param dms Vector of distance measures for all clusters
#' @return A numeric value representing the probability of the distance measure relative to the sum of all distance measures
dm_to_prob <- function(dm, dms){
	dm/sum(dms)
}


#' Function to perform K-means clustering on quantile empirical cumulative distribution functions (qECDFs) of non-conformity scores
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param q Quantiles to use for the qECDFs, default is a sequence from 0.1 to 0.9 in steps of 0.1
#' @param m Number of clusters to form
#' @return A list of clusters, where each element is a vector of class labels assigned to that cluster
kmeans_cluster_qecdf <- function(ncs,class_vec, q = seq(0.1, 0.9, by=0.1), m){
	i <- NULL

	class_labels <- unique(class_vec)
	qecdfs <- foreach::foreach(i = 1:length(class_labels)) %do%
		stats::quantile(ncs[class_vec == class_labels[i]], probs = q, na.rm = TRUE)

	qecdfs <- do.call(rbind, qecdfs)

	clusters_kmeans <- stats::kmeans(qecdfs, centers = m, nstart = 10)$cluster

	return(split(class_labels, clusters_kmeans))
}

#' Function to find the coverage gap for a set of clusters
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param cluster Vector of cluster labels
#' @return A numeric value representing the maximum coverage gap between the clusters
coverage_gap_finder <- function(ncs, class_vec, cluster){
	if(length(cluster) == 1){
		return(0)
	}else{
		i <- j <- NULL
		epsilons <- foreach::foreach(i = 1:(length(cluster)-1),.final=unlist) %:%
			foreach::foreach(j = (i+1):length(cluster),.final = unlist) %do%
				suppressWarnings(
					stats::ks.test(x = ncs[class_vec == cluster[i]], y = ncs[class_vec == cluster[j]])$statistic
				)
		return(max(epsilons, na.rm = TRUE))
	}
}

#' Function to compute the Calinski-Harabasz index for a set of clusters
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param clusters List of clusters, where each element is a vector of class labels assigned to that cluster
#' @param q Quantiles to use for the qECDFs, default is a sequence from 0.1 to 0.9 in steps of 0.1
#' @return A numeric value representing the Calinski-Harabasz index for the clusters
ch_index <- function(ncs,class_vec, clusters, q = seq(0.1, 0.9, by=0.1)){
	if(length(clusters) == 1){
		return(0)
	}

	if(!is.list(clusters)){

	clusters <- attr(clusters,'clusters')
	}
	i <- NULL

	wcss <- foreach::foreach(i = 1:length(clusters), .final = unlist) %do%
		wcss_compute(ncs, class_vec, clusters[[i]], q = q)

	wcss <- sum(wcss, na.rm = TRUE)

	bcss <- bcss_compute(ncs, class_vec, clusters, q = q)

	m <- length(clusters)

	ch <- (bcss/(m-1)) / (wcss/(length(unlist(clusters))-m))
	return(ch)
}


#' Function to compute the within-cluster sum of squares (WCSS) for a set of clusters
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param cluster Vector of cluster labels
#' @param q Quantiles to use for the qECDFs, default is a sequence from 0.1 to 0.9 in steps of 0.1
#' @return A numeric value representing the WCSS for the cluster
wcss_compute <- function(ncs, class_vec, cluster, q = seq(0.1, 0.9, by=0.1)){
	i <- NULL
qs <- foreach::foreach(i = 1:length(cluster),.final = dplyr::bind_rows) %do%
		stats::quantile(ncs[class_vec == cluster[i]], probs = q, na.rm = TRUE)

mean_qs <- stats::quantile(ncs[class_vec %in% cluster], probs = q, na.rm = TRUE)

	return(sum((t(qs) - mean_qs)^2, na.rm = TRUE))

}

#' Function to compute the between-cluster sum of squares (BCSS) for a set of clusters
#' @param ncs Vector of non-conformity scores
#' @param class_vec Vector of class labels
#' @param clusters List of clusters, where each element is a vector of class labels assigned to that cluster
#' @param q Quantiles to use for the qECDFs, default is a sequence from 0.1 to 0.9 in steps of 0.1
#' @return A numeric value representing the BCSS for the clusters
bcss_compute <- function(ncs, class_vec, clusters, q = seq(0.1, 0.9, by=0.1)){
	i <- NULL
	qs <- foreach::foreach(i = 1:length(clusters),.final = dplyr::bind_rows) %do%
		stats::quantile(ncs[class_vec %in% clusters[[i]]], probs = q, na.rm = TRUE)

	mean_qs <- stats::quantile(ncs, probs = q, na.rm = TRUE)

	bcss <- foreach::foreach(i = 1:length(clusters),.final = unlist) %do%{
		length(clusters[[i]])*sum((qs[i,] - mean_qs)^2, na.rm = TRUE)
}
	return(sum(bcss, na.rm = TRUE))
}



