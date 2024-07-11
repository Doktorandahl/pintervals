





#' Absolute Error
#'
#' @param pred a numeric vector of predicted values
#' @param truth a numeric vector of true values
#'
#' @return a numeric vector of absolute errors
#'
abs_error <- function(pred, truth){
	return(abs(pred-truth))
}

squared_error <- function(pred, truth){
	return((pred-truth)^2)
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
#'
#' @return a tibble with the predicted values and the lower and upper bounds of the prediction intervals
#'
grid_finder <- function(y_min,y_max,ncs,ncs_function,y_hat, alpha, min_step = NULL, grid_size = NULL,return_min_q=FALSE,weighted_cp = FALSE, calib = NULL){
	i <- NA
	if(is.null(grid_size)){
		pos_vals <- seq(from=y_min,to=y_max,by=min_step)
		if(length(pos_vals)>10000){
			warning('Grid size with set step size is large, consider adjusting min_step or using grid_size instead of min_step if the search is too slow')
		}
	}else{
		pos_vals <- seq(from=y_min,to=y_max,length.out=grid_size)
	}

	if(weighted_cp){
	out <- foreach::foreach(i = 1:length(y_hat)) %do%
								 	grid_inner(ncs_function(y_hat[i],pos_vals),y_hat[i],ncs,pos_vals,alpha,return_min_q, weights = weights_calculator(y_hat[i], calib))
	}else{
		out <- foreach::foreach(i = 1:length(y_hat)) %do%
			grid_inner(ncs_function(y_hat[i],pos_vals),y_hat[i],ncs,pos_vals,alpha,return_min_q)
	}
	return(dplyr::bind_rows(out))
}



#' Inner function for grid search
#'
#' @param hyp_ncs vector of hypothetical non-conformity scores
#' @param y_hat predicted value
#' @param ncs vector of non-conformity scores
#' @param pos_vals vector of possible values for the lower and upper bounds of the prediction interval
#' @param alpha confidence level
#'
#' @return a numeric vector with the predicted value and the lower and upper bounds of the prediction interval
grid_inner <- function(hyp_ncs,y_hat,ncs,pos_vals,alpha,return_min_q=FALSE, weights = NULL){
	if(!is.null(weights)){
	ncs <- sample(ncs, size = length(ncs), replace = TRUE, prob = weights)
	}
	if(sum(hyp_ncs<stats::quantile(ncs,1-alpha))==0){
		if(return_min_q){
			min_q <- stats::ecdf(ncs)(min(hyp_ncs))
			return(c(pred = as.numeric(y_hat), lower_bound = NA_real_, upper_bound = NA_real_, min_q = min_q))
		}else{
		return(c(pred = as.numeric(y_hat), lower_bound = NA_real_, upper_bound = NA_real_))
		}
	}else{
		lb <- min(pos_vals[hyp_ncs<stats::quantile(ncs,1-alpha)])
		ub <- max(pos_vals[hyp_ncs<stats::quantile(ncs,1-alpha)])

		if(return_min_q){

			min_q <- stats::ecdf(ncs)(min(hyp_ncs))

			return(c(pred = as.numeric(y_hat), lower_bound = lb, upper_bound = ub, min_q = min_q))
		}else{
			return(c(pred = as.numeric(y_hat), lower_bound = lb, upper_bound = ub))
		}
	}
}

weights_calculator <- function(y_hat, calib){
	similarity <- abs(y_hat-calib)
	weights <- 1-((similarity-min(similarity))/(max(similarity)-min(similarity)))
	weights <- length(calib)*weights/sum(weights)
	return(weights)

}

bootstrap_inner <- function(pred, error, nboot, alpha, lower_bound, upper_bound){
	i <- NA
	boot_error <- sample(error, size = nboot, replace = TRUE)
	boot_pred <- pred + boot_error
	lb <- as.numeric(stats::quantile(boot_pred, alpha/2))
	ub <- as.numeric(stats::quantile(boot_pred, 1-alpha/2))
	if(lb < lower_bound){
		lb <- lower_bound
	}
	if(ub > upper_bound){
		ub <- upper_bound
	}

	return(c(pred = as.numeric(pred), lower_bound = lb, upper_bound = ub))
}


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

minq_to_alpha <- function(minq, alpha){
	minq[which(minq > alpha)] <- alpha
	return(minq)
}

flatten_cp_bin_intervals <- function(lst,
																		 treat_noncontiguous = c('narrowest', 'most_conformal','full')){


	pred <- lst[[1]]$pred
	lower_bound <- foreach::foreach(i = 1:length(lst), .final = unlist) %do% lst[[i]]$lower_bound
	lower_bound <- matrix(lower_bound, nrow = length(pred), ncol = length(lst), byrow = FALSE)
	upper_bound <- foreach::foreach(i = 1:length(lst), .final = unlist) %do% lst[[i]]$upper_bound
	upper_bound <- matrix(upper_bound, nrow = length(pred), ncol = length(lst), byrow = FALSE)

	if(treat_noncontiguous == 'full'){
	lower_bound <- apply(lower_bound, 1, min, na.rm = TRUE)
	lower_bound[which(is.infinite(lower_bound))] <- NA

	upper_bound <- apply(upper_bound, 1, max, na.rm = TRUE)
	upper_bound[which(is.infinite(upper_bound))] <- NA
	return(tibble::tibble(pred = pred, lower_bound = lower_bound, upper_bound = upper_bound))
	}else if(treat_noncontiguous == 'narrowest'){
		empirical_lower_bounds <- apply(lower_bound, 2, min, na.rm = TRUE)
		empirical_upper_bounds <- apply(upper_bound, 2, max, na.rm = TRUE)

		contiguous_intervals <- foreach::foreach(i = 1:length(pred),.final = bind_rows) %do%{
			contiguize_intervals(lower_bound[i,], upper_bound[i,], empirical_lower_bounds, empirical_upper_bounds)
		}
		return(dplyr::bind_cols(tibble::tibble(pred = pred), contiguous_intervals))
	}else if(treat_noncontiguous == 'most_conformal'){
		min_qs <- foreach::foreach(i = 1:length(lst), .final = unlist) %do% lst[[i]]$min_q
		min_qs <- matrix(min_qs, nrow = length(pred), ncol = length(lst), byrow = FALSE)
		bin_minq <- apply(min_qs, 1, which.min, na.rm = TRUE)
		lower_bound <- apply(lower_bound, 1, function(x) x[bin_minq])
		upper_bound <- apply(upper_bound, 1, function(x) x[bin_minq])

		return(tibble::tibble(pred = pred, lower_bound = lower_bound, upper_bound = upper_bound))



	}




	}

contiguize_intervals <- function(pot_lower_bounds,
																 pot_upper_bounds,
																 empirical_lower_bounds,
																 empirical_upper_bounds){

	intervals <- matrix(c(pot_lower_bounds,pot_upper_bounds, empirical_lower_bounds, empirical_upper_bounds), nrow = length(pot_lower_bounds))

	i <- 1
	while(i < nrow(intervals)){
		if(intervals[i,2] == intervals[i, 4] & intervals[i+1, 3] == intervals[i+1, 1]){
			intervals[i,2] <- intervals[i+1,2]
			intervals[i,4] <- intervals[i+1,4]
			intervals <- intervals[-(i+1),]
		}else{
			i <- i + 1
		}
	}

	intervals[,5] <- intervals[,2] - intervals[,1]
	colnames(intervals) <- c('lower_bound', 'upper_bound', 'empirical_lower_bound', 'empirical_upper_bound', 'width')
	return(intervals[which.min(intervals[,5])[1],1:2])
}





