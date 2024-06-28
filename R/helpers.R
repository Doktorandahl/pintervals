





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
grid_finder <- function(y_min,y_max,ncs,ncs_function,y_hat, alpha, min_step = NULL, grid_size = NULL){
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
								 	grid_inner(ncs_function(y_hat[i],pos_vals),y_hat[i],ncs,pos_vals,alpha)

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
grid_inner <- function(hyp_ncs,y_hat,ncs,pos_vals,alpha){
	if(sum(hyp_ncs<stats::quantile(ncs,1-alpha))==0){
		return(c(pred = as.numeric(y_hat), lower_bound = NA_real_, upper_bound = NA_real_))
	}else{
		lb <- min(pos_vals[hyp_ncs<stats::quantile(ncs,1-alpha)])
		ub <- max(pos_vals[hyp_ncs<stats::quantile(ncs,1-alpha)])

		return(c(pred = as.numeric(y_hat), lower_bound = lb, upper_bound = ub))
	}
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
	nobs_per_value <- table(x)
	binsizes <- rep(target_num, nbins)
	binsizes2 <- rep(0, nbins)
	cutpoints <- rep(0, nbins-1)
	k <- 0
	while(!identical(binsizes,binsizes2) & k<10){
	for(i in 1:nbins){
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
	binsizes2 <- binsizes
	k <- k + 1
	}

if(!return_breaks){
	return(cut(x, breaks = c(-Inf,cutpoints, Inf), labels = FALSE))
}else{
	return(c(-Inf,cutpoints, Inf))
}
	}







