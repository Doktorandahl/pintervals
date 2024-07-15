#' Bin-conditional bootstrap prediction intervals
#'
#' @description
#' This function computes bootstrapped prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using bin-conditional bootstrapped prediction errors. The prediction errors to bootstrap from are computed using either a calibration set with predicted and true values or a set of pre-computed prediction errors from a calibration dataset or other data which the model was not trained on (e.g. OOB errors from a model using bagging). The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals.
#'
#' Currently not working as intended. May be removed in future versions.
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values
#' @param calib_truth A numeric vector of true values in the calibration partition. Only required if calib is a numeric vector
#' @param calib_bins A vector of bin identifiers for the calibration set
#' @param breaks A vector of break points for the bins to manually define the bins. If NULL, lower and upper bounds of the bins are calculated as the minimum and maximum values of each bin in the calibration set. Must be provided if calib_bins or nbins are not provided, either as a vector or as the last column of a calib tibble.
#' @param nbins Automatically chop the calibration set into nbins based on the true values with approximately equal number of observations in each bin. Must be provided if calib_bins or breaks are not provided.
#' @param calib_bin_type A string specicying whether the bins are based on the predicted values ('prediction') or the true values ('truth'). Default is 'prediction'. Ignored if calib_bins is provided.
#' @param error_type The type of error to use for the prediction intervals. Can be 'raw' or 'absolute'. If 'raw', bootstrapping will be done on the raw prediction errors. If 'absolute', bootstrapping will be done on the absolute prediction errors with random signs. Default is 'raw'
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param n_bootstraps The number of bootstraps to perform. Default is 1000
#' @param lower_bound Optional minimum value for the prediction intervals. If not provided, the minimum (true) value of the calibration partition will be used
#' @param upper_bound Optional maximum value for the prediction intervals. If not provided, the maximum (true) value of the calibration partition will be used
#' @param right Parameter passed to cut function to determine which side of the bin interval is closed. Default is TRUE
pinterval_boot_bins <- function(pred,
															 calib,
															 calib_truth = NULL,
															 calib_bins = NULL,
															 breaks = NULL,
															 nbins = NULL,
															 calib_bin_type = c('prediction', 'truth'),
															 error_type = c('raw','absolute'),
															 alpha = 0.1,
															 n_bootstraps=1000,
															 lower_bound = NULL,
															 upper_bound = NULL,
															 right = TRUE){

	i <- NA

	if(!is.numeric(pred)){
		stop('pred must be a single number or a numeric vector')
	}

	if(is.numeric(calib) & is.null(calib_truth)){
		stop('If calib is numeric, calib_truth must be provided')
	}

	if(!is.numeric(calib) && ncol(calib)<2){
		stop('calib must be a numeric vector or a 2 or 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and (optionally) the third column being the bin values if bin structure is not provided in argument bins')
	}

	if((is.null(breaks)) && is.null(nbins) && (is.null(calib_bins) | ncol(calib_bins)!=3)){
		stop('If breaks for bins or nbins are not provided, bins for the calibration set must be provided as a vector or a as the last column of the calib if calib is a tibble or matrix')
	}

	calib_bin_type <- match.arg(calib_bin_type,c('prediction', 'truth'))


	if(!is.numeric(calib)){
		calib_org <- calib
		if(is.matrix(calib)){
			calib <- as.numeric(calib_org[,1])
			calib_truth <- as.numeric(calib_org[,2])
			if(is.null(calib_bins)){
				calib_bins <- as.numeric(calib_org[,3])
			}
		}else{
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
			if(is.null(calib_bins)){
				calib_bins <- as.numeric(calib_org[[3]])
			}
		}
	}

	if(is.null(calib_bins)){
		if(!is.null(breaks)){
			if(calib_bin_type == 'prediction'){
				calib_bins <- cut(calib,breaks = breaks,labels = FALSE,right = right)

				if(breaks[1] != -Inf){
					breaks <- c(-Inf,breaks)
				}
				if(breaks[length(breaks)] != Inf){
					breaks <- c(breaks,Inf)
				}
			}else if(calib_bin_type == 'truth'){
			calib_bins <- cut(calib_truth,breaks = breaks,labels = FALSE,right = right)
			breaks <- foreach::foreach(i = 1:length(unique(calib_bins)),.final = unlist) %do%
				max(calib[calib_bins==i])

			breaks <- c(-Inf,breaks[1:(length(breaks)-1)],Inf)
			}
			if(!is.null(nbins)){
				warning('Both breaks and nbins provided, using breaks')
			}
		}else{
			if(!(nbins == round(nbins)) || length(nbins)!=1){
				stop('nbins must be a single integer value')
			}
			if(calib_bin_type == 'prediction'){
				calib_bins <- bin_chopper(calib,nbins = nbins)
				breaks <- foreach::foreach(i = 1:length(unique(calib_bins)),.final = unlist) %do%
					max(calib[calib_bins==i])
				breaks <- c(-Inf,breaks[1:(length(breaks)-1)],Inf)
			}else if(calib_bin_type == 'truth'){
			calib_bins <- bin_chopper(calib_truth,nbins = nbins)
			breaks <- foreach::foreach(i = 1:length(unique(calib_bins)),.final = unlist) %do%
				max(calib[calib_bins==i])
			breaks <- c(-Inf,breaks[1:(length(breaks)-1)],Inf)
			}
		}
	}



	nobs_bins <- as.numeric(table(calib_bins))

	pred_bins <- cut(pred,breaks = breaks,labels = FALSE,right = TRUE)

	if(any(nobs_bins*alpha/2<1)){
		warning('Some bins have few observations use when bootstrapping the prediction intervals at the specified alpha level. Consider using a larger calibration set, a smaller number of bins, or a smaller alpha level or interpret with caution')
	}

	bin_labels <- sort(unique(calib_bins))

	if(length(bin_labels)<2){
		stop('Calibration set must have at least two bins. For continuous bootstrap prediction intervals without bins, use pinterval_bootstrap() instead of pinterval_boot_bin()')
	}


	boot_intervals <- foreach::foreach(i = 1:length(bin_labels)) %do%
		suppressWarnings(pinterval_bootstrap(pred = pred,
																	 lower_bound = lower_bound,
																	 upper_bound = upper_bound,
																	 calib = calib[calib_bins==bin_labels[i]],
																	 calib_truth = calib_truth[calib_bins==bin_labels[i]],
																	 alpha = alpha))

	boot_intervals <- foreach::foreach(i = 1:length(pred)) %do%
		boot_intervals[[pred_bins[i]]][i,]


	if(is.null(lower_bound)){
		lb <- -Inf
	}else{
		lb <- lower_bound
	}
	if(is.null(upper_bound)){
		ub <- Inf
	}else{
		ub <- upper_bound
	}

	boot_intervals <- dplyr::bind_rows(boot_intervals)
	boot_intervals <- boot_intervals %>%
		dplyr::mutate(lower_bound = dplyr::case_when(.data$lower_bound<lb ~ lb,
																					TRUE ~ .data$lower_bound),
								 upper_bound = dplyr::case_when(.data$upper_bound>ub ~ ub,
																					TRUE ~ .data$upper_bound))

	return(boot_intervals)





}
