#' Bin-conditional conformal prediction intervals for continuous predictions
#'
#'@description
#'This function calculates bin-conditional conformal prediction intervals with a confidence level of 1-alpha for a vector of (continuous) predicted values using inductive conformal prediction on a bin-by-bin basis. The intervals are computed using either a calibration set with predicted and true values or a set of pre-computed non-conformity scores from the calibration set. In addition the function requires either a set of breaks or a vector of bin identifiers for the calibrations set, either as a standalone vector or as the third column of the calibration dataset if the calibration data is provided as a tibble. The function returns a tibble containing the predicted values along with the lower and upper bounds of the prediction intervals. Bin-conditional conformal prediction intervals are useful when the prediction error is not constant across the range of predicted values and ensures that the coverage is (approximately) correct for each bin under the assumption that the non-conformity scores are exchangeable within each bin.
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 or 3 column tibble or matrix with the first column being the predicted values and the second column being the truth values and (optionally) the third column being the bin values if bins are not provided as a standalone vector or if breaks are not provided
#' @param calib_truth A numeric vector of true values in the calibration partition
#' @param calib_bins A vector of bin identifiers for the calibration set
#' @param breaks A vector of break points for the bins to manually define the bins. If NULL, lower and upper bounds of the bins are calculated as the minimum and maximum values of each bin in the calibration set. Must be provided if calib_bins or nbins are not provided, either as a vector or as the last column of a calib tibble.
#' @param nbins Automatically chop the calibration set into nbins based on the true values with approximately equal number of observations in each bin. Must be provided if calib_bins or breaks are not provided.
#' @param alpha The confidence level for the prediction intervals. Must be a single numeric value between 0 and 1
#' @param ncs_function A function or a character string matching a function that takes two arguments, a vector of predicted values and a vector of true values, in that order. The function should return a numeric vector of nonconformity scores. Default is 'absolute_error' which returns the absolute difference between the predicted and true values.
#' @param ncs An optional numeric vector of pre-computed nonconformity scores from a calibration partition. If provided, calib will be ignored. If provided, bins must be provided in calib_bins as well.
#' @param min_step The minimum step size for the grid search. Default is 0.01. Useful to change if predictions are made on a discrete grid or if the resolution of the interval is too coarse or too fine.
#' @param grid_size Alternative to min_step, the number of points to use in the grid search between the lower and upper bound. If provided, min_step will be ignored.
#' #' @param right Logical, if TRUE the bins are right-closed (a,b] and if FALSE the bins are left-closed `[ a,b)`. Only used if breaks or nbins are provided.
#'
#' @return A tibble with the predicted values, the lower and upper bounds of the prediction intervals
#' @export
#'
#' @examples
pinterval_cp_bins = function(pred,
									 calib = NULL,
									 calib_truth = NULL,
									 calib_bins = NULL,
									 breaks = NULL,
									 nbins = NULL,
									 alpha = 0.1,
									 ncs_function = c('absolute_error','squared_error'),
									 ncs = NULL,
									 min_step = 0.01,
									 grid_size = NULL,
									 right = TRUE){

	i <- NA

	if(!is.numeric(pred)){
		stop('pred must be a numeric vector')
	}

	if(is.null(calib) & is.null(ncs)){
		stop('Either calib or ncs must be provided')
	}

	if(is.null(ncs)){
		if(is.numeric(calib) & is.null(calib_truth)){
			stop('If calib is numeric, calib_truth must be provided')
		}
		if(!is.numeric(calib) && ncol(calib)<2){
			stop('calib must be a numeric vector or a 2 or 3 column tibble or matrix with the first column being the predicted values, the second column being the truth values, and (optionally) the third column being the bin values if bin structure is not provided in argument bins')
		}

		if((is.null(breaks)) && is.null(nbins) && (is.null(calib_bins) || (!is.numeric(calib) && ncol(calib)!=3))){
			stop('If breaks for bins or nbins are not provided, bins for the calibration set must be provided as a vector or a as the last column of the calib if calib is a tibble or matrix')
		}
	}else{
		if(!is.numeric(ncs)){
			stop('ncs must be a numeric vector')
		}
		if(!is.null(calib)){
			warning('ncs provided, calib will be ignored')
		}

	if(is.null(calib_bins)){
		stop('calib_bins must be provided when ncs is provided')
	}
	}


	if(!is.numeric(alpha) || alpha<=0 || alpha>=1 || length(alpha)!=1){
		stop('alpha must be a single numeric value between 0 and 1')
	}

	if(is.character(ncs_function)){
		ncs_function <- match.arg(ncs_function, c('absolute_error','squared_error'))
	}

	if(ncs_function == 'absolute_error'){
		ncs_function <- abs_error
	}else if(ncs_function == 'squared_error'){
		ncs_function <- squared_error
	}else if(is.character(ncs_function)){
		ncs_function <- match.fun(ncs_function)
	}else if(!is.function(ncs_function) & is.null(ncs)){
		stop('ncs_function must be a function or a character string matching a function if ncs is not provided. Note, the ncs_function must take two arguments, a vector of predicted values and a vector of true values, in that order')
	}

	if((is.null(breaks))){
		warning('No explicit bin structure provided, breaks are calculated based on the calibration set')
	}

	if(!is.null(breaks) & !is.null(calib_bins)){
		warning("If breaks are provided, calib_bins will be ignored")
	}

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
			calib_bins <- cut(calib_truth,breaks = breaks,labels = FALSE,right = right)
			if(!is.null(nbins)){
				warning('Both breaks and nbins provided, using breaks')
			}
		}else{
			if(!(nbins == round(nbins)) || length(nbins)!=1){
				stop('nbins must be a single integer value')
			}
			calib_bins <- bin_chopper(calib_truth,nbins = nbins)
		}
	}

	nobs_bins <- as.numeric(table(calib_bins))

	if(any(nobs_bins*alpha/2<1)){
		warning('Some bins have too few observations to calculate prediction intervals at the specified alpha level. Consider using a larger calibration set or a smaller alpha level')
	}

		bin_labels <- unique(calib_bins)
		if(length(bin_labels)<2){
			stop('Calibration set must have at least two bins. For continuous prediction intervals without bins, use pinterval_cp_cont() instead of pinterval_cp_bins()')
		}

		lower_bounds <- foreach::foreach(i = bin_labels,.final = unlist) %do% min(calib_truth[calib_bins==i])
		upper_bounds <- foreach::foreach(i = bin_labels,.final = unlist) %do% max(calib_truth[calib_bins==i])


	if(is.null(ncs)){
		ncs <- ncs_function(calib,calib_truth)
	}


	cp_intervals <- foreach::foreach(i = 1:length(bin_labels)) %do%
		suppressWarnings(pinterval_cp_cont(pred = pred,
																	 lower_bound = lower_bounds[i],
																	 upper_bound = upper_bounds[i],
																	 ncs = ncs[calib_bins==bin_labels[i]],
																	 alpha = alpha, min_step = min_step))

	cp_intervals2 <- dplyr::bind_cols(cp_intervals,.name_repair = 'unique_quiet')

	cp_intervals2 <- cp_intervals2 %>% dplyr::rowwise() %>%
		dplyr::mutate(pred = unique(dplyr::c_across(dplyr::starts_with('pred'))),
					 lower_bound = min(dplyr::c_across(dplyr::starts_with('lower_bound')),na.rm = TRUE),
					 upper_bound = max(dplyr::c_across(dplyr::starts_with('upper_bound')),na.rm=TRUE)) %>%
					 	dplyr::ungroup() %>%
					 	dplyr::select(pred,lower_bound,upper_bound) %>%
		dplyr::mutate(lower_bound = dplyr::case_when(is.infinite(.data$lower_bound) ~ -Inf,
													TRUE ~ .data$lower_bound))


	return(cp_intervals2)
}
