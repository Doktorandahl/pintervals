#' Calibrate predictions using a calibration dataset
#'
#' @description
#' This function calibrates a vector of predicted values using a held-out calibration set of predictions and true outcomes. It supports both parametric calibration via generalized linear models (`glm`) and non-parametric calibration via isotonic regression. The function returns a vector of calibrated predictions.
#'
#'
#' @param pred Vector of predicted values
#' @param calib A numeric vector of predicted values in the calibration partition or a 2 or 3 column tibble or matrix with the first column being the predicted values and the second column being the truth values and the third, optional, column being the bin labels. If calib is a numeric vector, calib_truth and either calib_bins or breaks must be provided.
#' @param calib_truth A numeric vector of true values in the calibration partition
#' @param method Method for calibration, either 'glm' for generalized linear model or 'isotonic' for isotonic regression
#' @param family Family for the generalized linear model passed to `glm`. Default is 'gaussian' for a linear model. Other options include 'binomial', 'poisson', etc.
#' @param transform Optional transformation function, or string matching a function, to apply to the predicted values before calibration. If NULL, no transformation is applied. See details below.
#'
#' @details
#' This function calibrates model predictions to better match the distribution of observed outcomes, improving reliability especially when prediction error varies systematically. When `method = "glm"`, a generalized linear model is fit using the calibration set. When `method = "isotonic"`, a monotonic fit is estimated to map predicted values to expected outcomes. Platt scaling is a special case of calibration using `glm` with a logistic link function, which is useful for binary classification tasks.
#'
#' If a transformation function is provided via `transform`, it is applied to both the prediction vector (`pred`) and the calibration predictions (`calib`) before calibration is performed. This is useful in certain cases, for exampel if the predictions are on a probability scale, it can be beneficial to apply a logit transformation before calibration.
#'
#' @return A numeric vector of calibrated predictions, aligned with the input vector `pred`. If method = 'glm' is used, the regression coefficients, family, and transformation from the fitted model can be accessed as attributes of the returned vector.
#'
#' @export
#'
#' @examples
#' # Generate synthetic example
#' set.seed(123)
#' x <- runif(1000)
#' y <- rlnorm(1000, meanlog = 1 + x, sdlog = 0.5)
#' df <- tibble::tibble(x, y)
#'
#' # Split into training, calibration, and test sets
#' df_train <- df[1:500, ]
#' df_cal <- df[501:750, ]
#' df_test <- df[751:1000, ]
#'
#' # Fit a model and predict
#' model <- lm(log(y) ~ x, data = df_train)
#' pred_cal <- exp(predict(model, newdata = df_cal))
#' pred_test <- exp(predict(model, newdata = df_test))
#'
#' # Calibrate using isotonic regression
#' calibrated_iso <- calibrate_predictions(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   method = "isotonic"
#' )
#'
#' #' # Calibrate using GLM
#' calibrated_glm <- calibrate_predictions(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   method = "glm",
#'   family = "gaussian")
#'
#' # Calibrate using GLM with a log transform
#' calibrated_glm <- calibrate_predictions(
#'   pred = pred_test,
#'   calib = pred_cal,
#'   calib_truth = df_cal$y,
#'   method = "glm",
#'   family = "gaussian",
#'   transform = log
#' )
calibrate_predictions <- function(pred, calib, calib_truth = NULL,
																	method = c('glm','isotonic'),
																	family = 'gaussian',
																	transform = NULL){

	if(is.numeric(calib) && is.null(calib_truth)){
		stop('If calib is numeric, calib_truth must be provided')
	}

	if(!is.numeric(pred)){
		stop('pred must be a numeric vector')
	}

	if(!is.numeric(calib) && ncol(calib) != 2){
		stop('calib must be a numeric vector or a 2 column tibble or matrix with the first column being the predicted values and the second column being the truth values')
	}

	method <- match.arg(method, c('glm', 'isotonic'))

	if(is.null(calib_truth)){
		calib_org <- calib
		if(is.matrix(calib)){
			calib <- as.numeric(calib_org[,1])
			calib_truth <- as.numeric(calib_org[,2])
		}else{
			calib_truth <- as.numeric(calib_org[[2]])
			calib <- as.numeric(calib_org[[1]])
		}
	}

	if(method == 'glm'){
		if(is.null(transform)){
			transform <- function(x) x
		}else{
			if(is.character(transform)){
				transform <- match.fun(transform)
			} else if(!is.function(transform)){
				stop('transform must be a function or a character string matching a function')
			}

		}
		calib <- tibble::tibble(predicted = transform(calib), truth = calib_truth)
		model <- stats::glm(truth ~ predicted, data = calib, family = family)

		if(stats::coefficients(model)[2] < 0){
			warning("The calibration model has a negative slope, which indicates a negative relationship between the predicted and true values. This may indicate that the model predicts poorly or that the calibration set is not representative of the true distribution of outcomes.")
		}

		calibrated_predictions <- stats::predict(model, newdata = tibble::tibble(predicted = transform(pred)), type = 'response')
	} else if(method == 'isotonic'){
		isotonic_model <- stats::as.stepfun(stats::isoreg(calib, calib_truth))
		calibrated_predictions <- isotonic_model(pred)
	}

	attr(calibrated_predictions, "calibration_method") <- method

	if(method == 'glm'){
		attr(calibrated_predictions, "coefficients") <- model$coefficients
		attr(calibrated_predictions, "family") <- family
		if(!is.null(transform)){
			attr(calibrated_predictions, "transform") <- transform
		}
	}

	return(calibrated_predictions)



}



