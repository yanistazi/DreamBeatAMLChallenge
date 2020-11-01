#' Make prediction with a "\code{kmr}" object
#' 
#' Similar to other \code{predict} methods, this function predicts fitted values
#' from a fitted "\code{kmr}" object.
#' 
#' @param object Fitted "\code{kmr}" model object.
#' @param newx Matrix of new values for \code{x} or kernel matrix for new values
#'  crossing old values for \code{x}, at which predictions are to be made.
#' @param lambda Value(s) of the regularization parameter, a single scalar or a
#'   sequence of values. Default is 1.
#' @param ... Not used. Other arguments to predict.
#'   
#' @return A matrix of predicted values for the new samples (in rows) and all
#'   tasks (in columns), corresponding to the regularization parameter
#'   \code{lambda}. If \code{lambda} is a list of values, the function returns a
#'   list of matrices, corresponding to the predictions for the different values
#'   of \code{lambda} in the list.
#' 
#' @export
#' 
#' @seealso \code{\link{kmr}}
#' 

predict.kmr <- function(object, 
                        newx, 
                        lambda = 1, 
                        ...)
{
  # Kernel matrix between test and train sets
  K = switch(object$kx_type,
             "precomputed" = newx,
             "linear" = newx %*% t(object$x),
             "gaussian" = gausskernel(newx, object$kx_option$sigma, object$x))
  
  predict.kmr.singlelambda <- function(singlelambda) {
    S <- (object$Dx %*% t(object$Dt) + singlelambda)^(-1)
  
    # offset b
    b <- object$Ut %*% ( ( (t(S) %*% (object$gamma*object$gamma))^{-1} )*( ( t(S) * (t(object$Ut)%*%t(object$y)%*%object$Ux) ) %*% object$gamma ) )
  
    # weights alpha
    alpha <- object$Ux %*% ( S * ( t(object$Ux) %*% ( object$y - matrix(1,nrow=nrow(object$y),ncol=1) %*% t(b) ) %*% object$Ut ) ) %*% t(object$Ut)
  
    # KRR Prediction
    mypred <- K %*% alpha %*% object$Kt + matrix(1,nrow=nrow(newx),ncol=1)%*%t(b)
    return(mypred)
  }
  
  # Make predictions for all values of lambda
  pred <- lapply(lambda, predict.kmr.singlelambda)
  
  # Return the list of prediction matrices if there are several lambda's, or
  # just the prediction matrix if there is a single lambda
  if (length(pred) == 1) {
    return(pred[[1]])
  } else {
    return(pred)
  }
}