#' Plot the cross-validation curve produced by "\code{cv.kmr}"
#' 
#' Plots the cross-validation curve as a function of \code{lambda} values tested, 
#' as well as the optimal \code{lambda} selected by cross-validation for each task.
#' 
#' @param x Fitted "\code{cv.kmr}" model object.
#' @param ... Other graphical parameters to \code{matplot}.
#' 
#' @return A plot is produced, and nothing returned.
#' 
#' @export
#' 
#' @import graphics
#' 
#' @seealso \code{\link{cv.kmr}}
#' 

plot.cv.kmr <- function(x, 
                        ...)
{
  cvobj <- x
  matplot(log(cvobj$lambda), t(cvobj$meanCV), type = "l", lwd = 2, xlab = "log(lambda)", ylab = cvobj$type.measure, ...)
  points(log(cvobj$bestlambda), cvobj$bestCV, lwd = 2)
  invisible()
}
