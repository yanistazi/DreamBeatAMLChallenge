#' Evaluate prediction performance
#' 
#' Computes performance scores of predicted against true responses.
#' 
#' @param ypred Matrix or list of matrices (e.g. each corresponding to a tested lambda)
#'   of predicted reponses, of dimension \code{nobs x ntask}.
#' @param ytest Matrix of true reponses, of dimension \code{nobs x ntask}.
#' @param type.measure Measure type for evaluating performance. Possible options are 
#'   "\code{ci}" (concordance index), "\code{mse}" (mean squared error), 
#'   "\code{cor}" (pearson correlation). Default is "\code{ci}".
#' 
#' @return A matrix of mean performance scores (averaged over observations),
#'   of dimension \code{ntask x nlambda} where \code{nlambda} denotes the length of the list \code{ypred}.
#' 
#' @useDynLib kmr cidx
#' 
#' @note \code{ypred} and \code{ytest} are not interchangably equivalent for \code{type.measure="ci"},
#'   as \code{ytest} is treated as true responses and used to determine denominator in computing C-index.
#' 
#' @export
#' 
#' @importFrom stats cor
#' 

evalpred <- function(ypred, 
                     ytest, 
                     type.measure = c("ci", "mse", "cor"))
{
  type.measure <- match.arg(type.measure)
  
  errfun <- switch(type.measure,
              "ci" = 
                function(u) {
                  apply(rbind(u,ytest), 2, function(v) {
                    if (isTRUE(requireNamespace("Hmisc", quietly = TRUE))) {
                      Hmisc::rcorr.cens(v[1:nrow(u)], v[-(1:nrow(u))], outx = FALSE)[[1]]
                    } else {
                      .Call("cidx", v[1:nrow(u)], v[-(1:nrow(u))])
                    }
                  })
                },
              "mse" = 
                function(u) {
                  apply( ((u-ytest)^2), 2, mean )
                },
              "cor" = function(u) {
                diag(cor(u, ytest))
              }
  )
  
  if (!is.list(ypred))
    ypred <- list(ypred)
  
  return(sapply(ypred, errfun))
}
