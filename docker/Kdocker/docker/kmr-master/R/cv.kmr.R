#' Cross-validation for KMR
#' 
#' Does a k-fold cross-validation for \code{kmr}, and returns a fitted KMR model, CV performance scores and optimal values for the regularization parameter \code{lambda}.
#' 
#' @param x \code{x} input matrix as in \code{kmr}.
#' @param y \code{y} output matrix as in \code{kmr}.
#' @param kx_type Kernel type for observations as in \code{kmr}.
#' @param kx_option Optional list of parameters for the observation kernel as in \code{kmr}.
#' @param kt_type Kernel type for tasks as in \code{kmr}.
#' @param kt_option Optional list of parameters for the task kernel as in \code{kmr}.
#' @param lambda Sequence of (more than one) values for lambda that must be tested. Default is 10^(-5:5).
#' @param type.measure Measure type for evaluating performance. Possible options are 
#'   "\code{ci}" (concordance index), "\code{mse}" (mean squared error), 
#'   "\code{cor}" (pearson correlation). Default is "\code{ci}".
#' @param nfolds Number of folds for cross-validation. Default is 5.
#' @param nrepeats Number of times the k-fold cross-validation is performed. Default is 1.
#' @param seed A seed number for the random number generator (useful to have the same CV splits).
#' @param mc.cores Number of parallelable CPU cores to use.
#' 
#' @return An object (list) of class \code{"cv.kmr"}, which can then be used to make predictions for the different tasks on new observations, as a list containing the following slots:
#' \item{...}{Outputs of a CV-fitted KMR model as in \code{"kmr"}.}
#' \item{meanCV}{A matrix of CV performance scores of dim ntask x nlambda.}
#' \item{bestlambda}{A vector of lambdas of length ntask, each corresp to the underlying min CV score.}
#' \item{bestCV}{A vector of min CV performance scores of length ntask.}
#' \item{lambda}{Lambda sequence against which a model is tested.}
#' \item{type.measure}{Measure type for evaluating performance.}
#' 
#' @importFrom parallel mclapply
#' 
#' @export
#' 
#' @references 
#' Bernard, E., Jiao, Y., Scornet, E., Stoven, V., Walter, T., and Vert, J.-P. (2017). Kernel multitask regression for toxicogenetics. \href{https://doi.org/10.1101/171298}{bioRxiv-171298}.
#' 
#' @seealso \code{\link{kmr}}
#' 
#' @examples 
#' # Data
#' ntr <- 80
#' ntst <- 20
#' nt <- 50
#' p <- 20
#' xtrain <- matrix(rnorm(ntr*p),ntr,p)
#' xtest <- matrix(rnorm(ntst*p),ntst,p)
#' ytrain <- matrix(rnorm(ntr*nt),ntr,nt)
#' ytest <- matrix(rnorm(ntst*nt),ntst,nt)
#' 
#' # Train with Gaussian RBF kernel for x and multitask kernel for t
#' cvmo <- cv.kmr(x=xtrain, y=ytrain, kx_type="gaussian", kx_option=list(sigma=100), 
#'                kt_type="multitask", kt_option=list(alpha=0.5), type.measure="mse")
#' # Plot
#' plot(cvmo)
#' # Predict
#' cvpred <- predict(cvmo, xtest)
#' # Evaluate
#' evalpred(cvpred, ytest, "mse")
#' 

cv.kmr <- function(x, 
                   y, 
                   kx_type = c("linear", "gaussian", "precomputed"), 
                   kx_option = list(sigma = 1), 
                   kt_type = c("multitask", "empirical", "precomputed"), 
                   kt_option = list(alpha = 1, kt = NULL), 
                   lambda = 10^(-5:5), 
                   type.measure = c("ci", "mse", "cor"), 
                   nfolds = 5, 
                   nrepeats = 1, 
                   seed = 9182456, 
                   mc.cores = 1)
{
  kx_type <- match.arg(kx_type)
  kt_type <- match.arg(kt_type)
  type.measure <- match.arg(type.measure)
  N <- nrow(x)
  Nt <- ncol(y)
  Nl <- length(lambda)
  stopifnot(Nl > 1)
  
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  # Set random number generator seed
  set.seed(seed)
  
  ###  Make folds
  folds <- list()
  for (i in seq(nrepeats)) {
    folds <- c(folds, split(sample(seq(N)), rep(1:nfolds, length = N)))
  }
  nexp <- length(folds)
  
  ### Iterate over the folds
  resCV <- mclapply(seq(nexp) , function(iexp) {
    message('.', appendLF = F)
    
    # training samples for the fold
    mytrain <- seq(N)[-folds[[iexp]]]
    
    # test samples for the fold
    mytest <- folds[[iexp]]
    
    # Train the model and make the prediction
    if (kx_type == "precomputed") {
      xtrain <- x[mytrain,mytrain]
      xtest <- x[mytest,mytrain]
    } else {
      xtrain <- x[mytrain,]
      xtest <- x[mytest,]
    }
    
    # Train on the training set
    m <- kmr(xtrain, y[mytrain,,drop=F], kx_type, kx_option, kt_type, kt_option)
    # Predict on the test set
    ypred <- predict.kmr(m, xtest, lambda=lambda)
    ytest <- y[mytest,,drop=F]
    return(evalpred(ypred, ytest, type.measure))
  }, mc.cores = mc.cores)
  
  meanCV <- Reduce("+", resCV) / length(resCV)
  which.lambda <- switch(type.measure,
                         "ci" = which.max,
                         "mse"= which.min,
                         "cor" = which.max)
  ilambda <- apply(meanCV,1,which.lambda)
  bestlambda <- lambda[ilambda]
  bestCV <- meanCV[cbind(seq_along(ilambda),ilambda)]
  
  ### Train model on full data
  res <- kmr(x, y, kx_type, kx_option, kt_type, kt_option)
  
  res[['meanCV']] <- meanCV
  res[['bestlambda']] <- bestlambda
  res[['bestCV']] <- bestCV
  res[['lambda']] <- lambda
  res[['type.measure']] <- type.measure
  
  class(res) <- "cv.kmr"
  return(res)
}
