#' Gaussian RBF kernel
#' 
#' Computes the Gaussian RBF kernel k(x,x') = exp(-|x-x'|^2/(2*sigma^2)).
#' 
#' @param x Input matrix of covariates with samples in rows.
#' @param sigma Bandwidth of the Gaussian RBF kernel. Default is 1.
#' @param x2 (Optional) a second matrix of covariates with samples in rows.
#' 
#' @return If \code{x2} is provided, returns the kernel matrix crossing samples
#'   in \code{x} and samples in \code{x2}. If \code{x2} is not provided, returns the kernel
#'   Gram matrix of \code{x} versus itself.
#'   
#' @export
#' 

gausskernel <- function(x, 
                        sigma = 1, 
                        x2)
{
  if (missing(x2)) {
    x2 <- x
  }
  
  n <- nrow(x)
  m <- nrow(x2)
  i <- rep(seq_len(n), times=m)
  j <- rep(seq_len(m), each=n)
  u <- x[i, , drop = FALSE] - x2[j, , drop = FALSE]
  kernel <- matrix(exp(- rowSums(u ^ 2) / (2 * sigma^2)), n, m)
  dimnames(kernel) <- list(rownames(x),rownames(x2))
  
  return(kernel)
}
