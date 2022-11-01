# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to evaluate the joint density of the canonical correlations
#'
#' This function evaluates the joint density of the canonical correlations at a specific value \code{lambda}.
#' 
#' @param n sample size
#' @param p total number of variables
#' @param rp smaller of the size of the two groups of variables
#' @param a a scaling constant
#' @param lambda a \code{rp} length vector where the joint density is to be evaluated \eqn{1 \ge lambda[1]\ge lambda[2] \ge ... \ge lambda[rp] \ge 0} 
#' @return Scaled joined density of canonical correlations, evaluated at \eqn{lambda}.
#' @keywords internal
dCCA <- function(n, p, rp, a, lambda){
  c1 <- (p - 2*rp) * sum(log(lambda))
  c2 <- ((n - p - 2)/2) * sum(log(1-lambda^2))
  dip <- matrix(1, length(lambda), length(lambda))
  dip[lower.tri(dip)] <- stats::dist(lambda^2)
  c3 <- sum(log(dip))
  return(exp(a+c1+c2+c3))
}
