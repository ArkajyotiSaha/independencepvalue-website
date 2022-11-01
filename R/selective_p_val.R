# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Test independence of a data-dependent group of variables with the rest
#'
#' Given a covariance matrix `S` of `p` Gaussian variables and a grouping obtained
#' via thresholding absolute correlations at `1-c` using `block_diag()`, this 
#' function tests the null hypothesis of independence between two groups of 
#' Gaussian variables.
#' 
#' @param S a \eqn{p \times p} covariance matrix
#' @param CP a vector of length \eqn{p} with \eqn{i^{th}} element denoting the 
#' group \eqn{i^{th}} variable belongs to
#' @param k the group to be tested for independence with the remaining variables, i.e. \eqn{P = [i : CP[i]==k]}
#' @param n sample size
#' @param c a threshold
#' @param d0 a natural number; if the number of canonical correlations is greater than \code{d0}, Monte Carlo simulation will be used to approximate the p-value for computational convenience; default value is 5
#' @param maxeval the maximum number of function evaluations used to approximate the p-value using `selective_p_val_integrate()`; we recommend using a high value of this to obtain an approximation with high accuracy; default value is 10,000
#' @param tol the relative tolerance used to approximate the p-value using `selective_p_val_integrate()`; default value is 1e-05
#' @param mc_iter the number of Monte Carlo iterations used to approximate the p-value; we recommend using a high value of this to obtain an approximation with high accuracy; default value is 1,000
#' @return The selective p-value for the test of independence.
#' @examples
#' # Simulates a 10 x 5 X from N(0, I)
#' set.seed(1)
#' X <- matrix(rnorm(50), 10, 5)
#'
#' # Compute the correlation matrix of X.
#' corX <- cor(X)
#' # Use 'block_diag' to obtain any block diagonal structure
#' block_diag_structure <- block_diag(corX, c= 0.5)
#' # test for independence of the variables in group 1 with the remaining variables
#' selective_p_val(S=cov(X), n=10, CP=block_diag_structure, c=0.5, 
#' k=1, d0=5, tol = 1e-05, maxeval = 10000, mc_iter=100)
#' @export 
selective_p_val <- function(S, CP, k, n, c, d0 = 5, tol = 1e-05, maxeval = 1e5, mc_iter = 1000){
  test_hyp <- test_stat_CCA(S, CP, k)
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  if (p2 == 1 & p2 <= d0) {
    du <- selective_p_val_beta(S, CP, k, n, c, test_hyp)
  }
  else {
    du <- 0
    L <- form_L(test_hyp)
    g <- c(rep(c, 2 * p1 * p2), rep(0, p2), rep(1, p2))
    if (p2 <= d0) {
      du <- selective_p_val_integrate(n, L, g, test_hyp, tol, maxeval)
    }
    if (du <= 0 || du >= 1 || p2 > d0) {
      # use Monte Carlo approach
      du <- selective_p_val_MC(n, L, g, test_hyp, mc_iter)
    }
  }
  return(du)
}
