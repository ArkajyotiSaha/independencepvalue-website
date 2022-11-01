# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to test the independence of two pre-specified groups of variables
#'
#' Given a covariance matrix \eqn{S} of \eqn{p} Gaussian variables, and a pre-specified group
#' of variables \eqn{P}, this function tests the null hypothesis of independence between the groups of 
#' variables in \eqn{P} and \eqn{P^c}. Makes use of \code{test_stat_CCA()} and \code{sample_psi()}.
#' 
#' @param S a \eqn{p \times p} covariance matrix
#' @param CP a vector of length \eqn{p} with \eqn{i^{th}} element denoting the 
#' group \eqn{i^{th}} variable belongs to
#' @param k the group to be tested for independence with the remaining variables, i.e. \eqn{P = [i : CP[i]==k]}
#' @param n sample size
#' @param mc_iter the number of Monte Carlo iterations used to approximate the p-value; we recommend using a high value of this to obtain an approximation with high accuracy; default value is 1,000
#' @return The p-value for the test of independence. 
#' @examples
#' # Simulates a 10 x 3 X_1 from N(0, I)
#' set.seed(1)
#' X_1 <- matrix(rnorm(30), 10, 3)
#'
#' # Simulates a 10 x 2 X_2 from N(0, I) independently of X_1
#' set.seed(2)
#' X_2 <- matrix(rnorm(20), 10, 2)
#'
#' # Compute the covariance matrix of X = (X_1 X_2).
#' covX <- cov(cbind(X_1, X_2))
#' # tests for a difference in means between X_1 and X_2
#' classical_p_val(S=covX, CP=rep(1:2, times=c(3, 2)), k=1, n=10, mc_iter=100)
#' @export
classical_p_val <- function(S, CP, k, n, mc_iter= 1000){
  test_hyp <- test_stat_CCA(S, CP, k)
  p <- nrow(S)
  rp <- nrow(test_hyp$S22)
  if(rp == 1){
    classic_p_val <- 1 - stats::pbeta(1-test_hyp$statistic, (p - 1)/2, (n - p)/2) 
    # for r(p) = 1, we have more efficient way to compute p_value based on beta distribution.
  }
  if(rp > 1){
    sip <- future.apply::future_sapply(1:mc_iter, function(i) MC_function_classical(p, rp, n), future.seed = TRUE)
    #MC simulation to approximate the p-value. 
    classic_p_val <- mean(test_hyp$statistic >= sip)#computes p-value. The test statistic is smaller if it is away from null. 
  }
  return(classic_p_val)
}
