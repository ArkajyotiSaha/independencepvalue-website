# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function for Monte Carlo simulation for classical inference
#'
#' Samples from the joint distribution of the canonical correlations between two groups of independent Gaussian variables of size \eqn{p_1} and \eqn{p_2} using \code{sample_psi()}, and then computes the corresponding test statistic, which follows a Wilks' lambda distribution. 
#' 
#' @param p \eqn{p_1+p_2}
#' @param rp \eqn{min(p_1, p_2)}
#' @param n sample size
#' @return A sample from Wilks' lambda distribution.
#' @keywords internal
MC_function_classical <- function(p, rp, n){
  F_X_eigenvalues <- sample_psi(p, rp, n) # eigenvalues of (inv(W)T)
  while(any(F_X_eigenvalues < 0)){
    F_X_eigenvalues <- sample_psi(p, rp, n) # eigenvalues of (inv(W)T)
  }
  statistic <- 1/prod(1+F_X_eigenvalues)
  # prod(1 - lambda_i^2) = prod(1 - Psi_i/(1 + Psi_i)) = prod(1/(1 + Psi_i))
  return(statistic)
}
