# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function for Monte Carlo simulation for selective inference
#'
#' Samples from the joint distribution of the canonical correlations between two groups of independent Gaussian variables of size \eqn{p_1} and \eqn{p_2} using \code{sample_psi()}, and then computes the corresponding test statistic, which follows a Wilks' lambda distribution. Given a matrix \code{L} with \code{rp} columns, and a vector of length equal to the number of rows in \code{L}, checks if the sampled vector \eqn{\lambda} satisfies \eqn{\lambda : L\lambda \le g}
#' 
#' @param p \eqn{p_1 + p_2}
#' @param rp \eqn{min(p_1, p_2)}
#' @param n sample size
#' @param L a matrix with \code{rp} columns
#' @param g a vector of length equal to the number of rows in \code{L}
#' @return
#' \item{statistic}{Test statistic corresponding to simulated \eqn{\lambda}.}
#' \item{status}{A logical vector indicating if the simulated \eqn{\lambda} satisfies \eqn{L\lambda \leq g}.}
#' @keywords internal
MC_function_selective <- function(p, rp, n, L, g){
  if(nrow(L)!=length(g)){stop("error: number of rows of matrix L must be equal to the number of rows of vector g")}
  F_X_eigenvalues <- sample_psi(p, rp, n) # eigenvalues of (inv(W)T)
  while(any(F_X_eigenvalues<0)){
    F_X_eigenvalues <- sample_psi(p, rp, n) # eigenvalues of (inv(W)T)
  }
  statistic <- 1/prod(1+F_X_eigenvalues)
  # prod(1 - lambda_i^2) = prod(1 - Psi_i/(1 + Psi_i)) = prod(1/(1 + Psi_i))
  status <- all(L %*% sqrt(F_X_eigenvalues/(1+F_X_eigenvalues)) <= g)
  return(list(statistic=statistic, status=status))
}
