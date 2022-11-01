# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Sample from the distribution of eigenvalues of W * inv(T)
#' 
#' Samples from the \code{rp}-dimensional joint distribution of the eigenvalues of 
#' \eqn{WT^{-1}}, where \eqn{W} and \eqn{T} are independent Wisharts with dimensions specified 
#' in Prop 1(ii).  These are the \eqn{\Psi_i}, and taking \eqn{\sqrt{(\Psi_i/(1+\Psi_i))}}
#' gives a sample from the joint distribution of the canonical correlations 
#' between two groups of variables of size \eqn{p_1} and \eqn{p_2} under the null.
#' 
#' @param p \eqn{p_1+p_2}
#' @param rp \eqn{min(p_1, p_2)}
#' @param n sample size
#' @return A vector of length \code{rp} sampled from the joint distribution described
#' above.
#' @keywords internal
sample_psi <- function(p, rp, n) {
  tilde_W_X <- stats::rWishart(1, p - rp, diag(rp))#simulate W
  tilde_T_X <- stats::rWishart(1, n - (p - rp) - 1, diag(rp))#simulate T
  tilde_F_X <- tilde_W_X[,,1] %*% solve(tilde_T_X[,,1])# inv(W)T
  return(eigen(tilde_F_X)$values) # eigenvalues of (inv(W)T)
}
