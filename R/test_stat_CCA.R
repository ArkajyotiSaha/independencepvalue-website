# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Compute the test statistic as a function of canonical correlations. 
#' 
#' Given a sample covariance matrix \eqn{S} and a group of variables \eqn{P}, first
#' computes the cross-covariance matrix between the whitened variables: 
#' \eqn{S_{P, P^c}^W = S_{P, P}^{-0.5} S_{P, P^c} S_{P^c, P^c}^{-0.5}}. Next, 
#' computes the SVD of \eqn{S_{P, P^c}^W} and returns the test statistic, 
#' \eqn{S_{P, P}}, \eqn{S_{P^c, P^c}} and the canonical vectors 
#' of \eqn{S_{P, P^c}}.
#' 
#' @param S a \eqn{p \times p} sample covariance matrix
#' @param CP a vector of length \eqn{p} with \eqn{i^{th}} element denoting the 
#' group \eqn{i^{th}} variable belongs to
#' @param k the group to be tested for independence with the remaining variables, i.e. \eqn{P = i : CP[i]==k}
#' @return A list containing the following items:
#' \item{statistic}{Test statistic corresponding to \eqn{S} and group of variables \eqn{P}.}
#'
#' \item{S11}{\eqn{S_{P, P}} if \eqn{2|P| \ge p}, else \eqn{S_{P^c, P^c}}.}
#'
#' \item{S22}{\eqn{S_{P^c, P^c}} if \eqn{2|P| \ge p}, else \eqn{S_{P, P}}.}
#'
#' \item{left_SV}{Left canonical vectors of  \eqn{S_{P, P^c}}.}
#'
#' \item{right_SV}{Right canonical vectors of \eqn{S_{P, P^c}}.}
#'
#' @keywords internal
test_stat_CCA <- function(S, CP, k) {
  p <- nrow(S)
  ptemp <- sum(CP==k)
  if(2*ptemp >= p){
    S11_x <- as.matrix(S[which(CP==k), which(CP==k)])#S11
    S22_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S22
    S12_x <- as.matrix(S[which(CP==k), which(CP!=k)])#S_12 
  }
  if(2*ptemp < p){
    S11_x <- as.matrix(S[which(CP!=k), which(CP!=k)])#S11
    S22_x <- as.matrix(S[which(CP==k), which(CP==k)])#S22
    S12_x <- as.matrix(S[which(CP!=k), which(CP==k)])#S_12 
  }
  #ensures that p2 <= p1.
  S11_x_half <- amen::mhalf(S11_x)#S_11^{1/2}
  inv_S11_x_half <- solve(S11_x_half)#S_11^{-1/2}
  S22_x_half <- amen::mhalf(S22_x)#S_11^{1/2}
  inv_S22_x_half <- solve(S22_x_half)#S_11^{-1/2}
  tilde_S12_x <- inv_S11_x_half %*% S12_x %*% inv_S22_x_half#S_12^W = covariance matrix of whitened X_1 and X_2
  svdecom <- svd(tilde_S12_x)#compact SVD 
  singular_values <- svdecom$d#lambda
  test_stat <- prod(1-singular_values^2)#test statistic 
  L_x <- S11_x_half %*% svdecom$u
  R_x <- t(svdecom$v) %*% S22_x_half
  return(list(statistic=test_stat, S11 = S11_x, S22_x = S22_x, left_SV=L_x, right_SV=t(R_x)))
}
