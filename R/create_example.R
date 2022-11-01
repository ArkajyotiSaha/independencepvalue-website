# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Function to create a population correlation matrix for the vignette
#' 
#' This function creates a `p` x `p` matrix where the diagonal blocks are two p/2 equicorrelation matrix with parameter `a`. The all the elements in the offdiagonal blocks of the matrix are equal to `b`. `p` has to be an even number. 
#' 
#' @param p total number of variables, must be even
#' @param a parameter for the equicorrelation matrix on the diagonal blocks
#' @param b parameter for off-diagonal blocks
#' @return a `p` x `p` matrix.
#' @export
create_example <- function(p, a, b) {
  amat <- matrix(a,p/2,p/2)
  bmat <- matrix(b,p/2,p/2)
  Sigma <- rbind(cbind(amat,bmat),cbind(bmat,amat))
  diag(Sigma) <- 1
  Sigma
}
