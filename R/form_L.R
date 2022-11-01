# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Form L matrix described in Theorem 2
#' 
#' @param test_hyp output of `test_stat_CCA()`
#' 
#' @keywords internal
form_L <- function(test_hyp) {
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  
  tildeU <- diag(1/sqrt(diag(test_hyp$S11))) %*% test_hyp$left_SV
  if(p2 > 1){
    tildeV <- diag(1/sqrt(diag(test_hyp$S22))) %*% test_hyp$right_SV
  }
  if(p2 == 1){
    tildeV <- 1/sqrt(diag(test_hyp$S22)) * test_hyp$right_SV
  }
  mult_list <- future.apply::future_lapply(1:p2, function(j) {res<- data.table::copy(tildeU);  as.data.frame(collapse::setop(res, "*", tildeV[j,], rowwise = T))}, future.seed = TRUE)
  L1 <- data.table::rbindlist(mult_list)
  L <- matrix(0, (2 * p1 * p2 + 2 * p2), p2)
  L[1:(2 * p1 * p2),] <- as.matrix(rbind(L1, -L1))
  if(p2 > 1){
    index <- cbind(2 * p1 * p2 + c(1:(p2 - 1)), c(1:(p2 - 1)))
    L[index] <- -1
    index[,2] <- index[,2] + 1
    L[index] <- 1
  }
  L[2 * p1 * p2 + p2, p2] <- -1
  index <- cbind(2 * p1 * p2 + p2 + c(1:(p2)), c(1:(p2)))
  L[index] <- 1
  return(L)
}
