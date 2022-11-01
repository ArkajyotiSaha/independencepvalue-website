# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Compute selective p-value using Monte Carlo
#' 
#' @param n sample size
#' @param L matrix used to define conditioning set
#' @param g vector used to define conditioning set
#' @param test_hyp output of `test_stat_CCA()`
#' @param mc_iter number of Monte Carlo iterations
#' 
#' @keywords internal
selective_p_val_MC <- function(n, L, g, test_hyp, mc_iter) {
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  p <- p1 + p2
  sip <-
    future.apply::future_sapply(1:mc_iter, function(i)
      MC_function_selective(p, p2, n, L, g), future.seed = TRUE)
  zp <- sum(as.numeric(sip[2, ]))
  if (zp < 100) {
    sip <-
      future.apply::future_sapply(1:(min((mc_iter * 100 / zp), 100000)), function(i)
        MC_function_selective(p, p2, n, L, g), future.seed = TRUE)
  }
  mean(test_hyp$statistic >= sip[1, ][sip[2, ] == TRUE])
}
