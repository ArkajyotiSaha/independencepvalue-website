# Generated from create-independencepvalue.Rmd: do not edit by hand

#' Compute selective p-value using numerical integration
#' 
#' @param n sample size
#' @param L matrix used to define conditioning set
#' @param g vector used to define conditioning set
#' @param test_hyp output of `test_stat_CCA()`
#' @param maxeval the maximum number of function evaluations used to approximate the p-value using `selective_p_val_integrate()`; we recommend using a high value of this to obtain an approximation with high accuracy; default value is 10,000
#' @param tol the relative tolerance used to approximate the p-value using `selective_p_val_integrate()`; default value is 1e-05
#' 
#' @keywords internal
selective_p_val_integrate <- function(n, L, g, test_hyp, tol = 1e-05, maxeval = 1e4) {
  p1 <- nrow(test_hyp$S11)
  p2 <- nrow(test_hyp$S22)
  p <- p1 + p2

  du <- 0 # initialize du
  P <- rcdd::makeH(L, g, x = NULL) # create the Half space representation of the polytope
  PV_d <- rcdd::scdd(P) # compute the convex hull representation of the polytope.
  V_d <- as.matrix(PV_d$output[,-c(1, 2)])
  Pi <- volesti::Vpolytope(V = V_d) # create a convenient Vertex representation from the convex hull. 
  triang = try(geometry::delaunayn(Pi@V), silent = TRUE) # try the Delaunay triangulation
  if (!inherits(triang, 'try-error')) { # if the triangulation is successful
    prod_res <- 1
    for (i in 1:p2) {
      prod_res <-
        prod_res * gamma((n - i) / 2) / (gamma((n - p2 -  i) / 2) * gamma((p1 - i + 1) /
                                                                            2) * gamma((p2 - i + 1) / 2))
    }
    alpha <- log(pi ^ (p2 / 2) * 2 ^ p *  prod_res) # compute the constant
    f_tot <- function(x) {
      dt <- dCCA(n, p, p2, alpha, x)
      return(c(dt, dt * (prod(1 - x ^ 2) <= test_hyp$statistic))) # compute the function in the numerator
    }
    part_int <- function(i, Pi, triang, f_tot) { #code to perform integration
      if (stats::var(round(Pi@V[triang[i, ],][, 1], 5)) == 0) {
        Pi@V[triang[i, ],][, 1][length(Pi@V[triang[i, ],][, 1])] <-
          Pi@V[triang[i, ], ][, 1][length(Pi@V[triang[i, ],][, 1])] * (1 - 10 ^ (-5)) # Account for numerical instabilities on the triangulation by adding negligble deviation
      }
      return(SimplicialCubature::adaptIntegrateSimplex(f_tot, t(Pi@V[triang[i, ],]), fDim = 2, tol = tol, maxEvals = maxeval)$integral) # perform numerical integration on the simplices
    }
    par_I_tot_list <-
      future.apply::future_sapply(1:nrow(triang), part_int, Pi, triang, f_tot, future.seed =
                                    TRUE)
    du <- sum(par_I_tot_list[2, ]) / sum(par_I_tot_list[1, ])
    if (sum(par_I_tot_list[1, ]) == 0) { # if the integral in the denominator is effectively zero, set du to 0.
      du <- 0
    }
  }
  return(du)
}
