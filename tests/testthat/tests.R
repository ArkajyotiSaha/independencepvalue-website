# Generated from create-independencepvalue.Rmd: do not edit by hand  
testthat::test_that("test_stat_CCA() works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(
    round(test_stat_CCA(S=cov(X), CP=rep(1:2, times=c(3, 2)), k=1)$statistic, 2),
    0.62)
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  # testing group 2 should give identical results:
  testthat::expect_equal(
    round(test_stat_CCA(S=cov(X), CP=rep(1:2, times=c(3, 2)), k=2)$statistic, 2),
    0.62)
})

testthat::test_that("MC_function_classical() works", {
  set.seed(123)
  p <- 5
  n <- 20
  nsim <- 1e4
  from_wisharts <- sapply(1:nsim,
                          function(i) MC_function_classical(p = p, rp = 1, n = n))
  from_beta <- 1 - rbeta(nsim, (p - 1) / 2, (n - p) / 2)
  probs <- seq(0.05, 0.95, length = 10)
  qw <- quantile(from_wisharts,probs = probs)
  qb <- quantile(from_beta,probs = probs)
  testthat::expect_true( all(abs(qw - qb) < 0.01) )
})

testthat::test_that("classical_p_val() works", {
  testthat::skip_on_cran()
  nsim <- 1e3
  simulation_classical <- function(i){
    set.seed(i)
    X <- matrix(rnorm(50), 10, 5)
    return(classical_p_val(S=cov(X), CP=rep(1:2, times=c(3, 2)), k=1, n=10, mc_iter=100))
  }
  classical_pval <- future.apply::future_sapply(1:nsim, simulation_classical, future.seed = TRUE)
  probs <- seq(0.05, 0.95, length = 10)
  qw <- quantile(classical_pval, probs = probs)
  testthat::expect_true(all(abs(qw - probs) < 0.025) )
})

testthat::test_that("classical_p_val() works for r=1", {
  testthat::skip_on_cran()
  nsim <- 1e4
  simulation_classical <- function(i){
    set.seed(i)
    X <- matrix(rnorm(50), 10, 5)
    return(classical_p_val(S=cov(X), CP=rep(1:2, times=c(1, 4)), k=1, n=10, mc_iter=100))
  }
  classical_pval <- future.apply::future_sapply(1:nsim, simulation_classical, future.seed = TRUE)
  probs <- seq(0.05, 0.95, length = 10)
  qw <- quantile(classical_pval, probs = probs)
  testthat::expect_true(all(abs(qw - probs) < 0.005) )
})

testthat::test_that("block_diag() works", {
  set.seed(1)
  X <- matrix(rnorm(50), 10, 5)
  testthat::expect_equal(length(unique(block_diag(cor(X), c=0.5))), 3)
})

testthat::test_that("selective_p_val() works", {
  testthat::skip_on_cran()
  nsim <- 1e3
  set.seed(1)
  simulation_selective <- function(i){
    X <- matrix(rnorm(50), 10, 5)
    corX <- cor(X)
    block_diag_structure <- block_diag(corX, c=0.5)
    if(length(unique(block_diag_structure)) > 1){
      k0 <- sample(unique(block_diag_structure), 1)
      return(selective_p_val(S=cov(X), CP=block_diag_structure, k=k0, n=10, c=0.5, d0=2, maxeval = 1000, mc_iter=100))
    }
  }
  selective_pval <- sapply(1:nsim, simulation_selective)
  probs <- seq(0.05, 0.95, length = 10)
  qw <- quantile(unlist(selective_pval), probs = probs)
  testthat::expect_true(all(abs(qw - probs) < 0.025))
})

testthat::test_that("selective_p_val_beta() and selective_p_val_MC() produce similar results", {
  testthat::skip_on_cran()
  set.seed(4)
  X <- matrix(rnorm(600), 30, 20)
  corX <- cor(cbind(X))
  block_diag_structure <- block_diag(corX, c=0.3)
  table(block_diag_structure)
  # We test group 6 as r = 1 here. Setting d0 >= 1 will make selective_p_val() use selective_p_val_beta() for approximation
  p_beta <- selective_p_val(S=cov(X), CP=block_diag_structure, k=6, n=30, c=0.3, d0=5)
  set.seed(1)
  # Setting d0 = 0 forces selective_p_val() to use selective_p_val_MC() for approximation, even if r = 1
  p_MC <- selective_p_val(S=cov(X), CP=block_diag_structure, k=6, n=30, c=0.3, d0=0, mc_iter=5e4)
  testthat::expect_true(abs(p_beta - p_MC) < 0.002)
})

testthat::test_that("selective_p_val_integrate() and selective_p_val_MC() produce similar results", {
  testthat::skip_on_cran()
  set.seed(4)
  X <- matrix(rnorm(600), 30, 20)
  corX <- cor(cbind(X))
  block_diag_structure <- block_diag(corX, c=0.3)
  table(block_diag_structure)
  # We test group 2 as r = 2 here. Setting d0 >= 2 will make selective_p_val() use selective_p_val_integrate() for approximation. 
  p_integrate <- selective_p_val(S=cov(X), CP=block_diag_structure, k=2, n=30, c=0.3, d0=5, maxeval=2e5, mc_iter=0)
  set.seed(1)
  # Setting d0 = 0 forces selective_p_val() to use selective_p_val_MC() for approximation, even if r = 1
  p_MC <- selective_p_val(S=cov(X), CP=block_diag_structure, k=2, n=30, c=0.3, d0=0, mc_iter=1e5)
  testthat::expect_true(abs(p_integrate - p_MC) < 0.0015)
})

