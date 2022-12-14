% Generated by roxygen2: do not edit by hand
% Please edit documentation in create-independencepvalue.Rmd.
\name{sample_psi}
\alias{sample_psi}
\title{Sample from the distribution of eigenvalues of W * inv(T)}
\usage{
sample_psi(p, rp, n)
}
\arguments{
\item{p}{\eqn{p_1+p_2}}

\item{rp}{\eqn{min(p_1, p_2)}}

\item{n}{sample size}
}
\value{
A vector of length \code{rp} sampled from the joint distribution described
above.
}
\description{
Samples from the \code{rp}-dimensional joint distribution of the eigenvalues of
\eqn{WT^{-1}}, where \eqn{W} and \eqn{T} are independent Wisharts with dimensions specified
in Prop 1(ii).  These are the \eqn{\Psi_i}, and taking \eqn{\sqrt{(\Psi_i/(1+\Psi_i))}}
gives a sample from the joint distribution of the canonical correlations
between two groups of variables of size \eqn{p_1} and \eqn{p_2} under the null.
}
\keyword{internal}