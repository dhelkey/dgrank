% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperfunctions.R
\name{linCholSolver}
\alias{linCholSolver}
\title{Solve System of Equations w/ Cholesky}
\usage{
linCholSolver(R, y)
}
\arguments{
\item{R}{- Upper triangular matrix of dimension <n x p>; Cholesky decomposition of A.}

\item{y}{- Vector of length n.}
}
\description{
\code{linCholSolver} solves Ax = y for x.
It is first required to take the Cholesky decomposition,
obtaining A = t(R) %*% R. This must be performed outside of this function for speed.
Results should be checked against a known solver, as this function optimizes for speed.
}
