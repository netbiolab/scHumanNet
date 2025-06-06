% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{run_simplex_regression}
\alias{run_simplex_regression}
\title{Solves min_{X} (|| AX - B ||) s.t. simplex constraint}
\usage{
run_simplex_regression(A, B, computeXtX = FALSE)
}
\arguments{
\item{A}{Input matrix}

\item{B}{Input matrix}
}
\value{
X Solution
}
\description{
Solves min_{X} (|| AX - B ||) s.t. simplex constraint
}
\examples{
C = ACTION.out$C[[10]]
A = S_r \%*\% C
B = S_r
H = run_simplex_regression(A, B)
}
