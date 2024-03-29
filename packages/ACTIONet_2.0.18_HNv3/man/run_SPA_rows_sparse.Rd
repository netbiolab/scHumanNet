% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{run_SPA_rows_sparse}
\alias{run_SPA_rows_sparse}
\title{Runs Successive Projection Algorithm (SPA) to solve separable NMF}
\usage{
run_SPA_rows_sparse(A, k)
}
\arguments{
\item{A}{Input matrix}

\item{k}{Number of columns to select}
}
\value{
A named list with entries 'selected_columns' and 'norms'
}
\description{
Runs Successive Projection Algorithm (SPA) to solve separable NMF
}
\examples{
H = run_SPA(S_r, 10)
}
