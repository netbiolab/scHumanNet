% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reduction.R
\name{reduce.ace}
\alias{reduce.ace}
\title{Takes a `ACTIONetExperiment` object and adds the reduced kernel matrix}
\usage{
reduce.ace(
  ace,
  reduced_dim = 50,
  max_iter = 10,
  assay_name = "logcounts",
  reduction_slot = "ACTION",
  SVD_algorithm = 0,
  seed = 0
)
}
\arguments{
\item{ace}{ACTIONetExperiment object.}

\item{reduced_dim}{Dimension of SVD used for reducing kernel matrix}

\item{max_iter}{Number of SVD iterations}

\item{assay_name}{Name of assay to reduce.}

\item{reduction_slot}{Name of slot to store reduction.}
}
\value{
ACTIONetExperiment object with reduction in colMaps(ace).
}
\description{
Takes a `ACTIONetExperiment` object and adds the reduced kernel matrix
}
\examples{
ace = import.ace.from.10X(input_path)
ace = reduce.ace(ace)
}
