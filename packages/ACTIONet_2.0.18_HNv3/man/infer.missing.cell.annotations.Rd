% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_tools.R
\name{infer.missing.cell.annotations}
\alias{infer.missing.cell.annotations}
\title{Uses a variant of the label propagation algorithm to infer missing labels}
\usage{
infer.missing.cell.annotations(
  ace,
  initial_labels,
  iters = 3,
  lambda = 0,
  sig_threshold = 3
)
}
\arguments{
\item{ace}{Input results to be clustered
(alternatively it can be the ACTIONet igraph object)}

\item{initial_labels}{Annotations to correct with missing values (NA) in it.
It can be either a named annotation (inside ace$annotations) or a label vector.}

\item{double.stochastic}{Whether to densify adjacency matrix before running label propagation (default=FALSE).}

\item{max_iter}{How many iterative rounds of correction/inference should be performed (default=3)}

\item{adjust.levels}{Whether or not re-adjust labels at the end}
}
\value{
ace with updated annotations added to ace$annotations
}
\description{
Uses a variant of the label propagation algorithm to infer missing labels
}
\examples{
ace = infer.missing.cell.annotations(ace, sce$assigned_archetypes, 'updated_archetype_annotations')
}
