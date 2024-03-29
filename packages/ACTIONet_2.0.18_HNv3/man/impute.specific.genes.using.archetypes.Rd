% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputation.R
\name{impute.specific.genes.using.archetypes}
\alias{impute.specific.genes.using.archetypes}
\title{Imputing expression specificity of genes by interpolating over archetype profile}
\usage{
impute.specific.genes.using.archetypes(ace, genes)
}
\arguments{
\item{ace}{ACTIONet output}

\item{genes}{List of genes to impute}
}
\value{
A matrix of imputed expression values
}
\description{
Imputing expression specificity of genes by interpolating over archetype profile
}
\examples{
expression_imputed = impute.genes.using.archetype(ace, genes)
}
