% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{import.ace.from.10X.generic}
\alias{import.ace.from.10X.generic}
\title{Imports data from a 10X experiment folder and constructs an `SingleCellExeriment` object}
\usage{
import.ace.from.10X.generic(
  input_path,
  mtx_file = "matrix.mtx.gz",
  feature_metadata = "features.tsv.gz",
  sample_metadata = "barcodes.tsv.gz",
  sep = "\\t",
  use.names = TRUE,
  prefilter = FALSE,
  ...
)
}
\arguments{
\item{input_path}{Folder containing input files.}

\item{mtx_file}{Count file in Matrix Market format (default='matrix.mtx.gz').}

\item{feature_metadata}{Table of the same size as number of rows in the count matrix (default='features.tsv.gz').}

\item{sample_metadata}{Table of the same size as number of columns in the count matrix (default='barcodes.tsv.gz').}

\item{sep}{Column-separator used in the row/column annotations files (default='\\t').}

\item{prefilter}{Whether to prefilter rows/columns of input counts matrix. Must specify filtering parameters to pass to filter.ace().}
}
\value{
`SingleCellExeriment` object
}
\description{
Imports data from a 10X experiment folder and constructs an `SingleCellExeriment` object
}
\examples{
ace = import.ace.from.10X.generic(input_path, prefilter=TRUE, min_feats_per_cell = 500)
}
