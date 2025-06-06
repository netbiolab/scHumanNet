% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindAllHub.R
\name{FindAllHub}
\alias{FindAllHub}
\title{Find statistically significant hub genes in each scHumanNets}
\usage{
FindAllHub(
  net.list = NULL,
  centrality = "degree",
  q.method = "BH",
  threshold = 0.05
)
}
\arguments{
\item{net.list}{output of SortAddLLS. list of network dataframe}

\item{centrality}{default is "degree" centrality. Available options are: c("degree", "betweenness","closeness", "eigenvector"). Degree is calculated with sum of weights. Weights are uniform in betweenness and closeness centrality.}

\item{q.method}{default "BH", input one of the following c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")}

\item{threshold}{a threshold value to cut the significant hub genes. default 0.05}
}
\value{
a dataframe with Percentile Rank Centrality, gene, pvalue, qvalue and the celltype
}
\description{
Find statistically significant hub genes in each scHumanNets
}
\examples{
sorted.net.list <- SortAddLLS(Celltype.specific.networks, reference.network = graph.hn3)
sig.hub.df <- FindAllHub(net.list = sorted.net.list)
}
