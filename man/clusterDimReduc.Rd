% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterDimReduc.R
\name{clusterDimReduc}
\alias{clusterDimReduc}
\title{Performs LSI, clustering and UMAP}
\usage{
clusterDimReduc(
  archr.proj,
  varFeatures = 15000,
  resolution = c(0.2, 1, 2),
  ...
)
}
\arguments{
\item{archr.proj}{ArchRProject object}

\item{varFeatures}{Number of variable features to use for cluster}

\item{resolution}{Vector of clustering resolutions for multi-round clusters}

\item{...}{Additional arguments for addIterativeLSI}
}
\value{
None
}
\description{
Sets global variables about specific genome build and number of threads to utilze on the machine.
}
\keyword{internal}
