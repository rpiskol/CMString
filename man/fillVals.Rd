% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMString_helperFunctions.R
\name{fillVals}
\alias{fillVals}
\title{add missing genes to expression matrix. Gene expression for these genes in a particular sample will be filled with the average of gene expression across all genes in that sample}
\usage{
fillVals(sdat, missingGenes, addJitter = FALSE, verbose = FALSE)
}
\arguments{
\item{sdat}{gene expression matrix to which missing gene expression will be added (default: none)}

\item{missingGenes}{vector of gene IDs/symbols that will be added}

\item{addJitter}{boolean value specifying whether a small amount of jitter should be added to the expression of the supplemented genes (default: FALSE)}

\item{verbose}{specify whether messages should be printed (default: FALSE)}
}
\value{
matrix with added values
}
\description{
add missing genes to expression matrix. Gene expression for these genes in a particular sample will be filled with the average of gene expression across all genes in that sample
}
\author{
Robert Piskol
}
