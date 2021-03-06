% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRAPE_Basic_Functions.R
\name{makePairwiseOrderNames}
\alias{makePairwiseOrderNames}
\title{Make template names from gene names}
\usage{
makePairwiseOrderNames(path_genes)
}
\arguments{
\item{path_genes}{A vector of pathway gene names}
}
\value{
Names for the pairwise representation, of the form "gA<gB"
}
\description{
Takes in vector of pathway gene names, returns names corresponding to the pairwise binary representation
}
\examples{
path_genes <- c("gene_A","gene_B","gene_C","gene_D")
makePairwiseOrderNames(path_genes)
}
