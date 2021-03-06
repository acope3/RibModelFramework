% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getSynthesisRatePosteriorVarianceForGene}
\alias{getSynthesisRatePosteriorVarianceForGene}
\title{getSynthesisRatePosteriorVarianceForGene}
\arguments{
\item{samples}{number of samples over which to calculate variance}

\item{geneIndex}{corresponding index of gene in genome for which synthesis rate variance will be calculated. Should be a number between 1 and length(genome)}

\item{unbiased}{Should calculate variance using unbiased (N-1) or biased (N) correction}

\item{log_scale}{Calculate variance on log scale}
}
\value{
posterior mean synthesis rate for gene
}
\description{
Method of Parameter class (access via parameter$<function name>, where parameter is an object initialized by initializeParameterObject). Get synthesis rate variance for a gene
}
