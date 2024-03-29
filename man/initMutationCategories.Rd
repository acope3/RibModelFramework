% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{initMutationCategories}
\alias{initMutationCategories}
\title{initMutationCategories}
\arguments{
\item{files}{list of files containing starting values. Number of files should equal the number of categories.}

\item{numCategories}{number of mutation categories (should be less than or equal to number of mixtures)}

\item{fix}{Can use this parameter to fix mutation at current values (won't change over course of MCMC run)}
}
\description{
Initialize values for mutation CSP. File should be of comma-separated with header. Three columns should be of order Amino_acid,Codon,Value
}
