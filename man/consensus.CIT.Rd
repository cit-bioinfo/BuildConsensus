% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defineConsensus.R
\name{consensus.CIT}
\alias{consensus.CIT}
\title{Build a consensus classification using the CIT method}
\usage{
consensus.CIT(annot.class, outdir, CohenKappa.cut = seq(0.1, 0.7, 0.1))
}
\arguments{
\item{annot.class}{Dataframe of samples annotated according to the several 
classification systems to compare.}

\item{outdir}{Path to the directory where to store plots and results. 
A new directory will be created if the supplied path does not exist.}

\item{CohenKappa.cut}{Value of Cohen's Kappa to select significant associations.}
}
\value{
A list of length 3 containing the Cohen's Kappa measures for any pair 
of subtypes (CohenKappaMat), the corresponding network generated (graph), 
the cluster assignments (cl)
}
\description{
Compare and cluster various classification systems using Cohen's Kappa.
}
\examples{
\dontrun{
 #-- Using bladder cancer classes as example
 library(BuildConsensus)
 data(blca_class)
 
 cit_res <- consensus.CIT(blca_class, outdir = "CIT_res")
}
}
\author{
Aurélie Kamoun
}
\keyword{methods}
