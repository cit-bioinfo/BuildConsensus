\name{consensus.CIT}
\alias{consensus.CIT}
\title{consensus.CIT}
\description{Compare and cluster various classification systems using Cohen's Kappa.}
\usage{consensus.CIT(annot.class, outdir, CohenKappa.cut=.2)}
\arguments{
\item{annot.class}{Dataframe of samples annotated according to the several classification systems to compare.}
\item{outdir}{Path to the directory where to store plots and results. A new directory will be created if the supplied path does not exist.}
\item{CohenKappa.cut}{Value of Cohen's Kappa to select significant associations.}
}
\value{A list of length 3 containing the Cohen's Kappa measures for any pair of subtypes (CohenKappaMat), the corresponding network generated (graph), the cluster assignments (cl)
}
\author{Aurelie Kamoun}
\keyword{methods}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
