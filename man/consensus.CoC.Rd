\name{consensus.CoC}
\alias{consensus.CoC}
\title{consensus.CoC}
\description{Compare and cluster various classification systems using a "Cluster of Cluster" approach.}
\usage{consensus.CoC(annot.class, outdir, maxK=8)}
\arguments{
\item{annot.class}{Dataframe of samples annotated according to the several classification systems to compare.}
\item{outdir}{Path to the directory where to store plots and results. A new directory will be created if the supplied path does not exist.}
\item{maxK}{Integer value. maximum cluster number to evaluate.}
}
\value{A list of length maxK. Each element is a list containing consensusMatrix (numerical matrix), consensusTree (hclust), consensusClass (consensus class asssignments)
}
\author{Aurelie Kamoun}
\keyword{methods}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
