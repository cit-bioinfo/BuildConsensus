\name{consensus.MCL}
\alias{consensus.MCL}
\title{consensus.MCL}
\description{Compare and cluster various classification systems applying the same method as in Guinney et al (see reference).}
\usage{consensus.MCL(annot.class, I.values=3.2, outdir, resamp=.8, n.iter=1000, pval.cut=0.001)}
\arguments{
\item{annot.class}{Dataframe of samples annotated according to the several classification systems to compare.}
\item{I.values}{A vector of inflation factors values to use with MCL algorithm. The function running time linearly depends on the given number of values (about 3 minutes for each inflation factor value)}
\item{outdir}{Path to the directory where to store plots and results. A new directory will be created if the supplied path does not exist.}
\item{resamp}{Value between 0 and 1. Proportion of samples to use for MCL bootstrap iterations. Default is 0.8.}
\item{n.iter}{Number of MCL bootstrap iterations for each inflation factor value. Default is 1000.}
\item{pval.cut}{Cut-off for adjusted P-value to select significant network edges (Fisher test). Default is 0.001 }
}
\value{A list with the same length of I.values numerical vector. Each element is a list containing the cluster assignments (cl), the co-classification matrix (consensusMat), a subtype stability measure (subtype.stab), the mean weighted average silhouette width (wsil.mean)
}
\author{Aurelie Kamoun}
\keyword{methods}
\references{1.Guinney, J. et al. The consensus molecular subtypes of colorectal cancer. Nat Med advance online publication, (2015).}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
