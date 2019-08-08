
# BuildConsensus - Build a Consensus Classification from Multiple Classification Systems

## Overview

This package implements methods for building a consensus molecular classification using existing classification systems. It implements three methods: Cluster of Clusters (CoC), Markov Cluster Algorithm (MCL) and CIT, which uses Cohen's Kappa statistic.

The functions presented were used to build the Consensus Classification of Bladder Cancer [manuscript submitted], but can be used to build a consensus classification in many other settings.

## Installation

For reasons we that will be made clear below, we highly advise using a Debian-based Linux distribution to run this package. (*e.g.* Debian, Ubuntu, Linux Mint)

### MCL

To run the `consensus.MCL` function, you will need to install the MCL library. This can be done in Debian-based Linux distributions using the package manager:

```{bash}
sudo apt-get install mcl
```
This will enable the MCL R interface used by this package. For more information, please see: https://micans.org/mcl/

## GLPK for igraph

The consensus.CIT function uses the GLPK library. This library used to be packaged together with R's igraph library, but now requires separate installation.

If you don't have GLPK installed, please see the GLPK documentation. For Debian-based Linux, you can install GLPK using:

```{bash}
sudo apt-get install libglpk-dev glpk-utils
```
If you already have the igraph library installed in R, please remove it using `remove.packages("igraph")` and re-install `install.packages("igraph")` for the package to recompile.

## Package installation

The package can be installed directly from Github using `devtools` or `remotes` in R.

```{r}
devtools::install_github("cit-bioinfo/BuildConsensus", dependencies = TRUE, build_vignettes = TRUE)
```



