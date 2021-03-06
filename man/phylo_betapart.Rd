% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd.R
\name{phylo_betapart}
\alias{phylo_betapart}
\title{Phylogenetic beta diversity partition}
\usage{
phylo_betapart(comm = dat_1, tree)
}
\arguments{
\item{comm}{a site by species data frame, site names as row names}

\item{tree}{a phylogeny of class "phylo"}
}
\value{
a list of four: pairwise beta of jaccard and sorense; and multisite beta of jaccard and sorensen.
Pairwise beta has three distance matrix. For jaccard, phylo.beta.jtu is the turnover-fraction of Jaccard, phylo.beta.jne is the nestedness-fraction.
For sorensen, phylo.beta.sim is the turnover part measured as Simpson derived pairwise dissimilarity, phylo.beta.sne is the nestedness-fraction.
Similarly for the multisite version.
}
\description{
Calculate Phylogenetic beta diversity and its partition, adapted from betapart::phylo.belt.xx(). Since the pairwise and multisie version
share the same core computation process, it makes more sense to return both. Then we can choose which one to use.
}
