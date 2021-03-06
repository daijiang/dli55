% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pd.R
\name{unifrac2}
\alias{unifrac2}
\title{unifrac}
\usage{
unifrac2(comm, tree, comm_long)
}
\arguments{
\item{comm}{a site by sp data frame, row names are site names}

\item{tree}{a phylogeny with 'phylo' class}

\item{comm_long}{a long format of comm, can be missing}
}
\value{
a site by site distance object
}
\description{
calculate unifrac of pairwise site. This is based on picante::unifrac, but with phylocomr::ph_pd to calculate pd, which can improve speed dramatically.
}
