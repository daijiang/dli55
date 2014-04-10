\name{wide_by_site_rarefy}
\alias{wide_by_site_rarefy}
\title{\code{resample} to resample quadrats from survey data}
\usage{
wide_by_site_rarefy(df, n = 1000, q = 20)
}
\arguments{
  \item{df}{data frame of vegtation data. Should have site,
  quad, sp columns. It will call wide_by_site function
  inside.}

  \item{n}{number of ramdom sampling from quadrats of 2000s
  data for each site.}

  \item{q}{how many quadrats to sample out for 2000s data
  for all sites.}
}
\description{
\code{resample} to resample quadrats from survey data
}
\author{
Daijiang Li
}

