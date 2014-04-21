# Resample survey data
#' \code{wide_by_site_rarefy} to resample quadrats from survey data
#' 
#' @author Daijiang Li
#' 
#' @param df data frame of vegtation data. Should have site, quad, sp columns. 
#' It will call wide_by_site function inside.
#' @param n times of ramdom sampling from quadrats of 2000s data for each site.
#' @param q how many quadrats to sample out for 2000s data for all sites.
#' @export
#' @return a list of lists. A list of sites, each sites then has n ramdon rampled site by sp data frames 
#' (formed a list for that site).


wide_by_site_rarefy=function(df, n=1000, q = 20){ 
  df1=wide_by_site(df)
  aa=llply(df1, function (x) {
    rlply(.n = n,  x[sample(dim(x)[1], q),])
    }) 
  for (i in 1: length(aa)){
    names(aa[[i]])=rep(names(aa)[i], length(aa[[i]]))
  }
  aa
}