# Calculate C-score for each site
#' \code{cscore_site} to calculate C-score for each site.
#' 
#' @author Daijiang Li
#' 
#' @param df data frame of vegtation data. Should have site, quad, sp columns. 
#' This function will call wide_by_site function inside.
#' @param method methods for building null models. Check oecosium from vegan package for more options. "c0"
#' for fixed species incidence and equalprobalility for site (site species richness can vary); "quasiswap"
#' for fixed sp - fixed site,  suppose to be faster; "tswap" also for fixed-fixed but slower. Default is "c0".
#' @param nsimul how many null communities to use for simulation? Default is 1000.
#' @param rarefy Do you want to rarefy sites or not? Default is FALSE.
#' @param n used only if rarefy=TRUE. n times of ramdom sampling from quadrats of 2000s data for each site.
#' @param q used only if rarefy=TRUE. How many quadrats to sample out for 2000s data for all sites.
#' @export
#' @return If no rarefy, it returns a dataframe. It has "z", "means_sim", "pval", "cs_obs", "site" columns.
#' If rarefy=TRUE, it returns a list: detail and summary for all sites.
#' @details If no rarefy, this function uses wide_by_site() function to change long veg table into wide table for each site (returned
#' a list), then for each site of the list, calculate the standardized effect size z = (X_obs - X_meanSim) / sd_sim 
#' using null models. 
#' If rarefy=TRUE, it calls wide_by_site_rarefy() function to change long veg table into wide table for each site
#' n times
cscore_site=function(df, method="c0", nsimul=1000, rarefy=FALSE, n=1000,q=20){
  cscore_bipartite=function(df) C.score(df, normalise=F,na.rm=T)  # bipartite package,  df here is quad by sp matrix
  
  cs_sim=function(df){# vegan package, df here is quad by sp matrix
    oecosimu(df,nestfun=cscore_bipartite,method=method,nsimul=nsimul, 
             alternative="two.sided")[[2]][-c(4,5,7)]
    # return z, means of sim, pval, observed cscore
  }
 
  if(rarefy){
    ll=wide_by_site_rarefy(df, n=n, q=q)
    detail=ldply(ll, function(df){ # deal with all random sampling of one site
     aa=llply(df, cs_sim)
     bb=as.data.frame(matrix(unlist(aa), ncol=4, byrow=T), stringsAsFactors=F)
     bb$site=names(aa)
     names(bb)=c("z", "means_sim", "pval", "cs_obs", "site")
     bb
  }) # deal with all sites
    summary=ddply(detail, .(site), summarize, z=mean(z)) # get mean z for each site
    return(list(detail=detail, summary=summary))
    
  } else {
    
    aa=llply(wide_by_site(df), cs_sim)
    bb=as.data.frame(matrix(unlist(aa), ncol=4, byrow=T), stringsAsFactors=F)
    bb$site=names(aa)
    names(bb)=c("z", "means_sim", "pval", "cs_obs", "site")
    return(bb)
  }
}