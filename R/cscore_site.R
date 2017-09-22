# Calculate C-score for each site
#' \code{cscore_site} to calculate C-score for each site.
#' 
#' @author Daijiang Li
#' 
#' @param df data frame of vegtation data. Should have site, quad, sp columns. 
#' This function will call wide_by_site function inside.
#' @param method methods for building null models. Check oecosium from vegan package for more options. 'c0'
#' for fixed species incidence and equalprobalility for site (site species richness can vary); 'quasiswap'
#' for fixed sp - fixed site,  suppose to be faster; 'tswap' also for fixed-fixed but slower. Default is 'c0'.
#' @param nsimul how many null communities to use for simulation? Default is 1000.
#' @param rarefy Do you want to rarefy sites or not? Default is FALSE.
#' @param n used only if rarefy=TRUE. n times of ramdom sampling from quadrats of 2000s data for each site.
#' @param q used only if rarefy=TRUE. How many quadrats to sample out for 2000s data for all sites.
#' @export
#' @return If no rarefy, it returns a dataframe. It has 'z', 'means_sim', 'pval', 'cs_obs', 'site' columns.
#' If rarefy=TRUE, it returns a list: detail and summary for all sites.
#' @details If no rarefy, this function uses wide_by_site() function to change long veg table into wide table for each site (returned
#' a list), then for each site of the list, calculate the standardized effect size z = (X_obs - X_meanSim) / sd_sim 
#' using null models. 
#' If rarefy=TRUE, it calls wide_by_site_rarefy() function to change long veg table into wide table for each site
#' n times
cscore_site = function(df, method = "c0", nsimul = 1000, rarefy = FALSE, n = 1000, q = 20) {
    cscore_bipartite = function(df) bipartite::C.score(df, normalise = F, na.rm = T)  # bipartite package,  df here is quad by sp matrix
    
    cs_sim = function(df) {
        # vegan package, df here is quad by sp matrix
        vegan::oecosimu(df, nestfun = cscore_bipartite, method = method, nsimul = nsimul, 
                        alternative = "two.sided")[[2]][-c(4, 5, 7)]
        # return z, means of sim, pval, observed cscore
    }
    
    if (rarefy) {
        ll = wide_by_site_rarefy(df, n = n, q = q)
        detail = plyr::ldply(ll, function(df) {
            # deal with all random sampling of one site
            aa = plyr::llply(df, cs_sim)
            bb = as.data.frame(matrix(unlist(aa), ncol = 4, byrow = T), stringsAsFactors = F)
            bb$site = names(aa)
            names(bb) = c("z", "means_sim", "pval", "cs_obs", "site")
            bb
        })  # deal with all sites
        summary = plyr::ddply(detail, .(site), summarize, z = mean(z))  # get mean z for each site
        return(list(detail = detail, summary = summary))
        
    } else {
        
        aa = plyr::llply(wide_by_site(df), cs_sim)
        bb = as.data.frame(matrix(unlist(aa), ncol = 4, byrow = T), stringsAsFactors = F)
        bb$site = names(aa)
        names(bb) = c("z", "means_sim", "pval", "cs_obs", "site")
        return(bb)
    }
}

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


wide_by_site_rarefy = function(df, n = 1000, q = 20) {
    df1 = wide_by_site(df)
    aa = plyr::llply(df1, function(x) {
        plyr::rlply(.n = n, x[sample(dim(x)[1], q), ])
    })
    for (i in 1:length(aa)) {
        names(aa[[i]]) = rep(names(aa)[i], length(aa[[i]]))
    }
    aa
}

#' Change long table into wide table for each site
#'
#' Produce wide quad by species data frame for each site from long table in a list. 
#' Each element of the list is a quadrat by species data frame.
#' 
#' \code{wide_by_site} to change from long table to wide table for each site
#' 
#' @author Daijiang Li
#' @param df data frame of vegtation data. Should have site, quad, sp columns.
#' @return list. 
#' @export

# split sites into list, each element is a site by species data frame for one site
wide_by_site = function(df) {
    aa = plyr::dlply(df, .(site), function(x) {
        reshape2::dcast(x, quad ~ sp, value.var = "count", fun.aggregate = length)
    })
    
    aa = plyr::llply(aa, function(x) {
        rownames(x) = x$quad
        x$quad = NULL
        return(x)
    })
    aa
}
