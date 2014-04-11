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
wide_by_site=function(df){ 
  aa=dlply(df, .(site), function(x) dcast(x, quad~sp,value.var="count",
                                          fun.aggregate = length)) 
  aa=llply(aa, function(x) {
    rownames(x)=x$quad
    x$quad=NULL
    return(x)
    }) 
  aa
}
