#' @useDynLib dli55
#' @importFrom Rcpp sourceCpp
NULL

#' @importFrom ape read.tree write.tree drop.tip compute.brlen vcv.phylo vcv is.rooted 
#' @importFrom picante unifrac phylosor prune.sample 
#' @importFrom graphics hist par rect strwidth text
#' @importFrom dplyr %>% select filter arrange rename mutate left_join
NULL

#' Randomization tests
#' 
#' to perform randomization test between two vectors.
#' 
#' @author Daijiang Li
#' 
#' @param x a vector
#' @param y a vector
#' @param n the number of randomization, default is 1000.
#' @export
#' @return  a data frame with mean, rank, p.value, etc.
rand_test = function(x, y, n = 1000) {
    delta.obs = mean(x, na.rm = T) - mean(y, na.rm = T)
    # cat(delta.obs, '\n')
    pool = c(x, y)
    delta.null = plyr::rdply(.n = n, function() {
        x1.index = sample(length(pool), length(x), replace = F)
        x1 = pool[x1.index]
        y1 = pool[-x1.index]
        data.frame(delta.null = mean(x1, na.rm = T) - mean(y1, na.rm = T))
    })$delta.null
    obs.rank = rank(c(delta.obs, delta.null))[1]
    data.frame(mean_x = mean(x, na.rm = T), mean_y = mean(y, na.rm = T), 
               delta = delta.obs, rank = obs.rank, n = n, p.value = ifelse(obs.rank/(n + 
        1) < 0.5, obs.rank/(n + 1), 1 - obs.rank/(n + 1)))
}

# change a distance matrix or a matrix to a data frame

#' \code{dist_to_df} convert a distance class object or a matrix to a data frame.
#' 
#' @author Daijiang Li
#' 
#' @param x a object of class 'dist' or a matrix.
#' @export
#' @return  a data frame with three columns: site1, site2, distance.
dist_to_df = function(x) {
    mat = as.matrix(x)
    df = reshape2::melt(mat)
    df$Var1 = as.character(df$Var1)
    df$Var2 = as.character(df$Var2)
    df = as.data.frame(t(combn(colnames(mat), 2))) %>% rename(Var1 = V1, Var2 = V2) %>%
      dplyr::left_join(df, by = c("Var1", "Var2"))
    names(df) = c("site1", "site2", "distance")
    rownames(df) = 1:nrow(df)
    df
}

# function to add a column as row names, and remove it from columns
#' \code{var_to_rownames} convert a column of a data frame to be row name.
#' 
#' @author Daijiang Li
#' 
#' @param df a data frame.
#' @param var the column you want to put as row name, the name is quoted.
#' @export
#' @return  a data frame.
var_to_rownames = function(df, var = "site") {
    df = as.data.frame(df)
    row.names(df) = df[, var]
    df[, var] = NULL
    df
}

# function to remove sp that not observed at any site (site by sp matrix)
#' \code{rm_sp_noobs} remove species that not observed in any site.
#' 
#' @author Daijiang Li
#' 
#' @param df a data frame in wide form, i.e. site by species data frame, with site names as row name.
#' @export
#' @return  a data frame.
rm_sp_noobs = function(df) {
    if (any(colSums(df) == 0)) {
        df = df[, -which(colSums(df) == 0)]
    } else {
        df = df
    }
    df
}

# function to remove site that has no observations (site by sp matrix)
#' \code{rm_site_noobs} remove site that has no obsrevations of any species.
#' 
#' @author Daijiang Li
#' 
#' @param df a data frame in wide form, i.e. site by species data frame, with site names as row name.
#' @export
#' @return  a data frame.
rm_site_noobs = function(df) {
    if (any(rowSums(df) == 0)) {
        df = df[-which(rowSums(df) == 0), ]
    } else {
        df = df
    }
    df
}

# logit transformation
#' \code{logit_tran} to conduct logit-transformation.
#' 
#' @author Daijiang Li
#' 
#' @param x a vector of proportions (either in form of 0.20 or 20).
#' @param add_num the number to add for 0s, the dafult is 0.01.
#' @export
#' @return  a vector has been logit transformed.
logit_tran = function(x, add_num = 0.01) {
    if (any(x > 1)) 
        x = x/100  # convert to proportion
    log((x + add_num)/(1 - x + add_num))
}

# to make the first letter upper case
#' \code{simpleCap} to make the first letter to be upper case.
#' 
#' @author Daijiang Li
#' 
#' @param x a vector of species names.
#' @export
#' @return  a vector.
simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}

# Check plant species names for typos and synonyms
#' \code{check_names} to check plant species names in Wisconsin for typos and synonyms using UWSP herbarium as baseline.
#' 
#' @author Daijiang Li
#' 
#' @param sp.to.check A vector of plant species to check.
#' @param sp.dataset species names and their synonyms from UWSP herbarium website.
#' @examples
#' check_names(sp.to.check = c('Acer rubrum', 'Quercus abla'))
#' @return
#' \itemize{
#'   \item sp.origin the original species names.
#'   \item sp.spell the closest correct species names after check for typos, if any.
#'   \item sp.final the updated species names after check for synonyms, if any.
#'   \item spell.correct logical, are the species names spelled correctly?
#'   \item syc.correct logical, are the species names updated?
#' }
#' @export
check_names <- function(sp.to.check, sp.dataset = sp_syn_uwsp) {
    sp.dataset.list = unique(na.omit(sp.dataset$syn))
    # most close spelling in the dataset
    sp.spell.checked = sapply(sp.to.check, function(x) {
        # igonre xxx.spp xxx sp
        if (any(stringr::str_detect(pattern = ".*sp[0-9]?$|.*spp[0-9]?$|.*sp\\.$", string = x))) {
            x
        } else {
            sp.dataset.list[which.min(adist(x, sp.dataset.list))]
        }
    })
    # sp names used as main in webpages of UWSP herbarium
    sp.syc = sapply(sp.spell.checked, function(x) {
        if (any(stringr::str_detect(pattern = ".*sp[0-9]?$|.*spp[0-9]?$|.*sp\\.$", string = x))) {
            x
        } else {
            y = unique(na.omit(sp.dataset$sp[sp.dataset$syn == x]))
            if (length(y > 1)) {
                y = y[which.min(adist(x, y))]
            }
            y
        }
    })
    data.frame(sp.origin = sp.to.check, sp.spell = sp.spell.checked, 
               sp.final = sp.syc, spell.correct = sp.to.check == sp.spell.checked, 
               syc.correct = sp.spell.checked == sp.syc)
}

#' match taxa names with the Open Tree of Life
#' 
#' based on rotl::tnrs_match_names(), which only allow <= 250 species each time. This function will try to do batch match
#' 
#' @param taxa a vector of species names to match
#' @param n_per_req number of species to match per requery, must be <= 250
#' @return if the length of taxa is smaller than n_per_req, then a data frame will be returned. otherwise, a list
#' @export
#' 
tnrs_match_names_2 = function(taxa, n_per_req = 20, ...) {
    n = length(taxa)
    stopifnot(n_per_req <= 250)
    if (n < n_per_req) 
        return(rotl::tnrs_match_names(taxa, ...))
    x = data.frame(sp = taxa, nitem = c(rep(1:floor(n/n_per_req), each = n_per_req), 
                                        rep(ceiling(n/n_per_req), n - (n_per_req * floor(n/n_per_req)))), 
        stringsAsFactors = FALSE)
    out = vector("list", max(x$nitem))
    for (i in 1:(max(x$nitem))) {
        cat(i, " of ", max(x$nitem), "\n")
        out[[i]] = try(rotl::tnrs_match_names(x$sp[x$nitem == i], ...))
    }
    out
}
