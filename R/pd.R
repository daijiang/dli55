#' Faith's PD
#' 
#' Calculate faith's pd, this is a wrapper of phylocomr::ph_pd and picante::pd
#' 
#' @param comm a site by species data frame, site names as row names, only works for presence/absence data
#' @param tree a phylogeny of class "phylo"
#' @param include.root whether include root in picante::pd; phylocomr::ph_pd always include root
#' @param comm_long long format of comm, 3-columns: site, freq, sp; optional.
#' @return a data frame of PD for each site
#' @export
#'
pd2 = function(comm, tree, include.root = TRUE, comm_long){
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  class(tree) = "phylo"
  
  if (include.root) {
    # Make sure tree is rooted if needed
    if (!is.rooted(tree)) {
      stop("Rooted tree required to calculate PD with include.root=TRUE argument")
    }
    tree <- picante::node.age(tree)
    if(missing(comm_long)){
      comm_long = tibble::rownames_to_column(as.data.frame(comm), "site") %>% 
        tidyr::gather("sp", "freq", -site) %>% filter(freq > 0) %>% 
        arrange(site, sp) %>% select(site, freq, sp)
    }
    pdcomm = try(phylocomr::ph_pd(sample = comm_long, phylo = tree) %>% 
                   rename(pd.root = pd, site = sample)) # phylocom cannot ignore root
    
    if ("try-error" %in% class(pdcomm)) {
      cat("Phylocom has trouble with this phlyogney, switch to picante", "\n")
      pdcomm = picante::pd(comm, tree, include.root = TRUE) %>% 
        tibble::rownames_to_column("site") %>% rename(pd.root = PD)
    } 
    
    pdcomm = tibble::data_frame(site = row.names(comm)) %>% 
      dplyr::left_join(pdcomm, by = "site") # make sure the same 
  } else {
    pdcomm = tibble::data_frame(site = row.names(comm),
                                pd.uroot = PhyloMeasures::pd.query(tree, comm, standardize = FALSE)
    )
  }
 
  return(pdcomm)
}

#' calculate MPD and VPD (mean and variance of pairwise distance)
#' 
mvpd <- function(samp, dis, abundance.weighted=FALSE){
  N <- dim(samp)[1]
  mpd = numeric(N); vpd = numeric(N)
  for (i in 1:N) {
    # cat("row ", i)
    sppInSample <- names(samp[i, samp[i, ] > 0])
    if (length(sppInSample) > 1) {
      sample.dis <- dis[sppInSample, sppInSample]
      if (abundance.weighted) {
        sample.weights <- t(as.matrix(samp[i,sppInSample,drop=FALSE])) %*% as.matrix(samp[i,sppInSample,drop=FALSE])
        wm = weighted.mean(sample.dis,sample.weights)
        mpd[i] <- wm
        vpd[i] = sum(sample.weights[lower.tri(sample.weights)] * (sample.dis[lower.tri(sample.dis)] - wm)^2)
      }
      else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
        vpd[i] <- var(sample.dis[lower.tri(sample.dis)])
      }
    }
    else{
      mpd[i] = NA
      vpd[i] = NA
    }
  }
  data.frame(site = row.names(samp), mpd = mpd, vpd = vpd, stringsAsFactors = F)
}

#' calculate alpha phylogenetic diversity
#' 
#' A function to calculate a bunch of phylo diversity: Faith's PD, MPD, VPD (variance of pairwise distance), MNTD, PSV/PSE.
#' Results of MPD/MNTD from PhyloMeasures are equal with those from Phylocom/Picante with abundance.weight = FALSE
#' Results of PD from PhyloMeasures are unrooted; PD from Phylocom are rooted.
#' 
#' @param samp_wide wide version of samp_long, row.names are sites, colnames are species
#' @param tree a phylogeny with class of 'phylo'
#' @param samp_long a 3-column data frame, site, freq, sp
#' @param null.model whether to run null models for MPD and MNTD?
#' @param null.model.pd.root whether to run null models for rooted PD?
#' @param null.type.phylomeasures if null.model is TRUE, which null model to use for PhyloMeasure? 
#' @param null.type.phylocom if null.model is TRUE, which null model to use for Phylocom? 
#' - 0: phylogeny shuffle. 
#' - 1: maintain site richness and draw from species actually observed in sites. 
#' - 2: maintain site richness and draw from the phylogeny
#' - 3: independent swap. 
#' See ?phylocomr::ph_comstruct for details
#' @param null.type.picante if null.model is TRUE, which null model to use for Picante? Not used yet. 
#' @param n.item the number of randomization
#' @param abund.weight should abundance information used when calculating pd with Phylocom/Picante? Default is FALSE
#' @param verbose do you want to see relevant information?
#' @param vpd to calculate vpd (varance of pairwise distance) or not?
#' @return a data frame
#' @export
#' 
get_pd_alpha = function(samp_wide, tree, samp_long, 
                        null.model.phylomeasures = TRUE, null.model.phylocom = FALSE,
                        null.type.phylomeasures = "uniform", 
                        null.type.phylocom = 0, null.type.picante = "taxa.labels",
                        n.item = 999,
                        abund.weight = FALSE, 
                        verbose = TRUE, vpd = FALSE, ...){
  if(length(class(tree)) > 1 & "phylo" %in% class(tree)) class(tree) = "phylo"
  dist = cophenetic(tree)
  
  row.names(samp_wide) = stringr::str_trim(row.names(samp_wide)) %>% 
    tolower() %>% 
    stringr::str_replace_all(" ", "_") # phylocom will have trouble with space in site names
  
  if(missing(samp_long)){
    samp_long = tibble::rownames_to_column(as.data.frame(samp_wide), "site") %>% 
      tidyr::gather("sp", "freq", -site) %>% filter(freq > 0) %>% 
      arrange(site, sp) %>% select(site, freq, sp)
  }
  
  # faith pd
  # unrooted pd
    faith_pd = tibble::data_frame(site = row.names(samp_wide),
                       pd.uroot = PhyloMeasures::pd.query(tree = tree, matrix = samp_wide, standardize = FALSE)
    )
    if(null.model.phylomeasures){
      faith_pd$pd.uroot.z = PhyloMeasures::pd.query(tree = tree, matrix = samp_wide, 
                                                    standardize = TRUE, 
                                                    null.model = null.type.phylomeasures)
    }
  # rooted pd
    faith_pd2 = pd2(samp_wide, tree, include.root = TRUE) %>% select(site, pd.root)
    if(null.model.phylocom) {
      message("no analytical null model for rooted pd calculated with phylocom/picante yet.")
      pd_z = purrr::map(1:n.item, function(x){
        set.seed(x)
        if(verbose) cat("null model of pd", x, "\t")
        x2 = pd2(samp_wide, picante::tipShuffle(tree), include.root = TRUE) %>% select(site, pd.root)
        x2
      })
      names(pd_z) = paste0("rep_", 1:n.item)
      pd_z = dplyr::bind_rows(pd_z, .id = "repli")
      pd_z = dplyr::group_by(pd_z, site) %>% dplyr::summarise(ave_pd = mean(pd.root, na.rm = TRUE),
                                                       sd_pd = sd(pd.root, na.rm = TRUE))
      faith_pd2 = dplyr::left_join(faith_pd2, pd_z, by = "site") %>% 
        dplyr::mutate(pd.root.z = (pd.root - ave_pd)/sd_pd) %>% 
        dplyr::select(site, pd.root, pd.root.z)
    }

    out = left_join(faith_pd, faith_pd2, by = "site")
  # mpd and mntd
  if(abund.weight){ # PhyloMeasures cannot weight by abundance, use Phylocom or Picante
    if(null.model.phylocom){
      mpd_mntd = try(phylocomr::ph_comstruct(sample = samp_long, phylo = tree, 
                                         null_model = null.type.phylocom, 
                                         randomizations = n.item, 
                                         abundance = abund.weight))
      if("try-error" %in% class(mpd_mntd)){
        if(verbose) cat("Phylocom has trouble with this phlyogney, switch to picante", "\n")
        mpd_mntd = mvpd(samp_wide, dist, abundance.weighted = TRUE) %>% left_join(
          tibble::data_frame(site = row.names(samp_wide),
                             mntd = picante::mntd(samp_wide, dist, abundance.weighted = TRUE)
          ), by = "site")
        cat("No null model for mpd/mntd with picante yet", "\n")
      } else {
        mpd_mntd = select(mpd_mntd, plot, mpd, mntd, nri, nti) %>% 
          rename(site = plot) %>% 
          mutate(mpd.z = -1 * nri, mntd.z = -1 * nti) %>% 
          select(-nri, -nti)
      }
    } else {
      mpd_mntd = phylocomr::ph_comstruct(sample = samp_long, phylo = tree, 
                                         null_model = null.type.phylocom, 
                                         randomizations = 0, # phylocom has no turn off of null model
                                         abundance = abund.weight)
      if("try-error" %in% class(mpd_mntd)){
        if(verbose) cat("Phylocom has trouble with this phlyogney, switch to picante", "\n")
        mpd_mntd = mvpd(samp_wide, dist, abundance.weighted = TRUE) %>% left_join(
          tibble::data_frame(site = row.names(samp_wide),
                             mntd = picante::mntd(samp_wide, dist, abundance.weighted = TRUE)
          ), by = "site")
      } else {
        mpd_mntd = select(mpd_mntd, plot, mpd, mntd) %>% rename(site = plot)
      }
    }
    out = left_join(out, mpd_mntd, by = "site")
  } else { # no weight on abundance, i.e. presence/absence, use PhyloMeasures, 
    # results equal with phylocom/picante with abundace.weight = F
    mpd_s = tibble::data_frame(site = row.names(samp_wide),
                               mpd = PhyloMeasures::mpd.query(tree, samp_wide, standardize = FALSE)
    )
    if(null.model.phylomeasures){
      mpd_s$mpd.z = PhyloMeasures::mpd.query(tree = tree, matrix = samp_wide, 
                                             standardize = TRUE, 
                                             null.model = null.type.phylomeasures)
    }
    
    mntd_s = tibble::data_frame(site = row.names(samp_wide),
                                mntd = PhyloMeasures::mntd.query(tree, samp_wide, standardize = FALSE)
    )
    if(null.model.phylomeasures){
      mntd_s$mntd.z = PhyloMeasures::mntd.query(tree = tree, matrix = samp_wide, 
                                                standardize = TRUE, 
                                                null.model = null.type.phylomeasures)
    }
    out = left_join(out, mpd_s, by = "site") %>% left_join(mntd_s, by = "site")
  }

  # variance of pairwise distance
  if(vpd & !"vpd" %in% names(out)) {
    mvpd_s = mvpd(samp_wide, dist, abundance.weighted = abund.weight)
    out = left_join(out, select(mvpd_s, -mpd), by = "site")
  }
  
  # PSV
  # psr_s = psr(samp, tree, compute.var = FALSE) %>% mutate(site = row.names(samp)) %>% select(-SR)
  psv_s = picante::psv(samp_wide, tree, compute.var = FALSE) %>% 
    mutate(site = row.names(samp_wide)) %>% select(-SR) %>% rename(psv = PSVs)
  out = left_join(out, psv_s, by = "site")
  if(abund.weight){
    pse_s = picante::pse(samp_wide, tree) %>% 
      mutate(site = row.names(samp_wide)) %>% select(-SR) %>% rename(pse = PSEs)
    out = left_join(out, pse_s, by = "site")
  }
  
  # species richness
  out = vegan::specnumber(samp_wide) %>% as.data.frame() %>% setNames("sr") %>% 
    tibble::rownames_to_column("site") %>% 
    tibble::rowid_to_column() %>% 
    left_join(out, by = "site") %>% 
    select(site, rowid, sr, starts_with("pd"), mpd, mntd, psv, everything()) %>% 
    dplyr::tbl_df()
  out
}

# beta phylo diversity ----

#' Predicted PCD with species pool
#' 
#' expected PCD from one or two communities (depends on the species pool)
#' 
#' @param comm_old a site by species dataframe or matrix, with sites as rows and species as columns.
#' @param comm_new an optional second site by species data frame. It should have the same number of rows as comm_old.
#' @param tree the phylogeny for communities.
#' @param reps number of random draws, default is 1000 times.
#' @param cpp whether to use loops written with c++, default is TRUE
#' @return a list with species richness of the pool, expected PSV, PSV of the pool, and unique number of species richness across sites.
#' @export
#' 
pcd_pred = function(comm_old, comm_new = NULL, tree, reps = 10^3, cpp = TRUE) {
    # Make comm matrix a pa matrix
    comm_old[comm_old > 0] = 1
    if (!is.null(comm_new)) 
        comm_new[comm_new > 0] = 1
    
    if (is.null(comm_new)) {
        sp_pool = colnames(comm_old)
    } else {
        sp_pool = unique(c(colnames(comm_old), colnames(comm_new)))
    }
    
    # convert trees to VCV format
    if (is(tree)[1] == "phylo") {
        if (is.null(tree$edge.length)) {
            # If phylo has no given branch lengths
            tree = compute.brlen(tree, 1)
        }
        tree = drop.tip(tree, tip = tree$tip.label[!tree$tip.label %in% sp_pool])
        V = vcv.phylo(tree, corr = TRUE)
        comm_old = comm_old[, tree$tip.label[tree$tip.label %in% colnames(comm_old)]]
        if (!is.null(comm_new)) {
            comm_new = comm_new[, tree$tip.label[tree$tip.label %in% colnames(comm_new)]]
        }
        if (is.null(comm_new)) {
            sp_pool = colnames(comm_old)
        } else {
            sp_pool = unique(c(colnames(comm_old), colnames(comm_new)))
        }
    } else {
        V = tree
        V = V[sp_pool, sp_pool]
    }
    
    
    # m=number of communities; n=number of species; nsr=maximum sp richness value across all communities
    n = length(sp_pool)
    if (is.null(comm_new)) {
        nsr = unique(rowSums(comm_old))
    } else {
        nsr = unique(c(rowSums(comm_old), rowSums(comm_new)))
    }
    
    if (cpp) {
        SSii = predict_cpp(n = n, nsr, reps = reps, V = V)
    } else {
        SSii = vector("numeric", length(nsr))
        n1 = 2  # the number of n1 does not matter
        for (n2 in 1:length(nsr)) {
            temp = array(0, reps)
            for (t in 1:reps) {
                rp = sample(n)
                pick1 = rp[1:n1]
                
                rp = sample(n)
                pick2 = rp[1:nsr[n2]]
                
                C11 = V[pick1, pick1]
                C22 = V[pick2, pick2]
                C12 = V[pick1, pick2]
                
                invC22 = solve(C22)
                # S11 = C11 - C12%*%invC22%*%t(C12)
                S11 = C11 - tcrossprod(C12 %*% invC22, C12)
                SS11 = (n1 * sum(diag(S11)) - sum(S11))/(n1 * (n1 - 1))
                temp[t] = SS11
            }
            SSii[n2] = mean(temp)
        }
    }
    names(SSii) = as.character(nsr)
    
    SCii = 1 - (sum(V) - sum(diag(V)))/(n * (n - 1))
    
    return(list(nsp_pool = n, psv_bar = SSii, psv_pool = SCii, nsr = nsr))
}

#' pairwise site phylogenetic community dissimilarity (PCD) within a community
#' 
#' Calculate pairwise site PCD, users can specify expected values from \code{pcd_pred()}.
#'
#' @param comm site by species data frame or matrix, sites as rows.
#' @param tree a phylogeny for species
#' @param expectation nsp_pool, psv_bar, psv_pool, and nsr calculated from \code{pcd_pred()}.
#' @param cpp whether to use loops written with c++, default is TRUE
#' @param unif_dim the number of cells (nrow * ncol) of the comm, if it is less than unif_dim, then calculate unifrac and phylosor; these functions are very slow.
#' @param verbose do you want to see the progress?
#' @return a list of a variety of pairwise dissimilarities.
#' @export
#' @examples
#' x1 = pcd_pred(comm_old = read.csv('data/li_2015_old.csv', row.names = 1, check.names = F),
#'               comm_new = read.csv('data/li_2015_new.csv', row.names = 1, check.names = F),
#'                  tree = ape::read.tree('data/phy.tre'), reps = 100)
#' pcd2(comm = read.csv('data/li_2015_old.csv', row.names = 1, check.names = F),
#'       tree = ape::read.tree('data/phy.tre'), 
#'       expectation = x1)
pcd2 = function(comm, tree, expectation = NULL, cpp = TRUE, unif_dim = 1000, verbose = TRUE) {
    if (is.null(expectation)) {
        expectation = pcd_pred(comm_old = comm, tree = tree)
    }
    nsp_pool = expectation$nsp_pool
    nsr = expectation$nsr
    SSii = expectation$psv_bar
    SCii = expectation$psv_pool
    
    if (nrow(comm) * ncol(comm) < unif_dim) {
        unif = unifrac2(comm, tree)  # a DISSIMILAR distance matrix
        physor = 1 - picante::phylosor(comm, tree)  # a DISSIMILAR distance matrix
    } else {
        unif = physor = NA
    }
    
    # Make comm matrix a pa matrix
    comm[comm > 0] = 1
    
    # convert trees to VCV format
    if (is(tree)[1] == "phylo") {
        if (is.null(tree$edge.length)) {
            # If phylo has no given branch lengths
            tree = compute.brlen(tree, 1)
        }
        tree = prune.sample(comm, tree)
        V = vcv.phylo(tree, corr = TRUE)
        comm = comm[, tree$tip.label]
    } else {
        V = tree
        species = colnames(comm)
        preval = colSums(comm)/sum(comm)
        species = species[preval > 0]
        V = V[species, species]
        comm = comm[, colnames(V)]
    }
    
    if (!is.null(SSii) & length(SSii) < length(unique(rowSums(comm)))) {
        stop("The length of PSVbar is less than the unique number of species richness of the community.")
    }
    
    if (cpp) {
        xxx = pcd2_loop(SSii, nsr, SCii, as.matrix(comm), V, nsp_pool, verbose)
        PCD = xxx$PCD
        PCDc = xxx$PCDc
        PCDp = xxx$PCDp
        D_pairwise = xxx$D_pairwise
        dsor_pairwise = xxx$dsor_pairwise
    } else {
        # m=number of communities; n=number of species; nsr=maximum sp rich value across all communities
        m = dim(comm)[1]
        n = dim(comm)[2]
        
        # calculate PCD
        PCD = array(NA, c(m, m))
        PCDc = array(NA, c(m, m))
        PCDp = array(NA, c(m, m))
        D_pairwise = array(NA, c(m, m))
        dsor_pairwise = array(NA, c(m, m))
        for (i in 1:(m - 1)) {
            if(verbose) cat(i, " ")
            for (j in (i + 1):m) {
                # cat('i = ', i, ', j = ', j, '\n')
                pick1 = (1:n)[comm[i, ] == 1]
                pick2 = (1:n)[comm[j, ] == 1]
                
                n1 = length(pick1)
                n2 = length(pick2)
                
                C = V[c(pick1, pick2), c(pick1, pick2)]
                # cat('dim of C: ', dim(C), ' n1 = ', n1, ' n2 = ', n2, '\n')
                C11 = C[1:n1, 1:n1]
                C22 = C[(n1 + 1):(n1 + n2), (n1 + 1):(n1 + n2)]
                C12 = C[1:n1, (n1 + 1):(n1 + n2)]
                if (is.null(dim(C12))) {
                  if (is.null(dim(C22))) {
                    C12 = as.matrix(C12)
                  } else {
                    C12 = t(as.matrix(C12))
                  }
                }
                
                invC11 = solve(C11)
                # S22 = C22 - t(C12)%*%invC11%*%C12
                S22 = C22 - crossprod(C12, invC11) %*% C12
                
                invC22 = solve(C22)
                # S11 = C11 - C12%*%invC22%*%t(C12)
                S11 = C11 - tcrossprod(C12 %*% invC22, C12)
                if (n1 > 1) {
                  SC11 = (n1 * sum(diag(C11)) - sum(C11))/(n1 * (n1 - 1))
                  SS11 = (n1 * sum(diag(S11)) - sum(S11))/(n1 * (n1 - 1))
                } else {
                  SC11 = (n1 * sum(diag(C11)) - sum(C11))/(n1 * n1)
                  SS11 = (n1 * sum(diag(S11)) - sum(S11))/(n1 * n1)
                }
                if (n2 > 1) {
                  SC22 = (n2 * sum(diag(C22)) - sum(C22))/(n2 * (n2 - 1))
                  SS22 = (n2 * sum(diag(S22)) - sum(S22))/(n2 * (n2 - 1))
                } else {
                  SC22 = (n2 * sum(diag(C22)) - sum(C22))/(n2 * n2)
                  SS22 = (n2 * sum(diag(S22)) - sum(S22))/(n2 * n2)
                }
                
                D = (n1 * SS11 + n2 * SS22)/(n1 * SC11 + n2 * SC22)
                
                dsor = 1 - 2 * length(intersect(pick1, pick2))/(n1 + n2)
                
                pred.D = unname((n1 * SSii[as.character(n2)] + n2 * SSii[as.character(n1)])/(n1 * SCii + n2 * SCii))
                pred.dsor = 1 - 2 * n1 * n2/((n1 + n2) * nsp_pool)
                
                PCD[i, j] = D/pred.D
                PCDc[i, j] = dsor/pred.dsor
                PCDp[i, j] = PCD[i, j]/PCDc[i, j]
                D_pairwise[i, j] = D
                dsor_pairwise[i, j] = dsor
            }
        }
    }
    
    colnames(PCD) = rownames(comm)
    rownames(PCD) = rownames(comm)
    colnames(PCDc) = rownames(comm)
    rownames(PCDc) = rownames(comm)
    colnames(PCDp) = rownames(comm)
    rownames(PCDp) = rownames(comm)
    colnames(D_pairwise) = rownames(comm)
    rownames(D_pairwise) = rownames(comm)
    colnames(dsor_pairwise) = rownames(comm)
    rownames(dsor_pairwise) = rownames(comm)
    
    output = list(PCD = as.dist(t(PCD)), PCDc = as.dist(t(PCDc)), 
                  PCDp = as.dist(t(PCDp)), D_pairwise = as.dist(t(D_pairwise)), 
                  dsor_pairwise = as.dist(t(dsor_pairwise)),
                  uni_frac = unif, phy_sor = physor)
    
    # plyr::ldply(output, function(x){ dist_to_df(x) }) %>% rename(id = .id)
    output
}

#' unifrac
#' 
#' calculate unifrac of pairwise site. This is based on picante::unifrac, but with phylocomr::ph_pd to calculate pd, which can improve speed dramatically.
#' 
#' @param comm a site by sp data frame, row names are site names
#' @param tree a phylogeny with 'phylo' class
#' @param comm_long a long format of comm, can be missing
#' @return a site by site distance object
#' @export
#' 
unifrac2 <- function(comm, tree, comm_long) {
    if (is.null(tree$edge.length)) {
        stop("Tree has no branch lengths, cannot compute UniFrac")
    }
    
    if (!is.rooted(tree)) {
        stop("Rooted phylogeny required for UniFrac calculation")
    }
    
    class(tree) = "phylo"
    
    row.names(comm) = stringr::str_trim(row.names(comm)) %>% 
      tolower() %>% 
      stringr::str_replace_all(" ", "_") # phylocom will have trouble with space in site names
    
    if (missing(comm_long)) {
        comm_long = tibble::rownames_to_column(as.data.frame(comm), "site") %>% 
          tidyr::gather("sp", "freq", -site) %>% filter(freq > 0) %>% 
          arrange(site, sp) %>% select(site, freq, sp)
    }
    
    comm <- as.matrix(comm)
    s <- nrow(comm)
    phylodist <- matrix(NA, s, s)
    rownames(phylodist) <- rownames(comm)
    colnames(phylodist) <- rownames(comm)
    
    comm_comb <- matrix(NA, s * (s - 1)/2, ncol(comm))
    colnames(comm_comb) <- colnames(comm)
    
    i <- 1
    for (l in 1:(s - 1)) {
        for (k in (l + 1):s) {
            comm_comb[i, ] <- comm[l, ] + comm[k, ]
            i <- i + 1
        }
    }
    
    pdcomm = try(phylocomr::ph_pd(sample = comm_long, phylo = tree) %>% 
                   rename(PD = pd, site = sample))
    if ("try-error" %in% class(pdcomm)) {
        cat("Phylocom has trouble with this phlyogney, switch to picante", "\n")
        pdcomm = picante::pd(comm, tree, include.root = TRUE)
        pdcomm_comb <- picante::pd(comm_comb, tree)
    } else {
        comm_comb_long = tibble::rownames_to_column(as.data.frame(comm_comb), "site") %>% 
          tidyr::gather("sp", "freq", -site) %>% filter(freq > 0) %>% 
          arrange(site, sp) %>% select(site, freq, sp)
        pdcomm_comb = try(phylocomr::ph_pd(sample = comm_comb_long, phylo = tree) %>%
                            rename(PD = pd, site = sample) %>% arrange(site))
    }
    pdcomm = tibble::data_frame(site = row.names(comm)) %>% left_join(pdcomm, by = "site")  # make sure the same order
    
    i <- 1
    for (l in 1:(s - 1)) {
        pdl <- pdcomm[l, "PD"]
        for (k in (l + 1):s) {
            pdk <- pdcomm[k, "PD"]
            pdcomb <- pdcomm_comb[i, "PD"]
            pdsharedlk <- pdl + pdk - pdcomb
            phylodist[k, l] = purrr::as_vector((pdcomb - pdsharedlk)/pdcomb)
            i <- i + 1
        }
    }
    return(as.dist(phylodist))
}

#' Phylogenetic beta diversity partition
#' 
#' Calculate Phylogenetic beta diversity and its partition, adapted from betapart::phylo.belt.xx()
#' 
#' @param comm a site by species data frame, site names as row names
#' @param tree a phylogeny of class "phylo"
#' @param index.family "jaccard" or "sorensen"
#' @param pairwise calculate pairwise beta diversity or multisite beta diversity
#' @return a list of three distance matrix. For jaccard, phylo.beta.jtu is the turnover-fraction of Jaccard, phylo.beta.jne is the nestedness-fraction.
#' For sorensen, phylo.beta.sim is the turnover part measured as Simpson derived pairwise dissimilarity, phylo.beta.sne is the nestedness-fraction.
#' @export
#'
phylo_betapart = function(comm, tree, index.family = "jaccard", pairwise = TRUE){
  # adapted from betaprt::phylo.beta
  index.family <- match.arg(index.family, c("jaccard", "sorensen"))
  
  if (!is.matrix(comm)) {
    comm <- as.matrix(comm)
  }
  
  if(nrow(comm)<2) stop("Computing dissimilairty requires at least 2 communities", call. = TRUE)
  if (!is.numeric(comm)) stop("The data in 'comm' is not numeric.", call. = TRUE)
  xvals <- unique(as.vector(comm))
  if (any(!is.element(xvals, c(0, 1)))){
    warning("The community matrix contains values other than 0 and 1: trying to convert data to presence/absence.", call. = TRUE)
    comm[comm > 1] = 1
  }
    
  if (class(tree)!="phylo")
    stop("### invalid tree's format: \"phylo\" format required ###\n\n", call. = TRUE)
  if(any(!(colnames(comm)%in%tree$tip)))
    warning("At least one species in community matrix is not included in the tree" , call. = TRUE)
  
  ############ Paired matrix to distance matrix conversion (utility function) #######################
  
  dist.mat <- function(com, pair) {
    ncom <- nrow(com)
    distmat <- matrix(nrow=ncom, ncol=ncom, 0, dimnames = list(rownames(com), rownames(com)))
    for (i in 1:(ncom-1)){
      for(j in (i+1):ncom){
        distmat[i, j] = pair[paste(c(rownames(distmat)[i], colnames(distmat)[j]), collapse = "__")]
      }
    }
    # st <- c(0, cumsum(seq(ncom-1, 2))) + 1
    # end <- cumsum(seq(ncom-1, 1))
    # for (i in 1:(ncom-1)) distmat[i, (ncom:(seq(1, ncom)[i]))] = c(pair[end[i]:st[i]],0)
    distmat <- as.dist(t(distmat))
    return(distmat)
    
  } # end of function dist.mat
  
  # pariwise comparisons
  combin  <-  combn(nrow(comm), 2) # table with all pairs
  labcomb <-  apply(combin, 2, function(x) paste(rownames(comm)[x], collapse = "__"))
  colnames(combin) = labcomb
  
  pds <-  pd2(comm, tree, include.root = TRUE) # PD for each community of the community matrix
  pd_sites = pds$pd.root
  names(pd_sites) = pds$site
  com.tot.pair <- 1 * t(apply(combin, 2, function(x) (colSums(comm[x,]) > 0))) # which species in these pairs of comm
  pd.tot.pair0 <- pd2(com.tot.pair, tree, include.root = TRUE)  # PD of the two communities combined
  pd.tot.pair = pd.tot.pair0$pd.root
  names(pd.tot.pair) = pd.tot.pair0$site
  sum.pd.pair <- apply(combin, 2, function(x) sum(pd_sites[x])) # Sum of PD for each community, separetely
  
  min.not.shared <- apply(pd.tot.pair - t(combn(pd_sites, 2)), 1, min) # minimun (b,c)
  names(min.not.shared) = names(pd.tot.pair)
  max.not.shared <- apply(pd.tot.pair - t(combn(pd_sites, 2)), 1, max) # maximum (b,c)
  names(max.not.shared) = names(pd.tot.pair)
  sum.not.shared <- 2*pd.tot.pair - sum.pd.pair  # b+c
  shared <- pd.tot.pair - sum.not.shared    # a
  
  sumSi=sum(pd_sites)
  shared = dist.mat(comm, shared)
  sum.not.shared = dist.mat(comm, sum.not.shared)
  max.not.shared = dist.mat(comm, max.not.shared)
  min.not.shared = dist.mat(comm, min.not.shared)
  
  if(pairwise == FALSE){
    com.tot.multi <- 1 * t(as.matrix(colSums(comm)>0))
    rownames(com.tot.multi) = "all"
    pd.tot.multi <- as.numeric(pd2(com.tot.multi,tree)[,"PD"])  # PD of all communities combined
    St = pd.tot.multi
  }
  
  
  if(pairwise){
    switch(index.family, 
           sorensen = {
             phylo.beta.sim <- min.not.shared/(min.not.shared + shared)
             phylo.beta.sne <- ((max.not.shared - min.not.shared)/((2 * shared) + sum.not.shared)) * (shared/(min.not.shared + shared))
             phylo.beta.sor <- sum.not.shared/(2 * shared + sum.not.shared)
             out <- list(phylo.beta.sim = phylo.beta.sim, phylo.beta.sne = phylo.beta.sne, phylo.beta.sor = phylo.beta.sor)
           },
           jaccard = {
             phylo.beta.jtu <- (2 * min.not.shared)/((2 * min.not.shared) + shared)
             phylo.beta.jne <- ((max.not.shared - min.not.shared)/(shared + sum.not.shared)) * (shared/((2 * min.not.shared) + shared))
             phylo.beta.jac <- sum.not.shared/(shared + sum.not.shared)
             out <- list(phylo.beta.jtu = phylo.beta.jtu, phylo.beta.jne = phylo.beta.jne, phylo.beta.jac = phylo.beta.jac)
           }
    ) # end of switch
  } else { # multi sites
    switch(index.family, 
           sorensen = {
             phylo.beta.SIM <- sum(min.not.shared)/(sumSi - St + sum(min.not.shared))
             phylo.beta.SNE <- ((sum(max.not.shared) - sum(min.not.shared))/(2 * (sumSi - St) + sum(min.not.shared) + sum(max.not.shared))) * ((sumSi - St)/(sumSi - St + sum(min.not.shared)))
             phylo.beta.SOR <- (sum(min.not.shared) + sum(max.not.shared))/(2 * (sumSi - St) + sum(min.not.shared) + sum(max.not.shared))
             out <- list(phylo.beta.SIM = phylo.beta.SIM, phylo.beta.SNE = phylo.beta.SNE, phylo.beta.SOR = phylo.beta.SOR)
           },
           jaccard = {
             phylo.beta.JTU <- (2 * sum(min.not.shared))/((2 * sum(min.not.shared)) + sumSi - St)
             phylo.beta.JNE <- ((sum(max.not.shared) - sum(min.not.shared))/(sumSi - St + sum(max.not.shared) + sum(min.not.shared))) * ((sumSi - St)/(2 * sum(min.not.shared) + sumSi - St))
             phylo.beta.JAC <- (sum(min.not.shared) + sum(max.not.shared))/(sumSi - St + sum(min.not.shared) + sum(max.not.shared))
             out <- list(phylo.beta.JTU = phylo.beta.JTU, phylo.beta.JNE = phylo.beta.JNE, phylo.beta.JAC = phylo.beta.JAC)
           }
    )
  }
  return(out)
} 

#' calculate pairwise beta phylogenetic diversity
#' 
#' A function to calculate a bunch of pairwise beta phylo diversity: 
#' unifraction, MPD, MNTD, PCD
#' 
#' @param samp_wide wide version of samp_long, row.names are sites, colnames are species
#' @param tree a phylogeny with class of 'phylo'
#' @param samp_long a 3-column data frame, site, freq, sp
#' @param get.unif whether to calculate pairwise UniFrac
#' @param get.mpd whether to calculate pairwise mpd_beta
#' @param get.mntd whether to calculate pairwise mntd_beta
#' @param get.pcd calculate PCD or not? Can be time consuming.
#' @param null.model.phylocom whether to run null models for MPD and MNTD with Phylocom?
#' @param null.type.phylocom if null.model.phylocom is TRUE, which null model to use? See ?phylocomr::ph_comstruct for details
#' @param n.item the number of randomization with Phylocom
#' @param abund.weight should abundance information used when calculating MPD/MNTD with Phylocom? Note: mntd_beta with Phylocom seems not reliable.
#' @param null.model.mpd.phylomeasures whether to run null models for MPD with PhyloMeasures?
#' @param null.model.mntd.phylomeasures whether to run null models for MNTD with PhyloMeasures?
#' @param verbose do you want to see relevant information?
#' @param verbose.mntd.null do you want to see relevant information about null models of MNTD with PhyloMeasures?
#' @return a data frame
#' @export
#' 
get_pd_beta = function(samp_wide, tree, samp_long,
                       get.unif = TRUE, get.mpd = TRUE, get.mntd = TRUE, get.pcd = TRUE, 
                       null.model.phylocom = TRUE, null.type.phylocom = 0, 
                       n.item = 999, abund.weight = FALSE,
                       null.model.mpd.phylomeasures = TRUE,
                       null.model.mntd.phylomeasures = FALSE,
                       verbose = TRUE, verbose.mntd.null = FALSE, ...){
  if(length(class(tree)) > 1 & "phylo" %in% class(tree)) class(tree) = "phylo"
  dist = cophenetic(tree)
  
  row.names(samp_wide) = stringr::str_trim(row.names(samp_wide)) %>% 
    tolower() %>% 
    stringr::str_replace_all(" ", "_") # phylocom will have trouble with space in site names
  
  if(missing(samp_long) & abund.weight){
    samp_long = tibble::rownames_to_column(as.data.frame(samp_wide), "site") %>% 
      tidyr::gather("sp", "freq", -site) %>% filter(freq > 0) %>% 
      arrange(site, sp) %>% select(site, freq, sp)
  }
  
  if(get.unif){
    # unifrac: jaccard phylo dissimilarity, only works with presence/absence data
    # unif = unifrac2(samp_wide, tree, samp_long)
    # unif = as.matrix(unif)
    phy_beta = phylo_betapart(samp_wide, tree, index.family = "jaccard", pairwise = TRUE) # false won't work here
    if(verbose) cat("Done with phylo_betapart, Unifrac", "\n")
    unif = as.matrix(phy_beta$phylo.beta.jac)
    unif_turnover = as.matrix(phy_beta$phylo.beta.jtu)
    unif_nested = as.matrix(phy_beta$phylo.beta.jne)
  }
  
  # convert tibble to matrix
  to_m = function(tb){
    mpd_beta = as.data.frame(tb)
    row.names(mpd_beta) = mpd_beta$name; mpd_beta$name = NULL
    mpd_beta = as.matrix(mpd_beta)
    mpd_beta
  }
  
  # mpd_beta
  if(get.mpd){
    if(abund.weight){# weight with abundance, use Phylocom or picante
      mpd_beta_c = try(phylocomr::ph_comdist(samp_long, tree, rand_test = null.model.phylocom, 
                                             null_model = null.type.phylocom, randomizations = n.item, 
                                             abundance = TRUE))
      phylocom_trouble = "try-error" %in% class(mpd_beta_c)
      if(phylocom_trouble){
        if(verbose) cat("Phylocom has trouble with this phlyogney, switch to picante", "\n")
        mpd_beta = picante::comdist(samp_wide, dist, abundance.weighted = TRUE)
        mpd_beta = as.matrix(mpd_beta)
        message("null models for mpd_beta with picante is too slow, ignored.")
      } else {
        if(verbose) cat("Phylocom has no trouble with this phylogeny ", "\n")
        if(null.model.phylocom){
          mpd_beta = to_m(mpd_beta_c$obs)
          mpd_beta_z = to_m(mpd_beta_c$NRI_or_NTI) * (-1)
        } else {
          mpd_beta = to_m(mpd_beta_c)
        }
      }
    } else { # presence/absence, use PhyloMeasures
      mpd_beta = PhyloMeasures::cd.query(tree, samp_wide, standardize = FALSE)
      rownames(mpd_beta) = row.names(samp_wide)
      colnames(mpd_beta) = row.names(samp_wide)
      if(null.model.mpd.phylomeasures){
        mpd_beta_z = PhyloMeasures::cd.query(tree, samp_wide, standardize = TRUE)
        rownames(mpd_beta_z) = row.names(samp_wide)
        colnames(mpd_beta_z) = row.names(samp_wide)
      }
    }
    if(verbose) cat("Done with mpd_beta", "\n")
  }
  
  # mntd_beta
  if(get.mntd){
    if(abund.weight){# weight with abundance, use Phylocom or picante
      mntd_beta_c = try(phylocomr::ph_comdistnt(samp_long, tree, rand_test = null.model.phylocom, 
                                             null_model = null.type.phylocom, randomizations = n.item, 
                                             abundance = TRUE))
      phylocom_trouble = "try-error" %in% class(mntd_beta_c)
      if(phylocom_trouble){
        if(verbose) cat("Phylocom has trouble with this phlyogney, switch to picante", "\n")
        mntd_beta = picante::comdistnt(samp_wide, dist, abundance.weighted = TRUE)
        mntd_beta = as.matrix(mntd_beta)
        message("null models for mntd_beta with picante is too slow, ignored.")
      } else {
        if(verbose) cat("Phylocom has no trouble with this phylogeny ", "\n")
        if(null.model.phylocom){
          mntd_beta = to_m(mntd_beta_c$obs)
          mntd_beta_z = to_m(mntd_beta_c$NRI_or_NTI) * (-1)
        } else {
          mntd_beta = to_m(mntd_beta_c)
        }
      }
    } else { # presence/absence, use PhyloMeasures
      mntd_beta = PhyloMeasures::cdnt.averaged.query(tree, samp_wide)
      rownames(mntd_beta) = row.names(samp_wide)
      colnames(mntd_beta) = row.names(samp_wide)
      if(null.model.mntd.phylomeasures){
        mntd_beta_z = purrr::map(1:n.item, function(x){
          set.seed(x)
          if(verbose.mntd.null) cat("null model of mntd_beta", x, "\t")
          x2 = PhyloMeasures::cdnt.averaged.query(picante::tipShuffle(tree), samp_wide)
          rownames(x2) = row.names(samp_wide)
          colnames(x2) = row.names(samp_wide)
          x2
        })
        mntd_beta_z = array(unlist(mntd_beta_z), dim = c(nrow(samp_wide), nrow(samp_wide), n.item))
        mntd_zz_mean = apply(mntd_beta_z, MARGIN = c(1, 2), mean, na.rm = T)
        mntd_zz_sd = apply(mntd_beta_z, MARGIN = c(1, 2), sd, na.rm = T)
        diag(mntd_zz_sd) = 1
        mntd_beta_z = (mntd_beta - mntd_zz_mean)/mntd_zz_sd
      }
    }
    if(verbose) cat("Done with mntd_beta", "\n")
  }
  
  # pcd
  if(get.pcd){
    pcd_beta = pcd2(comm = samp_wide, tree, unif_dim = 0, verbose = verbose)
    PCD = as.matrix(pcd_beta$PCD)
    PCDc = as.matrix(pcd_beta$PCDc)
    PCDp = as.matrix(pcd_beta$PCDp)
    if(verbose) cat("Done with PCD", "\n")
  }
  
  # clean outputs
  out = tibble::as_data_frame(t(combn(row.names(samp_wide), 2))) %>% 
    rename(site1 = V1, site2 = V2)
  if(get.unif){
    out = mutate(out, unif = purrr::map2_dbl(.x = site1, .y = site2, ~unif[.x, .y]),
                 unif_turnover = purrr::map2_dbl(.x = site1, .y = site2, ~unif_turnover[.x, .y]),
                 unif_nested = purrr::map2_dbl(.x = site1, .y = site2, ~unif_nested[.x, .y]))
  }
  if(get.mpd){
    out = mutate(out, mpd_beta = purrr::map2_dbl(.x = site1, .y = site2, ~mpd_beta[.x, .y]))
  }
  if(get.mntd){
    out = mutate(out, mntd_beta = purrr::map2_dbl(.x = site1, .y = site2, ~mntd_beta[.x, .y]))
  }
  if(get.pcd){
    out = mutate(out, 
                 pcd_beta = purrr::map2_dbl(.x = site1, .y = site2, ~PCD[.x, .y]),
                 pcdc_beta = purrr::map2_dbl(.x = site1, .y = site2, ~PCDc[.x, .y]),
                 pcdp_beta = purrr::map2_dbl(.x = site1, .y = site2, ~PCDp[.x, .y]))
  }
  if(exists("mpd_beta_z")){
    out = mutate(out, 
                 mpd_beta_z = purrr::map2_dbl(.x = site1, .y = site2, ~mpd_beta_z[.x, .y]))
  }
  if(exists("mntd_beta_z")){
    out = mutate(out, 
                 mntd_beta_z = purrr::map2_dbl(.x = site1, .y = site2, ~mntd_beta_z[.x, .y]))
  }
  
  return(out)
}
