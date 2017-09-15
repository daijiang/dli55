Untitled
================

This package holds functions I wrote or collected that are useful for myself.

PCD
---

Adopted from the `picante` R package, wrote by Ives and Helmus 2010.

``` r
library(dli55)
x1.c = pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
              comm_new = read.csv("data/li_2015_new.csv", row.names = 1, check.names = F),
                  tree = ape::read.tree("data/phy.tre"), reps = 100)
x1.r = pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
                comm_new = read.csv("data/li_2015_new.csv", row.names = 1, check.names = F),
                tree = ape::read.tree("data/phy.tre"), reps = 100, cpp = F)
x1.c; x1.r
```

    ## $nsp_pool
    ## [1] 26
    ## 
    ## $psv_bar
    ##         1         6         7         5         2         4         3 
    ## 0.8122259 0.5436851 0.4795122 0.5791422 0.7079984 0.5898115 0.6668727 
    ##         9 
    ## 0.4523102 
    ## 
    ## $psv_pool
    ## [1] 0.8414645
    ## 
    ## $nsr
    ## [1] 1 6 7 5 2 4 3 9

    ## $nsp_pool
    ## [1] 26
    ## 
    ## $psv_bar
    ##         1         6         7         5         2         4         3 
    ## 0.7796509 0.5201671 0.4963931 0.5710461 0.6443410 0.6114685 0.7043220 
    ##         9 
    ## 0.4066287 
    ## 
    ## $psv_pool
    ## [1] 0.8414645
    ## 
    ## $nsr
    ## [1] 1 6 7 5 2 4 3 9

``` r
library(microbenchmark)
microbenchmark(pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
                        comm_new = read.csv("data/li_2015_new.csv", row.names = 1, check.names = F),
                        tree = ape::read.tree("data/phy.tre"), reps = 100, cpp = FALSE),
               pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
                        comm_new = read.csv("data/li_2015_new.csv", row.names = 1, check.names = F),
                        tree = ape::read.tree("data/phy.tre"), reps = 100, cpp = TRUE),
               times = 60)
```

    ## Unit: milliseconds
    ##                                                                                                                                                                                                                                               expr
    ##  pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1,      check.names = F), comm_new = read.csv("data/li_2015_new.csv",      row.names = 1, check.names = F), tree = ape::read.tree("data/phy.tre"),      reps = 100, cpp = FALSE)
    ##   pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1,      check.names = F), comm_new = read.csv("data/li_2015_new.csv",      row.names = 1, check.names = F), tree = ape::read.tree("data/phy.tre"),      reps = 100, cpp = TRUE)
    ##        min        lq      mean    median        uq       max neval cld
    ##  38.330547 39.922463 41.334938 41.088976 42.420504 47.929658    60   b
    ##   5.508976  5.704419  5.972851  5.815033  6.037773  8.334505    60  a

``` r
microbenchmark(pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
                    tree = ape::read.tree("data/phy.tre"), 
                    nsp_pool = x1.c$nsp_pool, 
                    PSV_bar = x1.c$psv_bar, 
                    PSV_pool = x1.c$psv_pool, 
                    nsr = x1.c$nsr, cpp = F),
               pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
                    tree = ape::read.tree("data/phy.tre"), 
                    nsp_pool = x1.c$nsp_pool, 
                    PSV_bar = x1.c$psv_bar, 
                    PSV_pool = x1.c$psv_pool, 
                    nsr = x1.c$nsr, cpp = T),
               times = 20)
```

    ## Unit: milliseconds
    ##                                                                                                                                                                                                                                      expr
    ##  pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),      tree = ape::read.tree("data/phy.tre"), nsp_pool = x1.c$nsp_pool,      PSV_bar = x1.c$psv_bar, PSV_pool = x1.c$psv_pool, nsr = x1.c$nsr,      cpp = F)
    ##  pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),      tree = ape::read.tree("data/phy.tre"), nsp_pool = x1.c$nsp_pool,      PSV_bar = x1.c$psv_bar, PSV_pool = x1.c$psv_pool, nsr = x1.c$nsr,      cpp = T)
    ##       min       lq     mean   median       uq      max neval cld
    ##  217.0081 219.3897 233.8237 223.5000 227.9881 414.5605    20   b
    ##  156.4941 161.4989 168.2967 164.5911 170.1828 236.1666    20  a
