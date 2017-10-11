Untitled
================

This package holds functions I wrote or collected that are useful for myself.

You can install it with `devtools::install_github("daijiang/dli55")`, but you will likely get compile errors. I have also iploaded a binary version, which you can install with `install.packages("https://raw.githubusercontent.com/daijiang/dli55/master/dli55_1.0.1.tgz", repos = NULL)`.

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
    ## 0.7656272 0.5062540 0.4515079 0.5928342 0.7105119 0.6037607 0.6264553 
    ##         9 
    ## 0.4061341 
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
    ## 0.7842456 0.5389647 0.4775832 0.5607200 0.7266255 0.6313489 0.6554613 
    ##         9 
    ## 0.4063451 
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
               times = 20)
```

    ## Unit: milliseconds
    ##                                                                                                                                                                                                                                               expr
    ##  pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1,      check.names = F), comm_new = read.csv("data/li_2015_new.csv",      row.names = 1, check.names = F), tree = ape::read.tree("data/phy.tre"),      reps = 100, cpp = FALSE)
    ##   pcd_pred(comm_old = read.csv("data/li_2015_old.csv", row.names = 1,      check.names = F), comm_new = read.csv("data/li_2015_new.csv",      row.names = 1, check.names = F), tree = ape::read.tree("data/phy.tre"),      reps = 100, cpp = TRUE)
    ##        min        lq      mean   median        uq       max neval cld
    ##  44.583163 48.606821 58.277472 49.85511 70.290546 104.35699    20   b
    ##   7.182271  7.620585  9.168821  8.90935  9.997768  14.14568    20  a

``` r
microbenchmark(pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
                    tree = ape::read.tree("data/phy.tre"), 
                    expectation = x1.c, cpp = F),
               pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),
                    tree = ape::read.tree("data/phy.tre"), 
                    expectation = x1.c, cpp = T),
               times = 20)
```

    ## Unit: milliseconds
    ##                                                                                                                                                         expr
    ##  pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),      tree = ape::read.tree("data/phy.tre"), expectation = x1.c,      cpp = F)
    ##  pcd2(comm = read.csv("data/li_2015_old.csv", row.names = 1, check.names = F),      tree = ape::read.tree("data/phy.tre"), expectation = x1.c,      cpp = T)
    ##       min       lq     mean   median       uq      max neval cld
    ##  276.6306 279.6848 291.9925 282.3650 291.1818 387.1725    20   b
    ##  199.4720 209.0499 232.6000 213.1885 238.1650 332.6463    20  a
