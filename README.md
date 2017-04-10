## Functions collected or wrote by Daijiang Li

This package contains R functions collected or wrote by Daijiang Li. By putting all functions into one package, I can keep and use them easily.

## Installation
To install this package, just run:

    library(devtools)
    install_github("dli55", "daijiang")
    
## Functions

- `cscore_site()`: calculate C-score for each site using null models.
- `check_names()`: check plant species names for typos and synonyms.
- `dist_to_df()`: convert a distance matrix to data frame.
- `logit_tran()`: logit transforme a vector.
- `multiplot()`: put multiple ggplot figures into one page.
- `panel.cor()`: used for pairs().
- `panel.hist()`: used for pairs().
- `rand_test()`: randomization tests.
- `rm_sp_noobs()`: remove species not observed from a site by species matrix.
- `simpleCap()`: convert the first letter to upper case.
- `var_to_rownames()`: move a column to be row name.
- `wide_by_site()`: change long table of vegetatino data into wide table for each site (quadrat by site matrix).
- `wide_by_site_rarefy()`: simialr as above, but randomly sample q quadrats out of sites' total quadrats n times. Returns a list of lists.

