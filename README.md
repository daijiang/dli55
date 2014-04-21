## Functions collected or wrote by Daijiang Li

This package contains R functions collected or wrote by Daijiang Li. By putting all functions into one package, I can keep and use them easily.

## Installation
To install this package, just run:

    library(devtools)
    install_github("dli55", "daijiang")
    
## Functions

- multiplot(): put multiple ggplot figures into one page.
- wide_by_site(): change long table of vegetatino data into wide table for each site (quadrat by site matrix).
- wide_by_site_rarefy(): simialr as above, but randomly sample q quadrats out of sites' total quadrats n times. Returns a list of lists.
- cscore_site(): calculate C-score for each site using null models.
    

