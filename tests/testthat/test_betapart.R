context("test partitioning of phylogenetic beta diversity")

x = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F)
x[x > 1] = 1
tree = ape::read.tree("../../data/phy.tre")
tol = 0.0001

test_that("testing phylo_betapart and betapart::phylo.beta.pair", {
  expect_equal(phylo_betapart(x, tree, "jaccard", pairwise = T), betapart::phylo.beta.pair(x, tree, "jaccard"), tolerance = tol)
  expect_equal(phylo_betapart(x, tree, "sorensen", pairwise = T), betapart::phylo.beta.pair(x, tree, "sorensen"), tolerance = tol)
  expect_equal(phylo_betapart(x, tree, pairwise = T)$phylo.beta.jac, unifrac2(x, tree), tolerance = tol, check.attributes = F)
  expect_equal(phylo_betapart(x, tree, pairwise = F), betapart::phylo.beta.multi(x, tree, "jaccard"), tolerance = tol)
  expect_equal(phylo_betapart(x, tree, pairwise = F, "sorensen"), betapart::phylo.beta.multi(x, tree), tolerance = tol)
})

# test_that("testing whether get_pd_beta works or not", {
#   expect_success(get_pd_alpha(x, tree, null.model = T, n.item = 10))
#   expect_success(get_pd_beta(x, tree, null.model = T, n.item = 10))
# })
