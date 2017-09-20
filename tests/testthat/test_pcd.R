context("test phylogenetic community dissimilarity, pcd")

test_that("testing pcd_pred, which calculate expectation of conditional PSV", {
  x1 = pcd_pred(comm_old = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
                comm_new = read.csv("../../data/li_2015_new.csv", row.names = 1, check.names = F),
                tree = ape::read.tree("../../data/phy.tre"), reps = 100)
  x2 = pcd_pred(comm_old = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
                      comm_new = read.csv("../../data/li_2015_new.csv", row.names = 1, check.names = F),
                      tree = ape::read.tree("../../data/phy.tre"), reps = 100, cpp = F)
  expect_type(x1, "list")
  expect_length(x1, 4)
  expect_type(x2, "list")
  expect_length(x2, 4)
  expect_equivalent(x1$nsp_pool, x2$nsp_pool)
  expect_equivalent(x1$psv_pool, x2$psv_pool)
  expect_equivalent(x1$nsr, x2$nsr)
  expect_equal(length(x1$psv_bar), length(x2$psv_bar))
})

test_that("testing pcd2, which calculate pairwise site dissimilarity", {
  x1 = pcd_pred(comm_old = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
                comm_new = read.csv("../../data/li_2015_new.csv", row.names = 1, check.names = F),
                tree = ape::read.tree("../../data/phy.tre"), reps = 100)
  x3 = pcd2(comm = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
        tree = ape::read.tree("../../data/phy.tre"),
        expectation = x1)
  x4 = pcd2(comm = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
            tree = ape::read.tree("../../data/phy.tre"),
            expectation = x1, cpp = F)
  expect_type(x3, "list")
  expect_type(x4, "list")
  expect_equivalent(x3, x4)
})

test_that("testing pcd2, without provide expectation", {
  x5 = pcd2(comm = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
            tree = ape::read.tree("../../data/phy.tre"))
  x6 = pcd2(comm = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
            tree = ape::read.tree("../../data/phy.tre"), cpp = F)
  expect_type(x5, "list")
  expect_type(x6, "list")
})

test_that("testing pcd2, expectation based on one community", {
  x1 = pcd_pred(comm_old = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
                tree = ape::read.tree("../../data/phy.tre"), reps = 100)
  x7 = pcd2(comm = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
            tree = ape::read.tree("../../data/phy.tre"),
            expectation = x1)
  x8 = pcd2(comm = read.csv("../../data/li_2015_old.csv", row.names = 1, check.names = F),
            tree = ape::read.tree("../../data/phy.tre"),
            expectation = x1, cpp = F)
  expect_type(x7, "list")
  expect_type(x8, "list")
  expect_equivalent(x7, x8)
})