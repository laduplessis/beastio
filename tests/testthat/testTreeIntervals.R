context("Tree Intervals")

test_that("getTreeIntervals works on a homochronous tree", {
  tree   <- ape::read.tree("../testTree1.tree")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,               c(0,0,0,0,0,0,11,23,25,30,40))
  expect_equal(as.numeric(result$nodetype), c(2,2,2,2,2,2, 1, 1, 1, 1, 1))
  expect_equal(result$nlineages,            c(0,1,2,3,4,5, 6, 5, 4, 3, 2))
  expect_equal(result$length,               c(0,0,0,0,0,0,11,12, 2, 5,10))
})

test_that("getTreeIntervals works on a heterochronous trees", {
  tree   <- ape::read.tree("../testTree2.tree")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,               c(0, 8,10,15,15,16,17,17,22))
  expect_equal(as.numeric(result$nodetype), c(2, 2, 2, 2, 1, 1, 2, 1, 1))
  expect_equal(result$nlineages,            c(0, 1, 2, 3, 4, 3, 2, 3, 2))
  expect_equal(result$length,               c(0, 8, 2, 5, 0, 1, 1, 0, 5))

  tree   <- ape::read.tree("../testTree3.tree")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,               c(0,0,6,6,11,20,23,25,28,30,40))
  expect_equal(as.numeric(result$nodetype), c(2,2,2,2, 1, 2, 1, 1, 2, 1, 1))
  expect_equal(result$nlineages,            c(0,1,2,3, 4, 3, 4, 3, 2, 3, 2))
  expect_equal(result$length,               c(0,0,6,0, 5, 9, 3, 2, 3, 2,10))

  tree   <- ape::read.tree("../testTree4.tree")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,               c(0,0,6,6,11,20,23,23,28,30,40))
  expect_equal(as.numeric(result$nodetype), c(2,2,2,2, 1, 2, 1, 1, 2, 1, 1))

  tree   <- ape::read.tree("../testTree5.tree")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,               c(0,0,6,6,11,20,23,28,30,32,40))
  expect_equal(as.numeric(result$nodetype), c(2,2,2,2, 1, 2, 1, 2, 1, 1, 1))
})

test_that("getTreeIntervals works on a tree with 2 nodes", {
  tree   <- ape::read.tree(text = "(1:2.0,2:3.0);")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,           c(0,1,3))
  expect_equal(as.numeric(result$nodetype), c(2,2,1))

  tree   <- ape::read.tree(text = "(1:2.0,2:2.0);")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,           c(0,0,2))
  expect_equal(as.numeric(result$nodetype), c(2,2,1))
})

test_that("getTreeIntervals works on a tree with 3 nodes", {
  tree   <- ape::read.tree(text = "(1:6.0,(2:3.0,3:4.0):7.0);")
  result <- getTreeIntervals(tree)
  expect_equal(result$height,           c(0,1,4,5,11))
  expect_equal(as.numeric(result$nodetype), c(2,2,1,2, 1))
})

test_that("getSamplingTimes works", {
  tree   <- ape::read.tree("../testTree3.tree")
  expect_equal(getSamplingTimes(tree),   c(0,0,6,6,20,28))

  tree   <- ape::read.tree("../testTree1.tree")
  expect_equal(getSamplingTimes(tree),   c(0,0,0,0,0,0))
})

test_that("getBranchingTimes works", {
  tree     <- ape::read.tree("../testTree3.tree")
  expected <- c(11,23,25,30,40)

  expect_equal(getBranchingTimes(tree, decreasing = FALSE), expected)
  expect_equal(getBranchingTimes(tree, decreasing = TRUE),  rev(expected))
})
