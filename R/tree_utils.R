###############################################################################
# Utilities for reading BEAST tree files and getting tree intervals
#
# TODO: Make it possible to only read in some number of trees or to buffer
# reading in big trees files
#

#' Read BEAST trees file
#'
#' This function reads in a \code{.trees} file written by BEAST,
#' BEAST2 or one of the associated packages. The input file
#' can contain a single tree or a set of posterior trees. The input
#' file must be in NEXUS format and must be terminated by "END;"
#'
#' In practice the function is a wrapper for \code{\link[ape]{read.nexus}},
#' with the ability to easily discard some proportion of trees as burnin-in.
#'
#' @section Warning:
#'   The function will attempt to read the whole trees file. If the
#'   file is large this may take a long time or run out of memory.
#'
#' @param filename The name of the trees file to read.
#' @param burnin Discard this proportion of trees at the start of the file.
#'
#' @return An object of class "phylo" or "multiphylo" containing one or many
#' trees.
#'
#' @seealso \code{\link[ape]{read.nexus}}, \code{\link[ape]{read.tree}}
#'
#' @examples
#' # Read trees file with 10% burnin
#' readTreeLog(filename)
#'
#' # Read trees file without a burnin
#' readTreeLog(filename, burnin=0)
#'
#'
#' @export
readTreeLog <- function(filename, burnin=0.1) {
  if (burnin > 1) {
      stop("Error: Burnin must be a proportion of samples between 0 and 1.")
  }
  treesfile <- ape::read.nexus(filename)
  n <- length(treesfile)
  return(treesfile[floor(burnin*n):n])
}


#' Internal function for depth first traversal of a tree in order to get
#' the node heights
treeDFS <- function(i, t, treetable) {

  t <- t + treetable$edgelen[i]
  treetable$heights[i] <- t

  children <- which(treetable$parents == i)
  if (length(children) == 0) {
    treetable$types[i] <- "sample"
  } else {
    treetable$types[i] <- "coalescent"
    for (c in children) {
      treetable <- treeDFS(c, t, treetable)
    }
  }

  return(treetable)
}

#' Internal function to get the number of lineages during an interval
getLineages <- function(types) {
  n <- length(types)
  increment <- rep(-1, n)
  increment[types == "sample"] <- 1

  return( c(0, cumsum(increment[1:(n-1)])) )
}

#' Intervals of a phylogenetic tree
#' 
#' Function to get the coalescent and sampling intervals of a binary tree
#' (the function also works for serially-sampled trees).
#' The function returns a table with each row representing one of the 
#' intervals in the tree. This can be used to get the vector of 
#' speciation/transmission and sampling times for a phylogenetic tree.
#' 
#' Each row of the output table represents an interval of the tree. Intervals
#' are between successive nodes in the tree, with time measured from the present
#' (\eqn{t = 0}) into the past. For an interval between time \eqn{x} and \eqn{y}, 
#' with \eqn{x < y}, the function returns
#' 
#' \itemize{
#'     \item The node number of the node at time \eqn{y} in \code{tree}.
#'     \item The node type of the node at time \eqn{y} ("coalescent" for 
#'           internal nodes and "sample" for leaves).
#'     \item The height of the node at the end of the interval (\eqn{y}).
#'     \item The length of the interval (\eqn{y - x}).
#'     \item The number of lineages in the tree during the interval.
#' }
#' 
#' Rows are ordered based on the end time of the intervals. If two nodes have the 
#' same height, nodes are ordered by type ("coalescent" < "sample") and then by 
#' node number.
#' 
#' To get the tree interval the function performs a depth first traversal of 
#' the tree.
#'
#' The results of this function are equivalent to \code{\link[TreeSim]{getx}} 
#' in the \code{TreeSim} package, but more verbose. 
#'
#'
#' @param tree An object of class "phylo" from ape
#' @param decreasing If FALSE intervals are ordered from the present into 
#'        the past. If TRUE intervals are ordered from the tMRCA to the 
#'        present.
#'        
#' @return A table with each row representing one of the tree intervals.
#'
#' @seealso \code{\link[TreeSim]{getx}}, \code{\link[ape]{branching.times}},
#'          \code{\link[ape]{coalescent.intervals}}
#'          
#' @examples
#'
#' @export
getTreeIntervals <- function(tree, decreasing=FALSE) {

  if (class(tree) != "phylo") {
      stop("Input tree is not of class \"phylo\".")
  }

  # Total number of nodes in the tree
  n <- tree$Nnode + length(tree$tip.label)

  # Initialise arrays
  parents <- numeric(n)
  edgelen <- numeric(n)
  heights <- numeric(n)
  #types   <- factor(rep(0,n), levels = c("sample", "coalescent"))
  types   <- ordered(rep(0,n), levels = c("coalescent", "sample"))
  treetable <- data.frame(parents, edgelen, heights, types)

  # Get parents and edge lengths
  for (i in 1:nrow(tree$edge)) {
    p <- tree$edge[i,1]
    c <- tree$edge[i,2]
    treetable$parents[c] <- p
    treetable$edgelen[c] <- tree$edge.length[i]
  }

  # Fill in heights and types
  rootnode  <- which(treetable$parents == 0)
  treetable <- treeDFS(rootnode, 0, treetable)

  # Flip and order
  tmrca <- max(treetable$heights)
  treetable$heights <- tmrca - treetable$heights

  ordering  <- order(treetable$heights, treetable$types, row.names(treetable), decreasing = decreasing)
  treetable <- treetable[ordering,]

  result <- data.frame(node      = row.names(treetable),
                       nodetype  = treetable$types,
                       height    = treetable$heights,
                       length    = c(0, diff(treetable$heights)),
                       nlineages = getLineages(treetable$types))
  return(result)
}

#' Sampling times of a phylogenetic tree
#' 
#' Return an ordered list of sampling times in a phylogenetic tree.
#' For ultrametric trees this function should return a list of 0's.
#' 
#' @inheritParams getTreeIntervals
#' @param ... Extra parameters to pass to \code{\link{getTreeIntervals}}.
#'
#' @return A list of times representing the heights of the leaf nodes in 
#'         the tree.
#'         
#' @seealso \code{\link{getTreeIntervals}}, \code{\link[TreeSim]{getx}}
#'
#' @export
getSamplingTimes <- function(tree, ...) {
  intervals <- getTreeIntervals(tree, ...)

  return( intervals$height[intervals$nodetype == "sample"])
}



#' Branching (coalescent) times of a phylogenetic tree
#' 
#' Return an ordered list of branching (coalescent) times in a phylogenetic tree.
#' For an ultrametric tree this should be equivalent to 
#' \code{\link[ape]{coalescent.intervals}}.
#' 
#' @inheritParams getTreeIntervals
#' @param ... Extra parameters to pass to \code{\link{getTreeIntervals}}.
#'
#' @return A list of times representing the heights of the leaf nodes in 
#'         the tree.
#'         
#' @seealso \code{\link{getTreeIntervals}}, \code{\link[TreeSim]{getx}}, 
#'          \code{\link[ape]{branching.times}}, \code{\link[ape]{coalescent.intervals}}
#'
#'
#' @export
getBranchingTimes <- function(tree, ...) {
  intervals <- getTreeIntervals(tree, ...)

  return( intervals$height[intervals$nodetype == "coalescent"])
}




