

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
#' @param filename The input trees file.
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
#' @todo Make it possible to skip a number of trees upfront and to only 
#' read in a fixed number of trees
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

#' Return a table of the heights and types of all nodes in a tree
#'
#' Internal nodes are labelled "coalescent" and leaves are lablled "sample"
#'
#' Use depth first search to get the tree intervals
#'
#' This can be used to get the vector of speciation/transmission and sampling times for a
#' phylogenetic tree.
#'
#' The results of this function are equivalent to TreeSim::getx()
#'
#' When sampling and coalescent events are at the same times the sampling events are considered first
#'
#' The row names are the node numbers
#'
#' @param tree An object of class "phylo" from ape
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
  types   <- factor(rep(0,n), levels = c("coalescent", "sample"))
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

  ordering  <- order(treetable$heights, decreasing = decreasing)
  treetable <- treetable[ordering,]

  result <- data.frame(node      = row.names(treetable), 
                       nodetype  = treetable$types, 
                       height    = treetable$heights,
                       length    = c(0, diff(treetable$heights)),
                       nlineages = getLineages(treetable$types))
  return(result)
}

#' Get sampling times of a tree
#'
#' @export
getSamplingTimes <- function(tree, ...) {
  intervals <- getTreeIntervals(tree, ...)

  return( intervals$height[intervals$nodetype == "sample"])
}

#' Get branching times of a tree
#'
#' For an ultrametric tree this should be equivalent to the ape::coalescent.intervals
#'
#' @export
getBranchingTimes <- function(tree, ...) {
  intervals <- getTreeIntervals(tree, ...)

  return( intervals$height[intervals$nodetype == "coalescent"])
}




