% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_utils.R
\name{getTreeIntervals}
\alias{getTreeIntervals}
\title{Intervals of a phylogenetic tree}
\usage{
getTreeIntervals(tree, decreasing = FALSE, raw = FALSE)
}
\arguments{
\item{tree}{An object of class "phylo" from ape}

\item{decreasing}{If FALSE intervals are ordered from the present into
the past. If TRUE intervals are ordered from the tMRCA to the
present.}

\item{raw}{If TRUE return the raw table before reordering and calculating
interval lengths or numbers of lineages. Mostly for debugging
purposes.}
}
\value{
A table with each row representing one of the tree intervals.
}
\description{
Function to get the coalescent and sampling intervals of a binary tree
(the function also works for serially-sampled trees).
The function returns a table with each row representing one of the
intervals in the tree. This can be used to get the vector of
speciation/transmission and sampling times for a phylogenetic tree.
}
\details{
Each row of the output table represents an interval of the tree. Intervals
are between successive nodes in the tree, with time measured from the present
(\eqn{t = 0}) into the past. For an interval between time \eqn{x} and \eqn{y},
with \eqn{x < y}, the function returns

\itemize{
    \item The node number of the node at time \eqn{y} in \code{tree}.
    \item The node type of the node at time \eqn{y} ("coalescent" for
          internal nodes and "sample" for leaves).
    \item The height of the node at the end of the interval (\eqn{y}).
    \item The length of the interval (\eqn{y - x}).
    \item The number of lineages in the tree during the interval.
}

Rows are ordered based on the end time of the intervals. If two nodes have the
same height, nodes are ordered by type ("coalescent" < "sample") and then by
node number.

To get the tree interval the function performs a depth first traversal of
the tree.

The results of this function are equivalent to \code{\link[TreeSim]{getx}}
in the \code{TreeSim} package, but more verbose.
}
\examples{

}
\seealso{
\code{\link{getBranchingTimes}}, \code{\link{getSamplingTimes}},
         \code{\link[TreeSim]{getx}}, \code{\link[ape]{branching.times}},
         \code{\link[ape]{coalescent.intervals}}
}
