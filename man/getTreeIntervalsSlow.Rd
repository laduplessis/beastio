% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_utils.R
\name{getTreeIntervalsSlow}
\alias{getTreeIntervalsSlow}
\title{Intervals of a phylogenetic tree}
\usage{
getTreeIntervalsSlow(tree, decreasing = FALSE, raw = FALSE)
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
\description{
Function to get the coalescent and sampling intervals of a binary tree
(the function also works for serially-sampled trees).
The function returns a table with each row representing one of the
intervals in the tree. This can be used to get the vector of
speciation/transmission and sampling times for a phylogenetic tree.
}
\details{
This is a deprecated version that is very slow on big trees (>10K nodes).
This version is kept because it is easier to read and debug.
Not using random access to arrays speeds up the initial setup.
}
\seealso{
\code{\link{getTreeIntervals}}
}
