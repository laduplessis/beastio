

#' Read in BEAST trees file
#'
#' The file can contain a single tree or a set of posterior trees. The
#' file must be in NEXUS format and must be terminated by END;
#'
#' @section Warning:
#'   The function will attempt to read the whole trees file. If the
#'   file is large this may take a long time or run out of memory.
#'
#' @param filename The trees file.
#' @param burnin Discard this proportion of samples at the start of the chain.
#'
#' @export
readTreefile <- function(filename, burnin=0.1) {
  if (burnin > 1) {
    stop("Error: Burnin must be a proportion of samples between 0 and 1.")
  }
  treesfile <- ape::read.nexus(filename)
  n <- length(treesfile)
  return(treesfile[floor(burnin*n):n])
}

