###############################################################################
# Utilities for reading BEAST logfiles and getting HPDs
#
# TODO: Port functions in the ESS class of BEAST2 to calculate convergence statistics
#
# TODO: Should coda be used instead of boa?


#' Read in BEAST logfile
#'
#' @param filename The logfile.
#' @param burnin Discard this proportion of samples at the start of the chain.
#' @param maxamples If > 0 stop after reading in this many lines
#'        (this option is only for testing and should generally not be used).
#'
#' @export
readLogfile <- function(filename, burnin=0.1, maxsamples=-1) {
  if (burnin > 1) {
      stop("Error: Burnin must be a proportion of samples between 0 and 1.")
  }
  logfile <- read.table(filename, sep="\t", header=TRUE, nrows=maxsamples)
  n <- nrow(logfile)
  return(logfile[floor(burnin*n):n,])
}



toCoda <- function(lf) {

  start <- lf$Sample[1]
  thin  <- lf$Sample[2]-lf$sample[1]
  rownames(lf) <- lf$Sample

  return(coda::mcmc(lf, start=start, thin=thin))
}

#' Extract all matching parameters from the logfile
#'
#' if par="R0" extract (R0s.1 R0s.2 R0s.3 etc.)
#'
#' @param par The string to match with the start of the parameters to extract.
#' @return A number of parameters from the logfile that match "^<par>"
#'
#' @export
getLogFileSubset <- function(logfile, par) {
  return(logfile[grepl(paste0('^',par), names(logfile))])
}


#' Calculate HPD of a parameter
#'
#' Get HPD of a posterior sample
#'
#' Uses Chen and Shao algorithm as implemented in boa package.
#'
#' @param data The samples from the posterior.
#' @param alpha The confidence level.
#' @param includleMedian If FALSE only return the HPD limits, and not the median as well.
#' @return c(lower, median, upper)
#'
#' @export
getHPD <- function(data, alpha=0.05, includeMedian=TRUE) {
  hpd <- boa::boa.hpd(data, alpha)
  if (includeMedian) {
    med <- median(data)
    out <- c("lower"=hpd[[1]], "median"=med, "upper"=hpd[[2]])
  } else {
    names(hpd) <- c("lower","upper")
    out <- hpd
  }
  return(out)
}


#' Calculate HPD of a set of parameters
#'
#' Get HPD of a matrix of values e.g. a skyline
#'
#' Uses Chen and Shao algorithm as implemented in boa package.
#'
#' @param data The samples from the posterior. Can be in a data frame or in a matrix.
#'             By default the function assumes that each row represents a sample from
#'             the posterior and each column a parameter.
#' @param margin If the data are formatted with samples in columns set margin to 1.
#' @param ...  Parameters to be passed to getHPD().
#' @param dataframe If TRUE return output as a data frame, otherwise a matrix.
#' @return 3xn matrix or data frame where each column is c(lower, median, upper) for
#'         each parameter in turn.
#'
#' @export
getMatrixHPD <- function(data, margin=2, dataframe=TRUE, ...) {
  out <- apply(data, margin, getHPD, ...)
  if (dataframe) {
    return(data.frame(out))
  } else {
    return(out)
  }
}


#' Calculate pairwise correlations between parameters
#'
#' Returns a matrix with the pairwise correlations between pairs of parameters.
#' This function is similar to \code{coda::corsscorr}, but unlike the
#' \code{coda} function other parameters can be passed to the correlation
#' function (so spearman or pearson correlation coefficients can be
#' calculated).
#'
#' @param data Data frame or matrix of the posterior samples of parameters.
#'   If data is a matrix, samples are rows and parameters are columns.
#'
pairwiseCorrelations <- function(data, ...) {
  return (cor(as.matrix(data, ...)))
}
