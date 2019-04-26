###############################################################################
# Utilities for reading BEAST logfiles and getting HPDs
#
# TODO: Port functions in the ESS class of BEAST2 to calculate convergence statistics
#       -> No, use coda instead and ignore the BEAST/BEAST2 functions
#
# TODO: Replace boa functions completely with coda functions


#' Read in a single BEAST logfile and return as a coda mcmc object
#'
#' @param filename The logfile.
#' @param burnin Discard this proportion of samples at the start of the chain.
#' @param maxamples If > 0 stop after reading in this many lines
#'        (this option is only for testing and should generally not be used).
#' @param as.mcmc If FALSE then return a data.frame, otherwise return an mcmc object
#'
read.log.single <- function(filename, burnin=0.1, maxsamples=-1, as.mcmc=TRUE) {
  if (burnin > 1) {
      stop("Error: Burnin must be a proportion of samples between 0 and 1.")
  }

  logfile <- read.table(filename, sep="\t", header=TRUE, nrows=maxsamples)
  n <- nrow(logfile)
  logfile <- logfile[floor(burnin*n):n,]

  if (as.mcmc == TRUE) {
      start <- logfile$Sample[1]
      thin  <- logfile$Sample[2]-logfile$Sample[1]
      rownames(logfile) <- logfile$Sample
      logfile$Sample    <- NULL

      return(coda::mcmc(logfile, start=start, thin=thin))
  } else {
      return(logfile)
  }
}

#' Read in a single BEAST logfile and return as a coda mcmc object.
#' If filenames contains more than one entry each log file is added as a separate chain and
#' a coda mcmc.list object is returned
#'
#' @param filename The logfile.
#' @param burnin Discard this proportion of samples at the start of the chain.
#' @param maxamples If > 0 stop after reading in this many lines
#'        (this option is only for testing and should generally not be used).
#' @param as.mcmc If FALSE then return a data.frame, otherwise return an mcmc object
#'
#' @export
read.log <- function(filenames, burnin=0.1, maxsamples=-1, as.mcmc=TRUE) {

  if (length(filenames) == 1) {
      return( read.log.single(filenames, burnin, maxsamples, as.mcmc) )
  } else {

      loglist <- list()
      for (filename in filenames) {
          loglist[[filename]] <- read.log.single(filename, burnin, maxsamples, as.mcmc)
      }

      if (as.mcmc == TRUE) {
          return (mcmc.list(loglist))
      } else {
          return (loglist)
      }
  }

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


plotPairwiseCorrelations <- function(data, ...) {

  ## put histograms on the diagonal
  panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))

    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE, breaks=20)
    breaks <- h$breaks
    nB <- length(breaks)

    y <- h$counts
    y <- y/max(y)
    #rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
    rect(breaks[-nB], 0, breaks[-1], y, col=pal.dark(cblue,0.5), border=pal.dark(cblue))
  }


  ## put correlations on the upper panels,
  ## with size proportional to the correlations.
  panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))

    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y)
    if (abs(r) < 0.5) {
      col <- pal.dark(cblue)
    } else
      if (abs(r) < 0.75) {
        col <- pal.dark(corange)
      } else {
        col <- pal.dark(cred)
      }
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)

    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r, col = col, font = 2)
  }


  panel.scatter <- function(x,y, alpha=0.2, ... ) {

    r <- cor(x,y)
    if (abs(r) < 0.5) {
      col <- pal.dark(cblue)
    } else
      if (abs(r) < 0.75) {
        col <- pal.dark(corange)
      } else {
        col <- pal.dark(cred)
      }
    panel.smooth(x, y, col = col, col.smooth = pal.dark(cred), ...)
  }

  data <- as.matrix(data)
  pairs(data, panel = panel.scatter, diag.panel = panel.hist, lower.panel = panel.cor,
        cex = 1.2, pch = 20, cex.labels = 1, font.labels = 2, labels=colnames(data), ...)

}
