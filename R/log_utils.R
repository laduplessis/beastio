###############################################################################
# Utilities for reading BEAST logfiles and getting HPDs
#
# TODO: Port functions in the ESS class of BEAST2 to calculate convergence statistics
#       -> No, use coda instead and ignore the BEAST/BEAST2 functions
#


#' Read in a single BEAST logfile and return as a coda mcmc object
#'
#' @param filename The logfile.
#' @param burnin Discard this proportion of samples at the start of the chain.
#' @param maxamples If > 0 stop after reading in this many lines
#'        (this option is only for testing and should generally not be used).
#' @param as.mcmc If FALSE then return a data.frame, otherwise return an mcmc object
#'
readSingleLog <- function(filename, burnin=0.1, maxsamples=-1, as.mcmc=TRUE) {
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
readLog <- function(filenames, burnin=0.1, maxsamples=-1, as.mcmc=TRUE) {

  if (length(filenames) == 1) {
      return( readSingleLog(filenames, burnin, maxsamples, as.mcmc) )
  } else {

      loglist <- list()
      for (filename in filenames) {
          loglist[[filename]] <- readSingleLog(filename, burnin, maxsamples, as.mcmc)
      }

      if (as.mcmc == TRUE) {
          return (coda::mcmc.list(loglist))
      } else {
          return (loglist)
      }
  }

}



#' Extract all matching parameters from the logfile
#' The object returned is always in the same format as logfile
#' (mcmc, mcmc.list, data.frame, or list<data.frame>)
#'
#' e.g if par="R0" extract (R0s.1 R0s.2 R0s.3 etc.)
#'
#' @param logfile The logfile object to subset
#' @param pattern The string to match with the start of the parameters to extract. Can be a regular expression
#' @param start   If TRUE then only the start of parameters are matched against pattern (pattern becomes "^<pattern>")
#' @return A logfile object with only those parameters that match pattern
#'
#' @export
getLogFileSubset <- function(logfile, pattern, start=FALSE) {
  if (start == TRUE) {
      pattern <- paste0("^",pattern)
  }

  if (coda::is.mcmc(logfile) || coda::is.mcmc.list(logfile)) {
      return(logfile[, grepl(pattern, varnames(logfile)), drop=TRUE])
  } else {

      if (is.data.frame(logfile)) {

        # Single data frame
        return(logfile[, grepl(pattern, names(logfile)), drop=TRUE])
      } else {

        # List of data frames
        result <- list()
        for (chain in names(logfile)) {
            result[[chain]] <- logfile[[chain]][, grepl(pattern, names(logfile[[chain]])), drop=TRUE]
        }
        return(result)
      }
  }
  #return(logfile[, grepl(paste0('^',par), names(logfile))])
}


#' Apply a function to an mcmc or mcmc.list object, e.g. median
#'
#' @param data The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#' @param fun  The function to apply to the MCMC sample, by default median
#' @param ...  Extra parameters to be passed to the function
#'
#' @export
applyMCMC <- function(data, fun=median, ...) {

  applySingleMCMC <- function(data) { apply(data, 2, fun, ...) }

  if (coda::is.mcmc(data)) {
    return(applySingleMCMC(data))
  } else
  if (coda::is.mcmc.list(data)) {

    result <- list()
    for (chain in chanames(data)) {
      result[[chain]] <- applySingleMCMC(data[[chain]])
    }
    return(result)

  } else {
    stop("Error: Not an mcmc or mcmc.list object")
  }

}

#' Calculate the HPD interval and median and return in a single object
#' Only works for mcmc and mcmc.list objects
#'
#' @param data The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#' @param prob Target probability content for the HPD interval (see ?HPDinterval)
#' @param ...  Extra parameters to be passed to the median function
#'
#' @export
getHPDMedian <- function(data, prob=0.95, ...) {

  t_hpd <- coda::HPDinterval(data, prob=prob)
  t_med <- applyMCMC(data, median, ...)

  pasteAndReorder <- function(hpd, med) { cbind(hpd, med)[,c(1,3,2)] }


  if (coda::is.mcmc(data)) {
    return( pasteAndReorder(t_hpd, t_med) )
  } else
  if (coda::is.mcmc.list(data)) {

    result <- list()
    for (chain in chanames(data)) {
      result[[chain]] <- pasteAndReorder(t_hpd[[chain]], t_med[[chain]])
    }
    return(result)

  } else {
    stop("Error: Not an mcmc or mcmc.list object")
  }

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
getHPD.boa <- function(data, alpha=0.05, includeMedian=TRUE) {
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
getMatrixHPD.boa <- function(data, margin=2, dataframe=TRUE, ...) {
  out <- apply(data, margin, getHPD, ...)
  if (dataframe) {
    return(data.frame(out))
  } else {
    return(out)
  }
}

#' Like constantPars, but for a single mcmc object
#' (private helper function)
constantParsMCMC <- function(data) {

  if (coda::is.mcmc(data)) {

      result <- c()
      for (par in varnames(data)) {
        if (all(data[ ,par, drop=FALSE] == data[1, par])) {
          result <- c(result, par)
        }
      }
      return(result)

  } else {
      return(NULL)
  }
}


#' Identify parameters in the mcmc chain that have not been sampled
#' i.e. parameters that have stayed constant across all chains
#'
#' If some parameters are constant on some chains but not on others
#' an error is thrown
#'
#' @param data The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#'
#' @export
constantPars <- function(data) {

  if (coda::is.mcmc(data)) {
      data <- coda::mcmc.list(data)
  }

  if (coda::is.mcmc.list(data)) {

      result <- constantParsMCMC(data[[1]])
      if (nchain(data) > 1) {
          for (i in 2:nchain(data)) {
              t <- constantParsMCMC(data[[i]])

              if (!all.equal(result,t)) {
                  stop(paste0("Error: Chain ",i," has a different set of constant parameters to chain 1. Possibly some parameters got stuck because of a poor proposal distribution."))
              }
          }
      }
      return(result)
  } else {
      stop("Error: Not an mcmc or mcmc.list object")
  }
}

#' The opposite to constantPars, that is returns all parameters
#' that have been sampled
#'
#' @param data The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#'
#' @export
sampledPars <- function(data) {
    return(setdiff(varnames(data), constantPars(data)))
}



#' Check if any parameters in an mcmc chain have an ESS value < cutoff
#' Constant (unsampled) parameters are ignored.
#'
#' @param data   The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#' @param cutoff The cutoff ESS value (default = 200)
#' @param value  If value == TRUE the parameters with ESS < cutoff and their ESS values are returned
#'               If value == FALSE only the column indices of those parameters are returned
#' @param ignored Parameters in the chain that should be ignored
#'
#' @export
checkESS <- function(data, cutoff = 200, value = TRUE, ignored = c()) {

  data <- as.mcmc(data)
  if (coda::is.mcmc(data) && !is.list(data)) {
      data <- data[, setdiff(sampledPars(data), ignored)]
      ess  <- coda::effectiveSize(data)

      if (value == TRUE) {
          return(ess[ess < cutoff])
      } else {
          return(which(ess < cutoff))
      }
  } else {
    stop("Error: Not an mcmc object")
  }
  #result <- lapply(lapply(mcmc.trace[, p], effectiveSize), function(x) x[x < 500])
}


#' Check that all parameters in an mcmc chain have an ESS value >= cutoff
#'
#' @export
checkESSl <- function(data, cutoff = 200, ignored = c()) {
  return (length(checkESS(data, cutoff, value = FALSE, ignored)) == 0)
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
