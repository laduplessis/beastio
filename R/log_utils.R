###############################################################################
# Utilities for reading BEAST log files and getting HPDs
#
# TODO: Port functions in the ESS class of BEAST2 to calculate convergence statistics
#       -> No, use coda instead and ignore the BEAST/BEAST2 functions
#

# require(coda)
# require(gplots)
# require(boa)


#' Read a single BEAST log file
#'
#' Read a single BEAST log file and return as a coda "mcmc" object.
#'
#' @param filename The name of the log file to read.
#' @param burnin Discard this proportion of samples at the start of the chain
#'        (if `burninAsSamples` == FALSE). Otherwise discard this many samples
#'        at the start of the chain.
#' @param maxsamples If > 0 stop after reading in this many lines
#'        (this option is only for testing and should generally not be used).
#' @param as.mcmc If FALSE then return an object of class "data.frame", else
#'        return an "mcmc" object
#' @param burninAsSamples if TRUE burnin is given as the number of samples,
#'        if FALSE burnin is a proportion (0 <= burnin < 1) of samples.
#'        (default = FALSE).
#'
#' @return An "mcmc" object (\code{\link[coda]{mcmc}}) or data frame
#'         (\code{\link{data.frame}}) object containing all of the parameters
#'         in the MCMC chain, with the burn-in discarded.
#'
#' @seealso \code{\link{readLog}}
#'
#' @examples
#'
readSingleLog <- function(filename, burnin=0.1, maxsamples=-1, as.mcmc=TRUE, burninAsSamples=FALSE) {
  if (!burninAsSamples && burnin > 1) {
      stop("Error: Burnin must be a proportion of samples between 0 and 1.", call. = FALSE)
  }

  if (burninAsSamples && burnin != round(burnin)) {
      stop("Error: Burnin must be an integer number of states.", call. = FALSE)
  }

  logfile <- read.table(filename, sep="\t", header=TRUE, nrows=maxsamples)
  n <- nrow(logfile)

  if (!burninAsSamples) {
      burnSamples <- floor(burnin*n)
  } else {
      burnSamples <- burnin+1
  }

  if (burnSamples >= n) {
      stop("Error: Discarding all samples in the log file.", call. = FALSE)
  }

  logfile <- logfile[burnSamples:n,]

  if (as.mcmc == TRUE) {
      if (is.null(logfile$Sample)) {
          if (is.null(logfile$state)) {
              logfile$Sample <- as.numeric(rownames(logfile))
          } else {
              logfile$Sample <- logfile$state
              logfile$state  <- NULL
          }
      }
    
      # Sometimes BEAST2 log files have an empty column at the end (lines end in tabs)
      # This needs to be removed before converting to mcmc object
      for (i in names(logfile)) {
           if (all(is.na(logfile[[i]]))) logfile[[i]] <- NULL
      }
    
      start <- logfile$Sample[1]
      thin  <- logfile$Sample[2]-logfile$Sample[1]
      rownames(logfile) <- logfile$Sample
      logfile$Sample    <- NULL

      return(coda::mcmc(logfile, start=start, thin=thin))
  } else {
      return(logfile)
  }
}

#' Read BEAST log files
#'
#' Read a single BEAST log file and return as a coda "mcmc" object.
#' If filenames contains more than one entry each log file is added as a separate chain and
#' a coda "mcmc.list" object is returned.
#'
#' @param filenames The name of the log file to read, or aternatively a vector or list of
#'        input log file names.
#' @inheritParams readSingleLog
#'
#' @return An "mcmc" object (\code{\link[coda]{mcmc}}) or data frame
#'         (\code{\link{data.frame}}) object containing all of the parameters
#'         in the MCMC chain, with the burn-in discarded. If more than one log file was
#'         given as input the output will be an "mcmc.list" object
#'         (\code{\link[coda]{mcmc.list}}) or a list of data frames
#'         (\code{\link{data.frame}}).
#'
#' @seealso \code{\link{readSingleLog}}
#'
#'
#' @examples
#'
#' @export
readLog <- function(filenames, burnin=0.1, maxsamples=-1, as.mcmc=TRUE, burninAsSamples=FALSE) {

  if (length(filenames) == 1) {
      return( readSingleLog(filenames, burnin, maxsamples, as.mcmc, burninAsSamples) )
  } else {

      loglist <- list()
      for (filename in filenames) {
          loglist[[filename]] <- readSingleLog(filename, burnin, maxsamples, as.mcmc, burninAsSamples)
      }

      if (as.mcmc == TRUE) {
          return (coda::mcmc.list(loglist))
      } else {
          return (loglist)
      }
  }

}

#' Save BEAST/Tracer log file
#'
#' What happens if it's a mcmc.list or a list of data frames?
#' What if input is just a data frame?
#'
#' @param logfile
#' @param filename
#' @param renumber Renumber states from 1 to nrow(logfile)
#' @param header Add a header with summary statistics to the file
#' @param model Add the model description from a logfile written in BEAST (not implemented)
#'
#' @export
writeLog <- function(logfile, filename, renumber=FALSE, header=TRUE, model="") {

    # Write header
    if (header) {
        headstr <- c(paste0("Log file saved by beastio v", packageVersion("beastio"), format(Sys.time(), " on %d %b %Y, %H:%M:%S")),
                    "(see https://github.com/laduplessis/beastio for more information)","","","summary:","")
        headstr <- c(headstr, capture.output( summary(logfile)))
        headstr <- paste("#",headstr)

        outfile <- file(filename, "wt")
        writeLines(headstr, outfile)
        close(outfile)
    }

    # Write log file
    if (renumber) {
        samples <- 0:(nrow(logfile)-1)
    } else {
        samples <- 1+thin(logfile)*(0:(nrow(logfile)-1))
    }
    logtable <- cbind(samples, data.frame(logfile))
    colnames(logtable)[1] <- "SAMPLE"
    write.table(logtable, filename, sep="\t", quote = FALSE, row.names = FALSE, col.names=TRUE, append = TRUE)

}


#' Extract a subset of parameters from a log file
#'
#' Extract all parameters from the log file that matches \code{pattern}.
#' e.g if par="R0" extract (R0s.1 R0s.2 R0s.3 etc.).
#' If the input log file contains more than one chain, the set of matching
#' parameters will be extracted from each chain.
#'
#'
#' @param logfile The logfile object to extract the subset from. Can be an
#'        "mcmc", "mcmc.list" or "data.frame" object, or a list of "data.frame"
#'        objects.
#' @param pattern The string to match with the names of the parameters to
#'        extract. This is treated as a regular expression.
#' @param start   If TRUE then only the start of parameter names are matched
#'        against pattern (\code{pattern} becomes "^<pattern>")
#'
#' @return The object returned is always of the same class as the input
#'         log file (\code{\link[coda]{mcmc}}, \code{\link[coda]{mcmc.list}},
#'         \code{\link{data.frame}} or a list of \code{\link{data.frame}}),
#'         but will only contain traces for the parameters matching
#'         \code{pattern}.
#'
#' @seealso \link{readLog}, \link{grep}
#'
#' @export
getLogFileSubset <- function(logfile, pattern, start=FALSE) {
  if (start == TRUE) {
      pattern <- paste0("^",pattern)
  }

  if (coda::is.mcmc(logfile) || coda::is.mcmc.list(logfile)) {
      return(logfile[, grepl(pattern, coda::varnames(logfile)), drop=TRUE])
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
#' @param data The object containing the MCMC sample - usually an object of class
#'        "mcmc" or "mcmc.list"
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
    stop("Error: Not an mcmc or mcmc.list object", call. = FALSE)
  }

}

#' Calculate the HPD interval and median and return in a single object
#'
#' Only works for mcmc and mcmc.list objects
#'
#' @param data The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#' @param prob Target probability content for the HPD interval (see ?HPDinterval)
#' @param ...  Extra parameters to be passed to the median function
#'
#' @export
getHPDMedian <- function(data, prob=0.95, ...) {

  t_hpd <- coda::HPDinterval(data, prob=prob)

  if (is.null(dim(data))) {
      t_med <- median(data)
  } else {
      t_med <- applyMCMC(data, median, ...)
  }

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
    stop("Error: Not an mcmc or mcmc.list object", call. = FALSE)
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

  if (!requireNamespace("boa", quietly = TRUE)) {
      stop("Error: Package \"boa\" needed for this function to work. Please install it.",
           call. = FALSE)
  }

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
  out <- apply(data, margin, getHPD.boa, ...)
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
      for (par in coda::varnames(data)) {
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
                  stop(paste0("Error: Chain ",i," has a different set of constant parameters to chain 1. Possibly some parameters got stuck because of a poor proposal distribution."),
                       call. = FALSE)
              }
          }
      }
      return(result)
  } else {
      stop("Error: Not an mcmc or mcmc.list object", call. = FALSE)
  }
}

#' The opposite to constantPars, that is returns all parameters
#' that have been sampled
#'
#' @param data The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#'
#' @export
sampledPars <- function(data) {
    return(setdiff(coda::varnames(data), constantPars(data)))
}


#' Overwrites coda::effectiveSize
#' Avoids numerical underflow by scaling parameters so mean is 1 before calculating
#' the effecive size
#'
#' @param data The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#' @param normalize If TRUE try to avoid numerical underflow, else just simply call
#'                  coda::effectiveSize
#'
#' @export
effectiveSize <- function(data, normalize=TRUE) {

  normalizeChain <- function(chain) {
    scaling <- 1/apply(chain, 2, mean)

    for (i in 1:coda::nvar(data)) {
      chain[, i] <- scaling[i] * chain[, i]
    }

    return(chain)
  }

  if (coda::is.mcmc.list(data)) {
    data <- lapply(data, normalizeChain)

  } else {
    data <- normalizeChain(as.mcmc(data))
  }

  return(coda::effectiveSize(data))
}


#' Plot barplot of ESS values
#'
#' @param data    The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#' @param ess     Output of running coda::effectiveSize(data) - this is optional and can be supplied
#'                to save time
#' @param cutoff  The cutoff ESS value for a parameter to be considered stationary (default = 200)
#' @param title   Title of the plot (x-label)
#' @param col1    Fill colour used for stationary parameters (ess >= cutoff)
#' @param col2    Fill colour used for non-stationary parameters (ess < cutoff)
#' @param border1 Border colour used for stationary parameters (ess >= cutoff)
#' @param border2 Border colour used for non-stationary parameters and labels (ess < cutoff)
#' @param labels  Label non-stationary parameters (ess < cutoff)
#' @param ...     Extra parameters passed to the \code{\link[gplots]{barplot2}}` function
#'
#' @export
plotESS <- function(data, ess=NULL, cutoff=200, title="",
                    col1="gray", col2="red2", border1="black", border2="darkred", labels=TRUE, ...) {

  if (!requireNamespace("gplots", quietly = TRUE)) {
    stop("Error: Package \"gplots\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (is.null(ess)) {
      ess <- effectiveSize(data)
  }

  nonstat      <- which(ess < cutoff)
  col          <- rep(col1, length(ess))
  col[nonstat] <- col2

  border          <- rep(border1, length(ess))
  border[nonstat] <- border2


  gplots::barplot2(ess, col=col, border=border, names.arg=NA, las=1,
                   ylab="ESS", xlab="", ...)
  mtext(side=1, text=title, line=1.5, cex.lab=0.8)
  abline(h=cutoff, lty=2, col=border2, lwd=1)

  if (length(nonstat) > 0 && labels) {
      nonstatx <- nonstat*1.2 - 0.5
      nonstatlab <- paste(names(ess)[nonstat], "=", round(ess[nonstat]))
      axis(3, at=nonstatx, labels=nonstatlab, pos=cutoff, col.axis=border2, lwd=0, las=2, cex.axis=0.8)
  }
}



#' Check if any parameters in an mcmc chain have an ESS value < cutoff
#' Constant (unsampled) parameters are ignored.
#'
#' @param data   The object containing the MCMC sample - usually of class "mcmc" or "mcmc.list"
#' @param cutoff The cutoff ESS value (default = 200)
#' @param value  If value == TRUE the parameters with ESS < cutoff and their ESS values are returned
#'               If value == FALSE only the column indices of those parameters are returned
#' @param ignored Parameters in the chain that should be ignored
#' @param plot   If plot == TRUE draw a barplot of the ESS values and mark parameters with ESS < cutoff
#' @param ...    Extra parameters to be passed to the plotESS() function
#'
#' @export
checkESS <- function(data, cutoff = 200, value = TRUE, ignored = c(), plot=FALSE, ...) {

  data <- as.mcmc(data)
  if (coda::is.mcmc(data) && !is.list(data)) {
      data <- data[, setdiff(sampledPars(data), ignored)]
      ess  <- effectiveSize(data)

      if (plot == TRUE) {
          plotESS(data, ess, cutoff = cutoff, ...)
      }

      if (value == TRUE) {
          return(ess[ess < cutoff])
      } else {
          return(which(ess < cutoff))
      }
  } else {
    stop("Error: Not an mcmc object", call. = FALSE)
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
