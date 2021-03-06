% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/log_utils.R
\name{getHPD.boa}
\alias{getHPD.boa}
\title{Calculate HPD of a parameter}
\usage{
getHPD.boa(data, alpha = 0.05, includeMedian = TRUE)
}
\arguments{
\item{data}{The samples from the posterior.}

\item{alpha}{The confidence level.}

\item{includleMedian}{If FALSE only return the HPD limits, and not the median as well.}
}
\value{
c(lower, median, upper)
}
\description{
Get HPD of a posterior sample
}
\details{
Uses Chen and Shao algorithm as implemented in boa package.
}
