#' Generate a sequence of months
#'
#' This function will return a sequence of the months between
#' the dates \code{from} and \code{to}.
#' The dates in the sequence that is returned are on the first day of
#' each month.
#'
#' @param from,to Start and end of the date sequence (these should be of
#'        class \code{\link{Date}}.
#' @param floorFrom If TRUE the first element in the sequence is the start of
#'   the first month before \code{from}, otherwise \code{from}.
#' @param ceilTo If TRUE the last element in the sequence is the start of
#'   the first month after \code{to}, otherwise \code{to}.
#'
#' @return
#'
#' @seealso
#'
#' @export
getMonths <- function(from, to, floorFrom=TRUE, ceilTo=TRUE) {

  start <- lubridate::floor_date(from, unit="month")
  end   <- lubridate::ceiling_date(to, unit="month")
  len   <- (start %--% end)/months(1)

  dateseq  <- start %m+% months(0:len)

  if (floorFrom != TRUE) {
      dateseq[1] <- from
  }
  if (ceilTo != TRUE) {
      dateseq[length(dateseq)] <- to
  }

  return(dateseq)
}


#' Generate a sequence of weeks
#'
#' This function will return all the weeks between \code{from} and \code{to}.
#' The dates in the sequence are the first day of each week.
#' The first day of the week can be set using the \code{week_start} parameter.
#'
#' @param from,to Start and end of the date sequence (these should be date
#'   objects).
#' @param week_start Specify the week start, default is Sunday (7)
#' @param floorFrom If TRUE the first element in the sequence is the start of
#'   the first month before \code{from}, otherwise \code{from}.
#' @param ceilTo If TRUE the last element in the sequence is the start of
#'   the first month after \code{to}, otherwise \code{to}.
#'
#' @export
getWeeks <- function(from, to, week_start = getOption("lubridate.week.start", 7),
                      floorFrom=FALSE, ceilTo=FALSE) {

  start <- lubridate::floor_date(from, unit="week", week_start=week_start)
  end   <- lubridate::ceiling_date(to, unit="week", week_start=week_start)
  len   <- (start %--% end)/weeks(1)

  dateseq  <- start + weeks(0:len)

  if (floorFrom != TRUE) {
    dateseq[1] <- from
  }

  if (ceilTo != TRUE) {
    if (dateseq[length(dateseq)-1] == to) {
      dateseq <- dataset[-length(dateseq)]
    } else {
      dateseq[length(dateseq)] <- to
    }
  }

  return(dateseq)
}
