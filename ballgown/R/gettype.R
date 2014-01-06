#' extract type of expression measurement from column names of expression tables
#'
#' @param x column name from which measurement type should be extracted
#' @return measurement type (e.g., \code{'FPKM'} or \code{'rcov'})
#' @details Mostly an internal function, but we imagine it could be useful to users as well, so we export it.  It returns the last section of a string, where sections are delimited by '.' -- this is useful in ballgown columns of the \code{texpr}, \code{eexpr}, \code{iexpr}, and \code{gexpr} tables are named using the scheme \code{<MEASUREMENT>.<SAMPLENAME>}.
#' @export

gettype = function(x) strsplit(x, split="\\.")[[1]][1]