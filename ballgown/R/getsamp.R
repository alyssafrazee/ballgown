#' extract sample name from column names of expression tables
#'
#' @param x column name from which sample should be extracted
#' @return sample name
#' @details Mostly an internal function, but we imagine it could be useful to users as well, so we export it.  It returns the last section of a string, where sections are delimited by '.' -- this is useful in ballgown columns of the \code{texpr}, \code{eexpr}, \code{iexpr}, and \code{gexpr} tables are named using the scheme \code{<MEASUREMENT>.<SAMPLENAME>}.  On a related note, you samples shouldn't be named using dots.
#' @export

getsamp = function(x) last(strsplit(x, split="\\.")[[1]])