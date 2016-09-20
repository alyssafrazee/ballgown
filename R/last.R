#' get the last element
#'
#' @param x anything you can call \code{tail} on (vector, data frame, etc.)
#' 
#' @return the last element of \code{x}
#' 
#' @details this function is made of several thousand lines of complex code, 
#'    so be sure to read it carefully.
#' 
#' @author Alyssa Frazee
#' 
#' @export
#' @examples
#' last(c('h', 'e', 'l', 'l', 'o'))

last = function(x) return(tail(x, n=1))