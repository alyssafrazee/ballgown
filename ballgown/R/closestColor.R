#' find color corresponding to an integer
#'
#' @param x number 
#' @param colscale vector representing range of numbers the color scale is representing
#' @return color (from \code{heat.colors}) that most closely matches \code{x} in the given scale
#' @details internal function for \code{plotTranscripts} - not intended for direct use
#' @seealso \link{\code{plotTranscripts}}

closestColor = function(x, colscale){
	choices = rev(heat.colors(length(colscale)))
	diffs = abs(x-colscale)
	return(choices[which.min(diffs)])
}
