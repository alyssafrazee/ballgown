closestColor = function(x, colscale){
    choices = rev(heat.colors(length(colscale)))
    diffs = abs(x-colscale)
    return(choices[which.min(diffs)])
}
