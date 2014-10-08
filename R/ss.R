ss = function(x, pattern, slot=1, ...){
    sapply(strsplit(x,pattern,...), "[", slot)  
}

strip_quotes = function(x){
    substr(x, 2, nchar(x)-1)
}
