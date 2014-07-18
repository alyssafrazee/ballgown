ss = function(x, pattern, slot=1, ...){
    sapply(strsplit(x,pattern,...), "[", slot)  
} 
