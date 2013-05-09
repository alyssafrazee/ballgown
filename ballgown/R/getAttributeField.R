getAttributeField <- function (x, field, attrsep = "; ") 
{
    s = strsplit(x, split = attrsep, fixed = TRUE)
    sapply(s, function(atts) {
        a = strsplit(atts, split = " ", fixed = TRUE)
        m = match(field, sapply(a, "[", 1))
        if (!is.na(m)) {
            rv = a[[m]][2]
        }
        else {
            rv = as.character(NA)
        }
        return(rv)
    })
}