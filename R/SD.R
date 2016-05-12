"SD" <- 
function (x, na.rm = TRUE) 
{
    if (is.matrix(x)) 
        apply(x, 2, SD, na.rm = na.rm)
    else if (is.vector(x)) 
        sqrt(var(x, na.rm = na.rm,use="pair"))
    else if (is.data.frame(x)) 
        apply(x,2, SD, na.rm = na.rm)
    else sqrt(var(as.vector(x), na.rm = na.rm,use="pair"))
}
