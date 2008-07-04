#corrected may 7, 2007 
skew <- 
function (x, na.rm = TRUE) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
        sdx <- SD(x,na.rm=na.rm)
        mx <- mean(x)
        skewer <- sum((x - mx)^3)/(length(x) * SD(x)^3) 
        } else {
    
    skewer <- rep(NA,dim(x)[2])
    if (is.matrix(x)) {mx <- colMeans(x,na.rm=na.rm)} else {mx <- mean(x,na.rm=na.rm)}
    sdx <- SD(x,na.rm=na.rm)
    for (i in 1:dim(x)[2]) {
    skewer[i] <- sum((x[,i] - mx[i])^3,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i]))) * sdx[i]^3)
            }
    }
    return(skewer)
}
