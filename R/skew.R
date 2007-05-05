#revised October 24 to run slightly faster by taking mean outside of sum
skew <- 
function (x, na.rm = TRUE) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
         
        mx <- mean(x)
        skewer <- sum((x - mx)^3)/(length(x) * sd(x)^3) 
        } else {
    
    skewer <- rep(NA,dim(x)[2])
    mx <- colMeans(x,na.rm=na.rm)
    
    for (i in 1:dim(x)[2]) {
    skewer[i]<- sum((x[,i] - mx[i])^3,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i])) * sd(x[,i], 
            na.rm = na.rm)^3))
            }
    }
    return(skewer)
}
