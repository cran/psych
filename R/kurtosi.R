"kurtosi" <- 
function (x, na.rm = TRUE) 
{
    if (length(dim(x)) == 0) {
        if (na.rm) {
            x <- x[!is.na(x)]
        			}
       if (is.matrix(x) ) { mx <- colMeans(x,na.rm=na.rm)} else {mx <- mean(x,na.rm=na.rm)}
       
         sdx <- SD(x,na.rm=na.rm)
        kurt <- sum((x - mx)^4)/(length(x) * SD(x)^4)  -3
        } else {
    
    kurt <- rep(NA,dim(x)[2])
  #  mx <- mean(x,na.rm=na.rm)
  if (is.matrix(x) ) { mx <- colMeans(x,na.rm=na.rm)} else {mx <- mean(x,na.rm=na.rm)}
       
    sdx <- SD(x,na.rm=na.rm)
    for (i in 1:dim(x)[2]) {
    kurt[i] <- sum((x[,i] - mx[i])^4,  na.rm = na.rm)/((length(x[,i]) - sum(is.na(x[,i]))) * sdx[i]^4)  -3
            }
    }
    return(kurt)
}
