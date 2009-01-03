"sim.dichot" <-
function (nvar = 72, nsub = 500, circum = FALSE, xloading = 0.6, 
    yloading = 0.6, gloading = 0, xbias = 0, ybias = 0,  
    low = 0, high = 0) 
{
    avloading <- (xloading + yloading)/2
    errorweight <- sqrt(1 - (avloading^2 + gloading^2))
    g <- rnorm(nsub)
    truex <- rnorm(nsub) * xloading + xbias
    truey <- rnorm(nsub) * yloading + ybias
    if (circum) {
        radia <- seq(0, 2 * pi, len = nvar + 1)
        rad <- radia[which(radia < 2 * pi)]
    }
    else rad <- c(rep(0, nvar/4), rep(pi/2, nvar/4), rep(pi, 
        nvar/4), rep(3 * pi/2, nvar/4)) 
    error <- matrix(rnorm(nsub * (nvar)), nsub)
    trueitem <- outer(truex, cos(rad)) + outer(truey, sin(rad))
    item <- gloading * g + trueitem + errorweight * error
    
   
    nvar2 <- nvar/2 
     iteml <- (item[,(1:nvar2)*2 -1] >= low)
     itemh <- (item[,(1:nvar2)*2] >= high)
     item <- cbind(iteml,itemh)+0
    

    return(item)
}
#revised October 2 to make difficulty and direction of factor loading unconfounded