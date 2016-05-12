"paired.r" <-
function(xy,xz,yz=NULL,n,n2=NULL,twotailed=TRUE) {
cl <- match.call() 
if (!is.null(yz)) {
       diff <- xy-xz
       determin=1-xy*xy - xz*xz - yz*yz + 2*xy*xz*yz
       av=(xy+xz)/2
       cube= (1-yz)*(1-yz)*(1-yz)
       t2 = diff * sqrt((n-1)*(1+yz)/(((2*(n-1)/(n-3))*determin+av*av*cube)))
       p <- pt(abs(t2),n-2,lower.tail=FALSE)
        if(twotailed) p <- 2*p
       value <- list(test="test of difference between two correlated  correlations",t=t2,p=p,Call=cl)
       } else {
        xy.z <- 0.5*log((1+xy)/(1-xy))
        xz.z <- 0.5*log((1+xz)/(1-xz))
        if(is.null(n2)) n2 <- n
        se.diff.r <- sqrt(1/(n-3) + 1/(n2-3))
        diff <- xy.z - xz.z
        z <- abs(diff/se.diff.r)
         p <- (1-pnorm(z ))
          if(twotailed) p <- 2*p
        value <- list(test="test of difference between two independent correlations",z=z,p=p,Call=cl)
     	}
     	class(value) <- c("psych","paired.r")
     	return(value)
     	}
     	

  
  