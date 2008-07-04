"rescale" <-  function(x,mean=100,sd=15,df=TRUE) {if(df) {x <- data.frame(t(t(scale(x))*sd+mean))
} else {x <- t( t(scale(x))*sd +mean)}
return(x)
}
#corrected April 3 to properly do matrix addition