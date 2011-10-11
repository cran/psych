"harmonic.mean" <- 
function(x,na.rm=TRUE) {if (is.null(nrow(x))) {1/mean(1/x,na.rm=na.rm) } else {
 1/(apply(1/x,2,mean,na.rm=na.rm))} }