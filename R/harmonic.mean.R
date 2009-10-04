"harmonic.mean" <- 
function(x,na.rm=TRUE) { 1/(mean(1/x,na.rm=na.rm)) }