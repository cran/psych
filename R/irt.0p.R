"irt.0p" <- function(items) {
possible <- dim(items)[2]
raw <- rowMeans(items,na.rm=TRUE) 
ave <- raw
valid <- rowSums(!is.na(items))
ave[(!is.na(ave))&(ave<.0001)] <- .5/(possible)
ave[(!is.na(ave))&(ave > .9999)] <-  (possible-.5)/possible
theta <- -log((1/ave) -1)
irt.0p <- matrix(c(raw,theta,valid),ncol=3)
colnames(irt.0p ) <- c("raw","theta0","valid")
return(irt.0p)
}