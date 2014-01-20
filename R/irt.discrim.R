"irt.discrim" <- 
function(item.diff,theta,items) {
#find the item discrimination parameter  (beta)
#find the item discrimination parameter -- optimized by item.discrim
 irt.item.discrim <- function(x,diff,theta,scores) {
  fit <- -1*(log(scores/(1+exp(x*(diff-theta))) + (1-scores)/(1+exp(x*(theta-diff)))))
  mean(fit,na.rm=TRUE)
  }
 nitems <- length(item.diff)
 discrim <- matrix(NaN,nitems,2)
 for (i in 1:nitems) {
    item.fit <- optimize(irt.item.discrim,c(-5,5),diff=item.diff[i],theta=theta,scores = items[,i])
   discrim[i,1] <- item.fit$minimum
   discrim[i,2] <- item.fit$objective}
  irt.discrim <- discrim
 }	