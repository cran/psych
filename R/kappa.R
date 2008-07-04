"wkappa" <- 
function(x,w=NULL) {
p <- dim(x)[2]
if (dim(x)[1]!= p) x <- table(x[,1],x[,2])
x <- as.matrix(x)
tot <- sum(x)
x <- x/tot  #convert to probabilities
rs <- rowSums(x)
cs <- colSums(x)
prob <- rs %*% t(cs)

po <- tr(x)
pc <- tr(prob)
kappa <- (po-pc)/(1-pc)
 
if(is.null(w)) { w <- matrix(0,ncol=p,nrow=p) 
                for (i in 1:p) {
                      for (j in 1:p) { w[i,j] <- .5 ^(abs(i-j)) }
}}

                
weighted.prob <- w*prob
weighted.obser <- w*x
#wkappa <- 1-sum(weighted.obser)/sum(weighted.prob)
wpo <- sum(weighted.obser)
wpc <- sum(weighted.prob)
wkappa <- (wpo-wpc)/(1-wpc)

return(list(kappa=kappa,weighted.kappa = wkappa))
}
