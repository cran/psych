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
                      for (j in 1:p) { w[i,j] <- 1- (abs(i-j))^2/9 }
}}

                
weighted.prob <- w*prob
weighted.obser <- w*x
#wkappa <- 1-sum(weighted.obser)/sum(weighted.prob)
wpo <- sum(weighted.obser)
wpc <- sum(weighted.prob)
wkappa <- (wpo-wpc)/(1-wpc)

return(list(kappa=kappa,weighted.kappa = wkappa))
}

"cohen.kappa" <- function(x, w=NULL,n.obs=NULL,alpha=.05) {
cl <- match.call()
p <- dim(x)[1]
len <- p
if ((dim(x)[2] == p) ||(dim(x)[2]  < 3))   {result <- cohen.kappa1(x, w=w,n.obs=n.obs,alpha=alpha) } else {
nvar <- dim(x)[2]
ck <- matrix(NA,nvar,nvar)
if(!is.null(colnames(x)) ){colnames(ck) <- rownames(ck) <- colnames(x)} else {colnames(ck) <- rownames(ck) <- paste("R",1:nvar,sep="") }
diag(ck) <- 1
result <- list(cohen.kappa=ck)
k <- 2

for (i in 2:nvar ) {
   for (j in 1:(i-1) ) {
   x1 <- data.frame(x[,i],x[,j])
   x1 <- na.omit(x1)
    result[[paste(colnames(ck)[j],rownames(ck)[i])]] <- cohen.kappa1(x1, w=w,n.obs=n.obs,alpha=alpha)
    ck[i,j] <- result[[k]]$kappa
    ck[j,i] <- result[[k]]$weighted.kappa
    k <- k + 1
   }
   }
    result[[1]] <- ck
   }
 class(result) <- c("psych","kappa")
 return(result)
 }
   

"cohen.kappa1" <- function(x, w=NULL,n.obs=NULL,alpha=.05) {
cl <- match.call()
p <- dim(x)[1]
len <- p


if (dim(x)[2]!= p) {
x1 <- x[,1]
x2 <- x[,2]
if(is.factor(x1) ) {  #this gets around a problem of tabling numbers as characters (bad idea) vs. tabling characters (good idea)
  x1 <- as.character(x[,1])
 x2 <- as.character(x[,2])} else {
  x1 <- x[,1]
  x2 <- x[,2]}
labels <- levels(as.factor(cbind(x1,x2)))
 len <- length(labels)
 x <- matrix(0,ncol=len,nrow=len)
 colnames(x) <- rownames(x) <- labels
 x1f <- factor(x1,levels=labels)
 x2f <- factor(x2,levels=labels)
 x <- table(x1f,x2f)
 

# #for (item in 1:p) {x[x1[item],x2[item]] <- x[x1[item],x2[item]] +1}
 }

x <- as.matrix(x)
tot <- sum(x)
x <- x/tot  #convert to probabilities
rs <- rowSums(x)
cs <- colSums(x)
prob <- rs %*% t(cs)

po <- tr(x) 
pc <- tr(prob)
kappa <- (po-pc)/(1-pc)  #(model - data)/(1-model)
 
if(is.null(w)) { w <- matrix(0,ncol=len,nrow=len) 
                 w[] <- abs((col(w) - row(w)))^2 #squared weights
                 w <- 1 - w/(len-1)^2}   #1 - squared weights/k
colnames(w) <- rownames(w) <- colnames(x)             
weighted.prob <- w * prob
weighted.obser <- w * x
wpo <- sum(weighted.obser)
wpc <- sum(weighted.prob)
colw <- colSums(w*cs)
roww <- colSums(w*rs)
if((!is.null(n.obs)) & (tot==1))  tot <- n.obs
I <- diag(1,len,len)
Vark <-  (1/(tot*(1-pc)^4))*    (tr(x * (I * (1-pc)   - (rs %+% t(cs ))*(1-po))^2 )  + (1-po)^2 * (sum(x * (cs %+% t(rs ))^2) - tr(x * (cs %+% t(rs ))^2))  -(po*pc - 2*pc +po)^2    )  
Varkw <-  (1/(tot*(1-wpc)^4))*  (sum(x * (w * (1-wpc)- (colw %+% t(roww ))*(1-wpo))^2 ) -(wpo*wpc - 2*wpc +wpo)^2    )   
if(tr(w) > 0) {wkappa <- (wpo-wpc)/(1-wpc) } else { wkappa <- 1- wpo/wpc}
if((!is.null(n.obs)) & (tot==1))  tot <- n.obs
bounds <- matrix(NA,2,3)
colnames(bounds) <- c("lower","estimate","upper")
rownames(bounds) <- c("unweighted kappa","weighted kappa")
bounds[1,2] <- kappa
bounds[2,2] <- wkappa
bounds[1,1] <- kappa + qnorm(alpha/2) * sqrt(Vark) 
bounds[1,3] <- kappa - qnorm(alpha/2) * sqrt(Vark)
bounds[2,1] <- wkappa + qnorm(alpha/2) * sqrt(Varkw) 
bounds[2,3] <- wkappa - qnorm(alpha/2) * sqrt(Varkw)
result <- list(kappa=kappa,weighted.kappa = wkappa,n.obs=tot,agree=x,weight=w,var.kappa =Vark, var.weighted = Varkw,confid=bounds,plevel=alpha,Call=cl)
class(result) <- c("psych","kappa")
return(result)
}