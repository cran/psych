#adapted from John Fox's Polychor
#polyc does all the work but does not work in cases of incomplete tables
#thus, the polychor function is used
"polyc" <- 
function(x,y=NULL,taux,tauy,correct=TRUE) {
 binBvn <- function (rho,rc,cc)    #adapted from John Fox's polychor
{ if (min(rc) < -9999) rc <- rc[-1]
  if (min(cc) < - 9999) cc <- cc[-1]
  if (max(rc) > 9999) rc <- rc[-length(rc)]
  if (max(cc)  > 99999) cc <- cc[-length(cc)]
  row.cuts <- c(-Inf,rc,Inf)
  col.cuts <- c(-Inf,cc,Inf)
  nr <- length(rc) + 1
  nc <- length(cc) + 1


    P <- matrix(0, nr,nc)
    R <- matrix(c(1,rho,rho,1),2,2)
   # diag(R) <- 1
    for (i in 1:nr) {
        for (j in 1:nc) {
            P[i, j] <- pmvnorm(lower = c(row.cuts[i], col.cuts[j]), 
                upper = c(row.cuts[i + 1], col.cuts[j + 1]), 
                corr = R)   #should we specify the algorithm to TVPACK or Miwa
        }}
    P   #the estimated n x n predicted by rho, rc, cc
}
 f <- function(rho,rc,cc) { 
      P <- binBvn(rho, rc, cc) 
       -sum(tab * log(P)) }  #the ML criterion to be minimized
      
  tab <- table(x,y)  
  tot <- sum(tab)
  tab <- tab/tot
  
  rho <- optimize(f,interval=c(-1,1),rc=taux, cc=tauy)
  result <- list(rho=rho$minimum,objective=rho$objective)
  return(result)
  }
  
  

#Basically just is used to find the thresholds and then does the polychoric r for a matrix
"polychoric" <- 
function(x,polycor=FALSE,ML = FALSE, std.err = FALSE) {
if(!require(mvtnorm) ) {stop("I am sorry, you must have mvtnorm installed to use polychoric")}
if(polycor && (!require(polycor))) {warning ("I am sorry, you must have  polycor installed to use polychoric with the polycor option")
 polycor <- FALSE}
 cl <- match.call() 
nvar <- dim(x)[2]
nsub <- dim(x)[1]
x <-as.matrix(x)
xt <- table(x)
nvalues <- length(xt)  #find the number of response alternatives 
if(nvalues > 10) stop("You have more than 10 categories for your items, polychoric is probably not needed")
xmin <- min(x,na.rm=TRUE)
xfreq <- apply(x- xmin+ 1,2,tabulate,nbins=nvalues)
n.obs <- colSums(xfreq)
xfreq <- t(t(xfreq)/n.obs)
tau <- qnorm(apply(xfreq,2,cumsum))[1:(nvalues-1),]  #these are the normal values of the cuts
if(!is.matrix(tau)) tau <- matrix(tau,ncol=nvar)
rownames(tau) <- names(xt)[1:(nvalues-1)]
colnames(tau) <- colnames(x)

 
mat <- matrix(0,nvar,nvar)
colnames(mat) <- rownames(mat) <- colnames(x)
x <- x - min(x,na.rm=TRUE) +1  #this is essential to get the table function to order the data correctly

for (i in 2:nvar) {
  for (j in 1:(i-1)) {
if(!polycor) {
poly <-  polyc(x[,i],x[,j],tau[,i],tau[,j])  
mat[i,j] <- mat[j,i] <- poly$rho } else {
#To use John Fox's version which requires the polycor package
  poly <- polychor(x[,i],x[,j],ML = ML,std.err = std.err)  #uses John Fox's function
  mat[i,j] <- mat[j,i] <- poly }
   }
   }
   diag(mat) <- 1
  tau <- t(tau)
  result <- list(rho = mat,tau = tau,n.obs=nsub,Call=cl) 
 class(result) <- c("psych","poly")
  return(result) 
  }
  
