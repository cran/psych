#adapted from John Fox's Polychor
#polyc does all the work but does not work in cases of incomplete tables
#thus, the polychor function is used
#moved these first two function out of the polyc function in the hope that they will be compiled just once and perhaps get a speed increase
#doesn't seem to make a difference although it does make the code a bit easier to read

 polyBinBvn <- function (rho,rc,cc)    #adapted from John Fox's polychor
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


 polyF <- function(rho,rc,cc,tab) { 
      P <- polyBinBvn(rho, rc, cc) 
       -sum(tab * log(P)) }  #the  criterion to be minimized

"polyc" <- 
function(x,y=NULL,taux,tauy,global=TRUE) {

  tab <- table(x,y)  
  tot <- sum(tab)
  tab <- tab/tot
  
 if(global) { rho <- optimize(polyF,interval=c(-1,1),rc=taux, cc=tauy,tab)#this uses the global taux and tauy
       } else { #use item row and column information for this pair, rather than global values
      #this seems to match the polycor function
      #the next five lines are taken directly from John Fox's polycor function
      zerorows <- apply(tab, 1, function(x) all(x == 0))
	zerocols <- apply(tab, 2, function(x) all(x == 0))
	zr <- sum(zerorows)
	#if (0 < zr) warning(paste(zr, " row", suffix <- if(zr == 1) "" else "s"," with zero marginal", suffix," removed", sep=""))
	zc <- sum(zerocols)
	#if (0 < zc) warning(paste(zc, " column", suffix <- if(zc == 1) "" else "s", " with zero marginal", suffix, " removed", sep=""))
	tab <- tab[!zerorows, ,drop=FALSE]  
	tab <- tab[, !zerocols, drop=FALSE] 
    csum <- colSums(tab)
    rsum <- rowSums(tab)
     cc <-  qnorm(cumsum(csum))[-length(csum)]
     rc <-  qnorm(cumsum(rsum))[-length(rsum)]
    rho <- optimize(polyF,interval=c(-1,1),rc=rc, cc=cc,tab)
        }
  result <- list(rho=rho$minimum,objective=rho$objective)
  return(result)
  }
  
  

#Basically just is used to find the thresholds and then does the polychoric r for a matrix
"polychoric" <- 
function(x,smooth=TRUE,global=TRUE,polycor=FALSE,ML = FALSE, std.err = FALSE) {
if(!require(mvtnorm) ) {stop("I am sorry, you must have mvtnorm installed to use polychoric")}
if(polycor && (!require(polycor))) {warning ("I am sorry, you must have  polycor installed to use polychoric with the polycor option")
 polycor <- FALSE}
 cl <- match.call() 
nvar <- dim(x)[2]
nsub <- dim(x)[1]
x <- as.matrix(x)
xt <- table(x)
#nvalues <- length(xt)  #find the number of response alternatives 
nvalues <- max(x,na.rm=TRUE) - min(x,na.rm=TRUE) + 1
if(nvalues > 8) stop("You have more than 8 categories for your items, polychoric is probably not needed")
xmin <- apply(x,2,function(x) min(x,na.rm=TRUE))  #allow for different minima
x <- t(t(x) - xmin +1)  #all numbers now go from 1 to nvalues
#xfreq <- apply(x- xmin + 1,2,tabulate,nbins=nvalues)
xfreq <- apply(x,2,tabulate,nbins=nvalues)
n.obs <- colSums(xfreq)
xfreq <- t(t(xfreq)/n.obs)
tau <- qnorm(apply(xfreq,2,cumsum))[1:(nvalues-1),]  #these are the normal values of the cuts
if(!is.matrix(tau)) tau <- matrix(tau,ncol=nvar)
rownames(tau) <- names(xt)[1:(nvalues-1)]
colnames(tau) <- colnames(x)

 
mat <- matrix(0,nvar,nvar)
colnames(mat) <- rownames(mat) <- colnames(x)
x <- x - min(x,na.rm=TRUE) +1  #this is essential to get the table function to order the data correctly

#cat("\n Finding Polychoric correlations\n" )
for (i in 2:nvar) {
progressBar(i^2/2,nvar^2/2,"Polychoric")
  for (j in 1:(i-1)) { 
 if(t(!is.na(x[,i]))%*% (!is.na(x[,j]))  > 2 ) {
if(!polycor) {
poly <-  polyc(x[,i],x[,j],tau[,i],tau[,j],global=global)  
mat[i,j] <- mat[j,i] <- poly$rho } else {
#To use John Fox's version which requires the polycor package
  poly <- polychor(x[,i],x[,j],ML = ML,std.err = std.err)  #uses John Fox's function
  mat[i,j] <- mat[j,i] <- poly }
  } else {mat[i,j] <- mat[j,i] <- NA}
   }
   } 
   diag(mat) <- 1
  if(any(is.na(mat))) {warning("some correlations are missing, smoothing turned off")
                        smooth <- FALSE}
                     
 if(smooth) {mat <- cor.smooth(mat) }
 tau <- t(tau)
  result <- list(rho = mat,tau = tau,n.obs=nsub,Call=cl) 
   flush(stdout())
 cat("\n") #put in to clear the progress bar
 flush(stdout())
 class(result) <- c("psych","poly")
  return(result) 
  }


