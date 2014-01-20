#Faster Polychoric uses tableF (a function to speed up 2 way tables of integers
#first, we introduce a new function to find integer tables without error checking

#if all data are integer then
#tableF is the fast version of table  
#it does no error checking and works only for two dimensional integer data
tableF <-
function(x,y) {
minx <- min(x,na.rm=TRUE)
maxx <- max(x,na.rm=TRUE)
miny <- min(y,na.rm=TRUE)
maxy <- max(y,na.rm=TRUE)
maxxy <- (maxx+(minx==0))*(maxy+(miny==0))
dims=c(maxx + 1 - min(1,minx),maxy+1 - min(1,minx))
bin <- x - minx+ (y-miny)*(dims[1])+ max(1,minx)
ans <- matrix(tabulate(bin,maxxy),dims)
ans
}

#adapted from John Fox's Polychor
#polyc does all the work but does not work in cases of incomplete tables
#thus, the polychor function is used
#moved these first two function out of the polyc function in the hope that they will be compiled just once and perhaps get a speed increase
#doesn't seem to make a difference although it does make the code a bit easier to read
#polychoric.mc  is added while we test it versus polychoric

 polyBinBvn.old <- function (rho,rc,cc)    #adapted from John Fox's polychor
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

 polyBinBvn<- function (rho,rc,cc)  {   #adapted from John Fox's polychor
      #recognizes that we don't need to calculate all cells because of degrees of freedom                               
 if (min(rc) < -9999) rc <- rc[-1]          
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
    for (i in 1:(nr-1)) {
        for (j in 1:(nc-1)) {
            P[i, j] <- pmvnorm(lower = c(row.cuts[i], col.cuts[j]), 
                upper = c(row.cuts[i + 1], col.cuts[j + 1]), 
                corr = R)   #should we specify the algorithm to TVPACK or Miwa
        }}
    P[1,nc] <- pnorm(rc[1]) - sum(P[1,1:(nc-1)] )
    P[nr,1] <- pnorm(cc[1]) - sum(P[1:(nr-1),1] )
    if(nr >2) {for (i in (2:(nr-1))) {P[i,nc] <- pnorm(rc[i]) -pnorm(rc[i-1])- sum(P[i,1:(nc-1)] ) }}
    if(nc >2) {for (j in (2:(nc-1))) {P[nr,j] <- pnorm(cc[j]) - pnorm(cc[j-1])-sum(P[1:(nr-1),j] ) }}
    if(nc>1)  P[nr,nc] <- 1- pnorm(rc[nr-1]) - sum(P[nr,1:(nc-1)]) 
    P   #the estimated n x n predicted by rho, rc, cc
}


 polyF <- function(rho,rc,cc,tab) { 
      P <- polyBinBvn(rho, rc, cc) 
       -sum(tab * log(P)) }  #the  criterion to be minimized
       
       
"wtd.table" <- function(x,y,weight) {
tab <- tapply(weight,list(x,y),sum,na.rm=TRUE,simplify=TRUE) #taken from questionr:wtd.table
tab[is.na(tab)] <- 0
return(tab)
}      
  
  "polyc" <-     #uses the tableF function instead of table
function(x,y=NULL,taux,tauy,global=TRUE,weight=NULL) {
  if(is.null(weight )) {tab <- tableF(x,y) }   else {tab <- wtd.table(x,y,weight)} 
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
  

"polychoric" <- 
function(x,smooth=TRUE,global=TRUE,polycor=FALSE,ML = FALSE, std.err = FALSE,weight=NULL,progress=TRUE,na.rm=TRUE,delete=TRUE) {
if(!require(parallel)) {message("polychoric requires the parallel package.")}
#declare these next two functions to be local inside of polychoric

myfun <- function(x,i,j) {polyc(x[,i],x[,j],tau[,i],tau[,j],global=global,weight=weight) }

matpLower <- function(x,nvar) {
k <- 1
il <- vector()
jl <- vector()
for(i in 2:nvar) {for (j in 1:(i-1)) {
   il[k] <- i
   jl [k] <- j
   k<- k+1}
   }
poly <- mcmapply(function(i,j) myfun(x,i,j) , il,jl) 
#now make it a matrix
mat <- diag(nvar)
mat[upper.tri(mat)] <- as.numeric(poly[1,]) #first row of poly is correlation, 2nd the fit
mat <- t(mat) + mat
diag(mat) <- 1
return(mat)
}


if(!require(mvtnorm) ) {stop("I am sorry, you must have mvtnorm installed to use polychoric")}
if(polycor && (!require(polycor))) {warning ("I am sorry, you must have  polycor installed to use polychoric with the polycor option")
 polycor <- FALSE}
 if(!is.null(weight)) {if(length(weight) !=nrow(x)) {stop("length of the weight vector must match the number of cases")}}
 cl <- match.call() 
nvar <- dim(x)[2]
nsub <- dim(x)[1]
x <- as.matrix(x)
xt <- table(x)
#nvalues <- length(xt)  #find the number of response alternatives 
nvalues <- max(x,na.rm=TRUE) - min(x,na.rm=TRUE) + 1
if(nvalues > 8) stop("You have more than 8 categories for your items, polychoric is probably not needed")
 #first delete any bad cases
  item.var <- apply(x,2,sd,na.rm=na.rm)
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0)  & delete) {
            for (baddy in 1:length(bad)) {warning( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
            nvar <- nvar - length(bad)
             }
xmin <- apply(x,2,function(x) min(x,na.rm=TRUE))  #allow for different minima
x <- t(t(x) - xmin +1)  #all numbers now go from 1 to nvalues
xmax <- apply(x,2,function(x)  max(x,na.rm=TRUE)) #check for different maxima
if (min(xmax) != max(xmax)) {global <- FALSE
                      message("The items do not have an equal number of response alternatives, global set to FALSE")}
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

mat <- matpLower(x,nvar)  #the local copy has the extra paremeters   #do the multicore version

   

  if(any(is.na(mat))) {warning("some correlations are missing, smoothing turned off")
                        smooth <- FALSE}
                     
 if(smooth) {mat <- cor.smooth(mat) }
 colnames(mat) <- rownames(mat) <- colnames(x)
 tau <- t(tau)
  result <- list(rho = mat,tau = tau,n.obs=nsub,Call=cl) 
  
 class(result) <- c("psych","poly")
  return(result) 
  }






#draft version to do one item at a time 
#not public 
#use polychor from John Fox to do the same
"polytab" <- 
function(tab) {
  tot <- sum(tab)
  tab <- tab/tot
  
  
 #use item row and column information for this pair, rather than global values
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
     
  result <- list(rho=rho$minimum,objective=rho$objective)
  return(result)
  }

