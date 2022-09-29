#Faster Polychoric uses tableF (a function to speed up 2 way tables of integers
#first, we introduce a new function to find integer tables without error checking

#if all data are integer then
#tableF is the fast version of table  
#it does no error checking and works only for two dimensional integer data
tableF <-
function(x,y) {
minx <- min(x,na.rm=TRUE)   #these probably could be found just once
maxx <- max(x,na.rm=TRUE)
miny <- min(y,na.rm=TRUE)
maxy <- max(y,na.rm=TRUE)
maxxy <- (maxx+(minx==0))*(maxy+(miny==0))
dims=c(maxx + 1 - min(1,minx),maxy+1 - min(1,minx))
bin <- x - minx+ (y-miny)*(dims[1])+ max(1,minx)
ans <- matrix(tabulate(bin,maxxy),dims)
ans
}


#perhaps even faster, but more importantly does not drop categories  - probably needs to be passed both x and y min and max
tableFast <-   #revised and preferred, but requires specifying the min and max
function(x,y,minx,maxx,miny,maxy) { #y and x can have separate min and max in the case of polydi,normally they are the same
maxxy <-  (maxx+(minx==0))*(maxy+(minx==0))
bin <- x-minx + (y-minx) *maxx+ 1
dims=c(maxx + 1 - min(1,minx),maxy+1 - min(1,miny))
ans <- matrix(tabulate(bin,maxxy),dims)
ans
}

#adapted from John Fox's Polychor
#polyc does all the work but does not work in cases of incomplete tables
#thus, the polychor function is used
#moved these first two function out of the polyc function in the hope that they will be compiled just once and perhaps get a speed increase
#doesn't seem to make a difference although it does make the code a bit easier to read
#polychoric.mc  is added while we test it versus polychoric

#  polyBinBvn.old <- function (rho,rc,cc)    #adapted from John Fox's polychor
# { if (min(rc) < -9999) rc <- rc[-1]          
#   if (min(cc) < - 9999) cc <- cc[-1]
#   if (max(rc) > 9999) rc <- rc[-length(rc)]
#   if (max(cc)  > 99999) cc <- cc[-length(cc)]
#   row.cuts <- c(-Inf,rc,Inf)
#   col.cuts <- c(-Inf,cc,Inf)
#   nr <- length(rc) + 1
#   nc <- length(cc) + 1
# 
# 
#     P <- matrix(0, nr,nc)
#     R <- matrix(c(1,rho,rho,1),2,2)
#    # diag(R) <- 1
#     for (i in 1:nr) {
#         for (j in 1:nc) {
#             P[i, j] <- pmvnorm(lower = c(row.cuts[i], col.cuts[j]), 
#                 upper = c(row.cuts[i + 1], col.cuts[j + 1]), 
#                 corr = R)   #should we specify the algorithm to TVPACK or Miwa
#         }}
#     P   #the estimated n x n predicted by rho, rc, cc
# }

 polyBinBvn<- function (rho,rc,cc)  {   #adapted from John Fox's polychor
      #but recognizes that we don't need to calculate all cells because of degrees of freedom                               
#  if ( min(rc,na.rm=TRUE) < -9999) rc <- rc[-1]          
#  if ( min(cc,na.rm=TRUE) < - 9999) cc <- cc[-1]
#  if (max(rc,na.rm=TRUE) > 9999)  rc <- rc[-length(rc)]
#  if (max(cc,na.rm=TRUE)  > 9999) cc <- cc[-length(cc)]
  row.cuts <- c(-Inf,rc,Inf)
  col.cuts <- c(-Inf,cc,Inf)
 # nr <- length(rc) + 1
 # nc <- length(cc) + 1   #replaced with next two lines 9/8/14
   nr <- length(row.cuts) -1
   nc <- length(col.cuts) -1   

    P <- matrix(0, nr,nc)
    R <- matrix(c(1,rho,rho,1),2,2)
   # diag(R) <- 1
    for (i in 1:(nr-1)) {
        for (j in 1:(nc-1)) {
            P[i, j] <- mnormt::sadmvn(lower = c(row.cuts[i], col.cuts[j]), 
                upper = c(row.cuts[i + 1], col.cuts[j + 1]), mean=rep(0,2),
                varcov = R)   #should we specify the algorithm to TVPACK or Miwa
        }}
    P[1,nc] <- pnorm(rc[1]) - sum(P[1,1:(nc-1)] )
    P[nr,1] <- pnorm(cc[1]) - sum(P[1:(nr-1),1] )
    if(nr >2) {for (i in (2:(nr-1))) {P[i,nc] <- pnorm(rc[i]) -pnorm(rc[i-1])- sum(P[i,1:(nc-1)] ) }}
    if(nc >2) {for (j in (2:(nc-1))) {P[nr,j] <- pnorm(cc[j]) - pnorm(cc[j-1])-sum(P[1:(nr-1),j] ) }}
    if(nc > 1)  P[nr,nc] <- 1- pnorm(rc[nr-1]) - sum(P[nr,1:(nc-1)]) 
    P   #the estimated n x n predicted by rho, rc, cc
}


#  polyF <- function(rho,rc,cc,tab) { 
#       P <- polyBinBvn(rho, rc, cc) 
#        -sum(tab * log(P)) }  #the  criterion to be minimized
#        
#revised 16/6/18 to cover the problem of 0 values in cells
polyF <- function(rho,rc,cc,tab) {  #doesn't blow up in the case of 0 cell entries  added 16/6/18
      P <- polyBinBvn(rho, rc, cc) 
     P[P <=0] <- NA   #added 18/2/9
     lP <- log(P)
      lP[lP == -Inf] <- NA
       lP[lP == Inf] <- NA
       -sum(tab * lP,na.rm=TRUE)  }  #the  criterion to be minimized
       
       
"wtd.table" <- function(x,y,weight) {
tab <- tapply(weight,list(x,y),sum,na.rm=TRUE,simplify=TRUE) #taken from questionr:wtd.table
tab[is.na(tab)] <- 0
return(tab)
}  

    
 
#modified 10/8/14 to create missing values when there are no cell entries
#modified 3/6/14 to create missing values when the data are hopeless
#modified 06/2/18 for the case of empty cells
 "polyc" <-     #uses the tableFast function instead of tableF
function(x,y=NULL,taux,tauy,global=TRUE,weight=NULL,correct=correct,gminx,gmaxx,gminy,gmaxy) {
  if(is.null(weight )) {tab <- tableFast(x,y,gminx,gmaxx,gminy,gmaxy)
   }   else {tab <- wtd.table(x,y,weight)}  #need to specify minx and maxx somehow
   fixed <- 0 
  tot <- sum(tab)
  if(tot ==0) {result <- list(rho=NA,objective=NA,fixed=1)
               return(result)} #we have no data for this cell   05/02/18
  tab <- tab/tot
    
 if(global) { rho <- optimize(polyF,interval=c(-1,1),rc=taux, cc=tauy,tab)#this uses the global taux and tauy
       } else { #use item row and column information for this pair, rather than global values
      #this seems to match the polycor function
      #the next five lines are adapted directly from John Fox's polycor function
    if(!is.na(sum(tab)) )  {   #this checks for completely missing data
   	 	zerorows <- apply(tab, 1, function(x) all(x == 0))
		zerocols <- apply(tab, 2, function(x) all(x == 0))
		zr <- sum(zerorows)
		zc <- sum(zerocols)
		tab <- tab[!zerorows, ,drop=FALSE]  
		tab <- tab[, !zerocols, drop=FALSE] 
		
		if(correct > 0) {if(any(tab[]==0)) {fixed <- 1
                 tab[tab==0] <- correct/tot }} #moved from below 16.6.22  which introduced a bug, moved to here 7/24/20
                 
                 
    	csum <- colSums(tab)
    	rsum <- rowSums(tab)
    	#if(correct > 0) tab[tab==0] <- correct/tot
   		if(min(dim(tab)) < 2) {rho <- list(objective = NA) } else {
    	 	 cc <-  qnorm(cumsum(csum)[-length(csum)])
     		 rc <-  qnorm(cumsum(rsum)[-length(rsum)])
    	 	rho <- optimize(polyF,interval=c(-1,1),rc=rc, cc=cc,tab)
    	 	
        	}
        } else { rho <- list(objective = NA, rho= NA)}}
     if(is.na(rho$objective)) {result <- list(rho=NA,objective=NA,fixed=fixed) } else {
            result <- list(rho=rho$minimum,objective=rho$objective,fixed=fixed)}
   
  return(result)
  }
##########

#We have dropped option to use John Fox's polycor package, so we don't need the options
#function(x,smooth=TRUE,global=TRUE,polycor=FALSE,ML = FALSE, std.err = FALSE,weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,delete=TRUE) {
#12/25/19 added the ability to have a y set of  variables as well
"polychoric" <- 
function(x,y=NULL,smooth=TRUE,global=TRUE,polycor=FALSE,ML = FALSE, std.err = FALSE,weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,delete=TRUE,max.cat=8)  {
#function(x,smooth=TRUE,global=TRUE,polycor=FALSE,weight=NULL,correct=.5,progress=TRUE,na.rm=TRUE,delete=TRUE) {
#if(!require(parallel)) {message("polychoric requires the parallel package.")}
#declare these next two functions to be local inside of polychoric
#The polycor paramater was dropped because it was not being used.  But, several programs are using it.
if(polycor) message("The polycor option has been removed from the polychoric function in the psych package.  Please fix the call.")
if(ML) message("The ML option has been removed from the polychoric function in the psych package.  Please fix the call.")
if(std.err) message("The std.error option has been removed from the polychoric function in the psych package.  Please fix the call.")



myfun <- function(x,y,i,j,gminx,gmaxx,gminy,gmaxy) {polyc(x[,i],x[,j],tau[,i],tau[,j],global=global,weight=weight,correct=correct,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy) }
myfuny <- function(x,y,i,j,gminx,gmaxx,gminy,gmaxy,tauy) {polyc(x[,i],y[,j],tau[,i],tauy[,j],global=global,weight=weight,correct=correct,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy) }

matpLower <- function(x,nvar,gminx,gmaxx,gminy,gmaxy) {
k <- 1
il <- vector()
jl <- vector()
for(i in 2:nvar) {for (j in 1:(i-1)) {
   il[k] <- i
   jl [k] <- j
   k<- k+1}
   }
poly <- mcmapply(function(i,j) myfun(x,y=NULL,i,j,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy) , il,jl) 
#poly <- mapply(function(i,j) myfun(x,y=NULL,i,j,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy) , il,jl) 
#debugging, we turn off the mcmapply function and do it by hand
# browser()
# ppl <- list()
# for (i in 2:nvar) {for (j in 1:(i-1)) {ppl[[i+j]] <- myfun(x,i,j,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy) } }

#now make it a matrix
mat <- diag(nvar)
if(length(dim(poly)) == 2) {
mat[upper.tri(mat)] <- as.numeric(poly[1,]) #first row of poly is correlation, 2nd the fit
 mat <- t(mat) + mat
 fixed <- as.numeric(poly[3,])
diag(mat) <- 1 
fixed <- sum(fixed) 
if((fixed > 0) && ( correct > 0)) { warning(fixed ," cells were adjusted for 0 values using the correction for continuity. Examine your data carefully.")}

return(mat)} else {
warning("Something is wrong in polycor ")
return(poly)
#never actually gets here
stop("we need to quit because something was seriously wrong.  Please look at the results")} 
}

matpxy <- function(x,y,nvar,nvar.y,gminx,gmaxx,gminy,gmaxy,tauy) {


if(!is.matrix(tauy)) tauy <- matrix(tauy,ncol=nvar.y)
#make the lists that are passed to myfuny
k <- 1
il <- vector()
jl <- vector()
for(i in 1:nvar) {for (j in 1:(nvar.y)) {
   il[k] <- i
   jl [k] <- j
   k<- k+1}
   }
 
   
#Now, we apply the function myfuny and return the polychoric of the x by y
#poly <- mapply(function(i,j) myfuny(x,y,i,j,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy,tauy) , il,jl) 
poly <- mcmapply(function(i,j) myfuny(x,y,i,j,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy,tauy) , il,jl) 
 
#now make it a matrix
mat <-  matrix(NA,ncol=nvar,nrow=nvar.y)
if(length(dim(poly)) == 2) {
mat<-  as.numeric(poly[1,]) #first row of poly is correlation, 2nd the fit

 fixed <- as.numeric(poly[3,])

fixed <- sum(fixed) 
if((fixed > 0) && ( correct > 0)) { warning(fixed ," cells were adjusted for 0 values using the correction for continuity. Examine your data carefully.")}

return(mat)} else {
warning("Something is wrong in polycor ")
return(poly)
#never actually gets here
stop("we need to quit because something was seriously wrong.  Please look at the results")} 
}

#the main function starts here

 if(!is.null(weight)) {if(length(weight) !=nrow(x)) {stop("length of the weight vector must match the number of cases")}}
 cl <- match.call() 
 
nvar <- dim(x)[2]
nsub <- dim(x)[1]
if((prod(dim(x)) == 4) | is.table(x))  {result <- polytab(x,correct=correct)
                     print("You seem to have a table, I will return just one correlation.") } else {  #the main function
x <- as.matrix(x)
if(!is.numeric(x)) {x <- matrix(as.numeric(x),ncol=nvar)
    message("Converted non-numeric input to numeric")}
    

nvalues <- max(x,na.rm=TRUE) - min(x,na.rm=TRUE) + 1
if(nvalues > max.cat) stop("You have more than", max.cat," categories for your items, polychoric is probably not needed")
 #first delete any bad cases
  item.var <- apply(x,2,sd,na.rm=na.rm)
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0)  & delete) {
            for (baddy in 1:length(bad)) {message( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
            nvar <- nvar - length(bad)
             }
 

 xmin <- apply(x,2,function(x) min(x,na.rm=TRUE))
#if(global)  { xmin <- min(xmin)} 
 xmin <- min(xmin)
x <- t(t(x) - xmin +1)  #all numbers now go from 1 to nvalues
 
  gminx <- gminy <- 1  #allow for different minima if minmax is null
  xmax <- apply(x,2,function(x)  max(x,na.rm=TRUE))
#if(global) xmax <- max(xmax) 
if (min(xmax) != max(xmax)) {global <- FALSE
                      warning("The items do not have an equal number of response alternatives, global set to FALSE.")}
    
 xmax <- max(xmax)  #don't test for globality xmax
 gmaxx <- gmaxy <- xmax #check for different maxima

xfreq <- apply(x,2,tabulate,nbins=nvalues)
n.obs <- colSums(xfreq)
xfreq <- t(t(xfreq)/n.obs)
tau <- qnorm(apply(xfreq,2,cumsum)[1:(nvalues-1),])  #these are the normal values of the cuts   #moved ) to after ]
if(!is.matrix(tau)) tau <- matrix(tau,ncol=nvar)
#rownames(tau) <- levels(as.factor(x))[1:(nvalues-1)]  #doesn't work if one response is missing
rownames(tau) <- 1:(nvalues -1) 
colnames(tau) <- colnames(x)

mat <- matrix(0,nvar,nvar)
colnames(mat) <- rownames(mat) <- colnames(x)
#x <- x - min(x,na.rm=TRUE) +1  #this is essential to get the table function to order the data correctly -- but we have already done it
if(is.null(y)) {
mat <- matpLower(x,nvar,gminx,gmaxx,gminy,gmaxy)  #the local copy has the extra paremeters   #do the multicore version


  if(any(is.na(mat))) {message("some correlations are missing, smoothing turned off")
                        smooth <- FALSE}

 if(smooth) {mat <- cor.smooth(mat) }

 colnames(mat) <- rownames(mat) <- colnames(x)
 tau <- t(tau)
 tauy<- NULL
 } else { #process the x * y data
   ymin <- apply(y,2,function(x) min(x,na.rm=TRUE))
   ymin <- min(ymin,na.rm=TRUE)
   
   nvar.y <- NCOL(y)
   y <- t((t(y) - ymin +1))  #all numbers go from 1 to ymax +1
   ymax <- apply(y,2,function(x) max(x,na.rm=TRUE))
   ymax <- max(ymax,na.rm=TRUE)
   gminy <- 1
   gmaxy <- ymax 
   nvaluesy <- ymax - ymin +1
   yfreq <- apply(y,2,tabulate,nbins=nvaluesy)
 	n.obs.y <- colSums(yfreq)
	yfreq <- t(t(yfreq)/n.obs.y)
	tauy <- qnorm(apply(yfreq,2,cumsum))[1:(nvalues-1),]
	if(!is.matrix(tauy)) tauy <- matrix(tauy,ncol=nvar.y) 
	rownames(tauy) <- 1:(nvalues-1)
	colnames(tauy) <- colnames(y)
	
   mat <- matpxy(x,y,nvar,nvar.y,gminx,gmaxx,gminy,gmaxy,tauy)
   mat <- matrix(mat,ncol=nvar,nrow=nvar.y)
   colnames(mat )<- colnames(x)
   rownames(mat) <- colnames(y)
   tauy <- t(tauy)
   tau <- t(tau)
   }
  result <- list(rho = mat,tau = tau,tauy = tauy,n.obs=nsub,Call=cl) 
  
 class(result) <- c("psych","poly")
 }
  return(result) 
  }
  
  #####



#use polychor from John Fox to do the same
#matches polychor output perfectly if correct=FALSE
#polytab is for just one pair of variables 
"polytab" <- 
function(tab,correct=TRUE) {
  
  tot <- sum(tab)
  tab <- tab/tot
  if(correct > 0) tab[tab==0] <- correct/tot 
 #use item row and column information for this pair, rather than global values  
    csum <- colSums(tab)
    rsum <- rowSums(tab)
     cc <-  qnorm(cumsum(csum[-length(csum)]))
     rc <-  qnorm(cumsum(rsum[-length(rsum)]))
     rho <- optimize(polyF,interval=c(-1,1),rc=rc, cc=cc,tab)

  result <- list(rho=rho$minimum,objective=rho$objective,tau.row=rc,tau.col =cc)
  return(result)
  }

##########################################################################################################

#9/6/14  to facilitate mixed cor  we find polytomous by dichotomous correlations
#4/08/17  fixed to not do table(p) or table(d)  
#has a problem if we are correcting 0 values 
"polydi" <- function(p,d,taup,taud,global=TRUE,ML = FALSE, std.err = FALSE,weight=NULL,progress=TRUE,na.rm=TRUE,delete=TRUE,correct=.5) {

#if(!require(parallel)) {message("polychoric requires the parallel package.")}
#declare these next two functions to be local inside of polychoric

myfun <- function(x,i,j,correct,taup,taud,gminx,gmaxx,gminy,gmaxy,np) {polyc(x[,i],x[,j],taup[,i],taud[1,(j-np)],global=global,weight=weight,correct=correct,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy) } #global changed to true 16/6/19  and set back again to global=global on 09/07/17

matpLower <- function(x,np,nd,taup,taud,gminx,gmaxx,gminy,gmaxy) {
k <- 1
il <- vector()
jl <- vector()
for(i in 1:np) {for (j in 1:nd) {
   il[k] <- i
   jl [k] <- j
   k <- k+1}
   }
   
poly <- mcmapply(function(i,j) myfun(x,i,j,correct=correct,taup=taup,taud=taud,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy,np=np) , il,jl+np)  #the multicore version
#poly <- mapply(function(i,j) myfun(x,i,j,correct=correct,taup=taup,taud=taud,gminx=gminx,gmaxx=gmaxx,gminy=gminy,gmaxy=gmaxy,np=np) , il,jl +np)   #the normal version for debugging
#now make it a matrix
mat <- matrix(np,nd)
mat <- as.numeric(poly[1,]) #first row of poly is correlation, 2nd the fit
return(mat)
}


#if(!require(mnormt) ) {stop("I am sorry, you must have mnormt installed to use polychoric")}

 if(!is.null(weight)) {if(length(weight) !=nrow(x)) {stop("length of the weight vector must match the number of cases")}}
 cl <- match.call() 
np <- dim(p)[2]
nd <- dim(d)[2]
if(is.null(np)) np <- 1
if(is.null(nd)) nd <- 1
nsub <- dim(p)[1]
p <- as.matrix(p)
d <- as.matrix(d)
#pt <- table(p)    #why do we do this?
#nvalues <- length(xt)  #find the number of response alternatives 
nvalues <- max(p,na.rm=TRUE) - min(p,na.rm=TRUE) + 1
#dt <- table(d)
dmin <- apply(d,2,function(x) min(x,na.rm=TRUE))
dmax <- apply(d,2,function(x) max(x,na.rm=TRUE))
dvalues <- max(dmax-dmin) 


if(dvalues !=1) stop("You did not supply a dichotomous variable")
if(nvalues > 8) stop("You have more than 8 categories for your items, polychoric is probably not needed")
 #first delete any bad cases
  item.var <- apply(p,2,sd,na.rm=na.rm)
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0)  & delete) {
            for (baddy in 1:length(bad)) {message( "Item = ",colnames(p)[bad][baddy], " had no variance and was deleted")}
            p <- p[,-bad] 
            np <- np - length(bad)
             }
pmin <- apply(p,2,function(x) min(x,na.rm=TRUE))  #allow for different minima
#gminx <- min(pmin)
minx <- min(pmin)
p <- t(t(p) - pmin +1)  #all numbers now go from 1 to nvalues 
#p <- t(t(p) - gminx +1)  #all numbers now go from 1 to nvalues    but we should use global minimima

#gminy <- min(dmin)
miny <- min(dmin)
#d <-  t(t(d) - gminy +1)
d <-  t(t(d) - dmin +1)   #this allows a separate minimum for each d variable
gminx <- gminy <- 1    #set the global minima to 1

pmax <- apply(p,2,function(x)  max(x,na.rm=TRUE)) #check for different maxima
gmaxx <- max(pmax)

if (min(pmax) != max(pmax)) {global <- FALSE
          warning("The items do not have an equal number of response alternatives, I am setting global to FALSE")}
gmaxy <- max(apply(d,2,function(x) max(x,na.rm=TRUE)))                      
#xfreq <- apply(x- xmin + 1,2,tabulate,nbins=nvalues)
pfreq <- apply(p,2,tabulate,nbins=nvalues)
n.obs <- colSums(pfreq)
pfreq <- t(t(pfreq)/n.obs)
taup <- as.matrix(qnorm(apply(pfreq,2,cumsum))[1:(nvalues-1),],ncol=ncol(pfreq))  #these are the normal values of the cuts
#if(!is.matrix(tau)) tau <- matrix(tau,ncol=nvar)

#rownames(taup) <- names(pt)[1:(nvalues-1)]
rownames(taup) <- paste(1:(nvalues-1))
colnames(taup) <- colnames(p)

dfreq <- apply(d,2,tabulate,nbins=2)
if(nd < 2) {n.obsd <- sum(dfreq) } else {n.obsd <- colSums(dfreq) } 
dfreq <- t(t(dfreq)/n.obsd)
taud <-  qnorm(apply(dfreq,2,cumsum)) 
 
mat <- matrix(0,np,nd)
rownames(mat) <- colnames(p)
colnames(mat)  <- colnames(d)
#x <- x - min(x,na.rm=TRUE) +1  #this is essential to get the table function to order the data correctly
x <- cbind(p,d)

mat <- matpLower(x,np,nd,taup,taud,gminx,gmaxx,gminy,gmaxy)  #the local copy has the extra paremeters   #do the multicore version

 
 mat <- matrix(mat,np,nd,byrow=TRUE)
 rownames(mat) <- colnames(p)
colnames(mat)  <- colnames(d)
 taud <- t(taud)
  result <- list(rho = mat,tau = taud,n.obs=nsub,Call=cl) 
 class(result) <- c("psych","polydi")
  return(result) 
  }
