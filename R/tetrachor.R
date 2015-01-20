#adapted from John Fox's Polychor
#this does all the work

 #the following two functions are called repeatedly by tetrac and are put here to speed up the process
 
# "tetraBinBvn.old" <-
#  function (rho,rc,cc)    #adapted from John Fox's polychor
# { row.cuts <- c(-Inf, rc, Inf)
#     col.cuts <- c(-Inf, cc, Inf)
#     P <- matrix(0, 2,2)
#     R <- matrix(c(1, rho, rho, 1), 2, 2)
#     for (i in 1:2) {
#         for (j in 1:2) {
#             P[i, j] <- pmvnorm(lower = c(row.cuts[i], col.cuts[j]), 
#                 upper = c(row.cuts[i + 1], col.cuts[j + 1]), 
#                 corr = R)
#         }}
#     P   #the estimated 2 x 2 predicted by rho, rc, cc
# }

#modified 5/8/14 to be consistent with call from tetraF
#no change in functionality, just more esthetic

#changed 10/16/14  to use sadmvn instead of mvtnorm
"tetraBinBvn" <-
 function (rho,rc,cc)    #adapted from John Fox's polychor
{ row.cuts <- c(-Inf, rc, Inf)
    col.cuts <- c(-Inf, cc, Inf)
    P <- matrix(0, 2,2)
    R <- matrix(c(1, rho, rho, 1), 2, 2)
    P[1,1] <- sadmvn(lower = c(row.cuts[1], col.cuts[1]), 
                upper = c(row.cuts[2], col.cuts[2]), mean=c(0,0),
                varcov = R)
    P[2,1] <- pnorm(rc) - P[1,1]
    P[1,2] <- pnorm(cc) - P[1,1]
    P[2,2] <-  1- pnorm(rc) - P[1,2]

    P   #the estimated 2 x 2 predicted by rho, rc, cc
}


"tetraF" <-
 function(rho,cc,rc,tab) { 
      P <- tetraBinBvn(rho, rc, cc) 
       -sum(tab * log(P)) }  #the ML criterion to be minimized

"tetrac" <- 
function(x,y=NULL,taux,tauy,i,j,correct=.5,global=TRUE,weight=NULL) {

      
 if(is.null(y)) {tab <- x} else {
 if(is.null(weight)) {tab <- tableF(x,y) }  else {tab <- wtd.table(x,y,weight)}  #switched to tableF for speed
 }
 
 #changed 9/8/14 
#if((length(tab) < 4) | (is.na(sum(tab)) | ((tab[1,1] + tab[1,2]) < 1) |  ((tab[2,1] + tab[2,2]) < 1)  | ((tab[1,1] + tab[2,1]) < 1) |  ((tab[2,1] + tab[2,2]) < 1))) {warning("For i = ", i," j = ",j, "  No variance for either i or j   rho set to NA")
 if((length(tab) < 4) | (is.na(sum(tab)) )) {warning("For i = ", i," j = ",j, "  No variance for either i or j   rho set to NA")          
      result <- list(rho=NA,tau=c(NA,NA),objective=NA)
          }  else   {
              
 if((sum(tab) > 1) && (min(tab) == 0) && (correct > 0)) {
    warning("For i = ", i," j = ",j, "  A cell entry of 0 was replaced with .5.  Check your data!")
    tab[tab==0] <- correct  #correction for continuity

    }
  ###### put in the weights here
  if(global) {cc <- taux
              rc <- tauy } else {
  tot <- sum(tab)
  tab <- tab/tot
  rc <- qnorm(colSums(tab))[1]
  cc <- qnorm(rowSums(tab))[1]
  } 
  rho <- optimize(tetraF,interval=c(-1,1),rc=rc,cc=cc,tab=tab)
  result <- list(rho=rho$minimum,tau=c(cc,rc),objective=rho$objective) }
  return(result)
  }
  
 "wtd.table" <- function(x,y,weight) {
tab <- tapply(weight,list(x,y),sum,na.rm=TRUE,simplify=TRUE) #taken from questionr:wtd.table
tab[is.na(tab)] <- 0
return(tab)
}    
 
  #repeatedly do the analysis to form a matrix of output 
 #added the pmin instead of min on Sept 10, 2013
"tetra.mat" <- 
function(x,y=NULL,correct=.5,smooth=TRUE,global=TRUE,weight=NULL) {

#the functions to do parallelism
myfun <- function(x,i,j) {if(t(!is.na(x[,i]))%*% (!is.na(x[,j]))  > 2 ) {
  tetra <-  tetrac(x[,i],x[,j],tau[i],tau[j],i,j,correct=correct,global=global,weight=weight)} else  {
  tetra <- NA}}


matpLower <- function(x,nvar) {
k <- 1
il <- vector()
jl <- vector()
for(i in 2:nvar) {for (j in 1:(i-1)) {
   il[k] <- i
   jl [k] <- j
   k<- k+1}
   }
tet <- mcmapply(function(i,j) myfun(x,i,j) , il,jl)
#now make it a matrix
mat <- diag(nvar)
mat[upper.tri(mat)] <- as.numeric(tet[1,]) #first row of poly is correlation, 2nd the fit
mat <- t(mat) + mat
diag(mat) <- 1
return(mat)
}


nvar <- dim(x)[2]
mx <- apply(x,2,function(x) min(x,na.rm=TRUE))
x <- t(t(x) - mx)
#x <- x -min(x,na.rm=TRUE) #in case the numbers are not 0,1  -- using pmin allows different minima for different variables
n.obs <- dim(x)[1]
if(is.null(y)) {
if(max(x,na.rm=TRUE) > 1) {stop("Tetrachoric correlations require dictomous data")}
if(is.null(weight)) {tau <- -qnorm(colMeans(x,na.rm=TRUE))} else
        {tau <- -qnorm(apply(x,2,function(y) weighted.mean(y,weight,na.rm=TRUE)))}   #weighted tau

mat <- matrix(0,nvar,nvar)
colnames(mat) <- rownames(mat) <- colnames(x)
names(tau) <- colnames(x)
#cat("\nFinding the tetrachoric correlations\n")
#for (i in 2:nvar) {
#progressBar(i^2/2,nvar^2/2,"Tetrachoric")
 # for (j in 1:(i-1)) {
 # if(t(!is.na(x[,i]))%*% (!is.na(x[,j]))  > 2 ) {
 # tetra <-  tetrac(x[,i],x[,j],tau[i],tau[j],i,j,correct=correct,global=global,weight=weight)
 #  mat[i,j] <- mat[j,i] <- tetra$rho} else {mat[i,j] <- mat[j,i] <- NA}
 #  }
  # }
  # diag(mat) <- 1
 mat <- matpLower(x,nvar)
 
   if(any(is.na(mat))) {warning("some correlations are missing, smoothing turned off")
                        smooth <- FALSE}
                     
  if(smooth) {mat <- cor.smooth(mat) }  #makes it positive semidefinite 
  colnames(mat) <- rownames(mat) <- colnames(x)  
  result <- list(rho = mat,tau = tau,n.obs=n.obs) } else {
  
      # the case of having a y variable
      my <- apply(x,2,function(x) min(x,na.rm=TRUE))
        y <- t(t(y) - my)
      #y <- y -min(y,na.rm=TRUE) #in case the numbers are not 0,1 
      if(is.matrix(y)) {ny <- dim(y)[2]
           tauy <- -qnorm(colMeans(y,na.rm=TRUE))
           n.obs.y <- dim(y)[1]
             } else {
             	ny <- 1
           		n.obs.y <- length(y)}
           	tauy <- -qnorm(mean(y,na.rm=TRUE))
           y <- as.matrix(y) 
           
      if(dim(x)[1] != n.obs.y)  {stop("x and y must have the same number of observations")}
      taux <- -qnorm(colMeans(x,na.rm=TRUE))
      
      nx <- dim(x)[2]
     
      mat <- matrix(0,nx,ny)
      colnames(mat) <- colnames(y)
      rownames(mat) <- colnames(x)
      for (i in 1:nx) {
        for (j in 1:ny) {tetra <-  tetrac(x[,i],y[,j],taux[i],tauy[j],correct=correct)
        mat[i,j] <- tetra$rho }
         }
      colnames(mat) <- rownames(mat) <- colnames(x)    
    result <-   list(rho = mat,tau = taux,tauy= tauy,n.obs=n.obs)
     }
 flush(stdout()) 
  cat("\n" )  #put in to clear the progress bar
  return(result) 
  }
  
 #convert comorbidity type numbers to a table
 pqr <- function(q1,q2=NULL,p=NULL) {
    if(length(q1) > 1) {
       q2 <- q1[2]
       p <- q1[3]
       q1 <- q1[1]}
   tab <- matrix(0,2,2)
   tab[1,1] <- p
   tab[2,1] <- q1-p
   tab[1,2] <- q2-p
   tab[2,2] <- 1-q1 - tab[1,2]
   return(tab)}

  
 #repeatedly do the analysis to form a matrix of output 
 #added the pmin instead of min on Sept 10, 2013
"tetra.mat.sc" <- 
function(x,y=NULL,correct=.5,smooth=TRUE,global=TRUE,weight=NULL) {



nvar <- dim(x)[2]
mx <- apply(x,2,function(x) min(x,na.rm=TRUE))
x <- t(t(x) - mx)
#x <- x -min(x,na.rm=TRUE) #in case the numbers are not 0,1  -- using pmin allows different minima for different variables
n.obs <- dim(x)[1]
if(is.null(y)) {
if(max(x,na.rm=TRUE) > 1) {stop("Tetrachoric correlations require dictomous data")}
if(is.null(weight)) {tau <- -qnorm(colMeans(x,na.rm=TRUE))} else
        {tau <- -qnorm(apply(x,2,function(y) weighted.mean(y,weight,na.rm=TRUE)))}   #weighted tau

mat <- matrix(0,nvar,nvar)

names(tau) <- colnames(x)
#cat("\nFinding the tetrachoric correlations\n")
for (i in 2:nvar) {
progressBar(i^2/2,nvar^2/2,"Tetrachoric")
  for (j in 1:(i-1)) {
  if(t(!is.na(x[,i]))%*% (!is.na(x[,j]))  > 2 ) {
  tetra <-  tetrac(x[,i],x[,j],tau[i],tau[j],i,j,correct=correct,global=global,weight=weight)
   mat[i,j] <- mat[j,i] <- tetra$rho} else {mat[i,j] <- mat[j,i] <- NA}
   }
   }
   diag(mat) <- 1
   if(any(is.na(mat))) {warning("some correlations are missing, smoothing turned off")
                        smooth <- FALSE}
                     
  if(smooth) {mat <- cor.smooth(mat) }  #makes it positive semidefinite 
   colnames(mat) <- rownames(mat) <- colnames(x)  
   
  result <- list(rho = mat,tau = tau,n.obs=n.obs) } else {
  
      # the case of having a y variable
      my <- apply(x,2,function(x) min(x,na.rm=TRUE))
        y <- t(t(y) - my)
      #y <- y -min(y,na.rm=TRUE) #in case the numbers are not 0,1 
      if(is.matrix(y)) {ny <- dim(y)[2]
           tauy <- -qnorm(colMeans(y,na.rm=TRUE))
           n.obs.y <- dim(y)[1]
             } else {
             	ny <- 1
           		n.obs.y <- length(y)}
           	tauy <- -qnorm(mean(y,na.rm=TRUE))
           y <- as.matrix(y) 
           
      if(dim(x)[1] != n.obs.y)  {stop("x and y must have the same number of observations")}
      taux <- -qnorm(colMeans(x,na.rm=TRUE))
      
      nx <- dim(x)[2]
     
      mat <- matrix(0,nx,ny)
      colnames(mat) <- colnames(y)
      rownames(mat) <- colnames(x)
      for (i in 1:nx) {
        for (j in 1:ny) {tetra <-  tetrac(x[,i],y[,j],taux[i],tauy[j],correct=correct)
        mat[i,j] <- tetra$rho }
         }
         
     

    result <-   list(rho = mat,tau = taux,tauy= tauy,n.obs=n.obs)
     }
 flush(stdout()) 
  cat("\n" )  #put in to clear the progress bar
  return(result) 
  }
  
 #convert comorbidity type numbers to a table
 pqr <- function(q1,q2=NULL,p=NULL) {
    if(length(q1) > 1) {
       q2 <- q1[2]
       p <- q1[3]
       q1 <- q1[1]}
   tab <- matrix(0,2,2)
   tab[1,1] <- p
   tab[2,1] <- q1-p
   tab[1,2] <- q2-p
   tab[2,2] <- 1 - q1 - q2 + p
   return(tab)}
   
 #the public function
 "tetrachoric" <- 
 function(x,y=NULL,correct=.5,smooth=TRUE,global=TRUE,weight=NULL,na.rm=TRUE,delete=TRUE) {

 
# if(!require(mnormt)) {stop("I am sorry, you must have mnormt installed to use tetrachoric")}
 cl <- match.call() 
 if (!is.matrix(x) && !is.data.frame(x)) {
  if (length(x) ==4) {x <- matrix(x,2,2) } else {
  if(length(x) ==3 ) {x <- pqr(x) } else {
  stop("Data must be either a 1 x 4 vector, a 2 x 2 matrix, a comorbidity table, or a data.frame/matrix of data")} 
  }}
  nvar <- dim(x)[2]
  n.obs <- dim(x)[1]
  if(!is.null(weight)) {if (length(weight)!= n.obs) stop("The number of weights must match the number of observations") }
 if (n.obs == nvar) {result <- tetrac(x,correct=correct,i=1,j=1,global=FALSE)} else {
 #first delete any bad cases
  item.var <- apply(x,2,sd,na.rm=na.rm)
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0)  & delete) {
            for (baddy in 1:length(bad)) {warning( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
            nvar <- nvar - length(bad)
 
            }
            
 # parallel is now built into the system, so we don't need this.
 # if(!require(parallel)) {warning("need parallel installed to take advantage of multiple cores.  Using single core version instead")
 #   result <- tetra.mat.sc(x,y=y,correct=correct,smooth=smooth,global=global,weight=weight)} else {
 result <- tetra.mat(x,y=y,correct=correct,smooth=smooth,global=global,weight=weight)}
 
 result$Call <- cl
 class(result) <- c("psych","tetra")
 return(result) 
 }
 #modified 1/14/14 to include the tableF function to double the speed for large problems
 
 "tetrachor" <- 
 function(x,correct=.5) {
 #if(!require(mnormt)) {stop("I am sorry, you must have mnormt installed to use tetrachor")}
 cl <- match.call() 
 if (!is.matrix(x) && !is.data.frame(x)) {
  if (length(x) ==4) {x <- matrix(x,2,2) } else {
  if(length(x) ==3 ) {x <- pqr(x) } else {
   stop("Data must be either a 1 x 4 vector, a 2 x 2 matrix, a comorbidity table, or a data.frame/matrix of data")} 
   }}
  nvar <- dim(x)[2]
 if (dim(x)[1] == nvar) {result <- tetrac(x,correct=correct)} else {
 result <- tetra.mat(x,correct=correct)}

 result$Call <- cl
 class(result) <- c("psych","tetra")
 return(result) 
 }
 
  #does the work
"biserialc" <-
function(x,y,i,j) {
cc <- complete.cases(x,y)
x <- x[cc]
y <- y[cc]
yf <- as.factor(y)
lev <- levels(yf)
if(length(lev)!=2) {#stop("y is not a dichotomous variable")
                   warning("For x = ",i, " y = ", j, " y is not dichotomous")
                    r <- NA} else {
ty <- table(y)  
 tot <- sum(ty)
 tab <- ty/tot
 if(length(tab) < 2) {r <- NA
                    warning("For x = ",i, " y = ", j, " no variance for y   r set to NA")} else { #this treats the case of no variance in the dichotmous variable
 zp <- dnorm(qnorm(tab[2]))
 hi <- mean(x[y==lev[2]],na.rm=TRUE)
 lo <- mean(x[y==lev[1]],na.ram=TRUE)
# r <- (hi - lo)*sqrt(prod(tab))/(sd(x,na.rm=TRUE))  #point biserial
 r <- (hi - lo)*(prod(tab))/(zp * sd(x,na.rm=TRUE))
 if(!is.na(r) && abs(r)  >1 ) {  if (r > 1)  {r <- 1 } else {r <- -1}   #in case we are correlating a dichotomous variable with itself 
          warning("For x = ",i, " y = ", j, " x seems to be dichotomous, not continuous")
  }}}
 return(r)
}
 
"biserial" <- 
function(x,y) { 
x <- as.matrix(x)
y <- as.matrix(y)
nx <- dim(x)[2]
ny <- dim(y)[2]
if(is.null(nx)) nx <- 1
if(is.null(ny)) ny <- 1
mat <- matrix(NaN,nrow=ny,ncol=nx)
colnames(mat) <- colnames(x)
rownames(mat) <- colnames(y)
#cat("\n Finding the biserial correlations\n")
for(i in 1:ny) {
progressBar(i*(i-1)/2,ny^2/2,"Biserial")
   for (j in 1:nx) {
    mat[i,j] <- biserialc(x[,j],y[,i],j,i)
    }}
  flush(stdout())
  cat("\n" )  #put in to clear the progress bar
   return(mat)
}


"polyserial" <-
function(x,y) {
# y <- matrix(y)
   min.item <- min(y, na.rm = TRUE)
   max.item <- max(y, na.rm = TRUE)
 if(is.null(ncol(y)))  {n.var <- 1
                        n.cases <- length(y)
                } else {n.var <- ncol(y)
        n.cases <- dim(y)[1]}
        dummy <- matrix(rep(min.item:max.item, n.var), ncol = n.var)
        colnames(dummy) <- names(y)
        xdum <- rbind(y, dummy)
        frequency <- apply(xdum, 2, table)
        frequency <- t(frequency - 1)  
        responses <- rowSums(frequency)
        frequency <- frequency/responses
        frequency <- t(apply(frequency,1,cumsum))
        len <- dim(frequency)[2]
        tau <- dnorm(qnorm(frequency[,-len,drop=FALSE]))
        stau <- rowSums(tau)
  
    rxy <- cor(x,y,use="pairwise")
    sdy <- apply(y,2,sd,na.rm=TRUE)
    rps <- t(rxy) * sqrt((n.cases-1)/n.cases) * sdy/stau
    rps[rps > 1.0] <- 1.0
    rps[rps < -1.0] <- -1.0
 return(rps)
 }

#modified November 28, 2014 to be slightly more aggressive about smoothing
#this is more similar to cov2cor(nearPD$mat)
"cor.smooth" <- function(x,eig.tol=10^-12) {
eigens <- try(eigen(x),TRUE)
if(class(eigens)== as.character("try-error")) {warning('I am sorry, there is something seriously wrong with the correlation matrix,\ncor.smooth failed to  smooth it because some of the eigen values are NA.  \nAre you sure you specified the data correctly?')
                                                     } else {
                                                        
if(min(eigens$values) < .Machine$double.eps)  {warning("Matrix was not positive definite, smoothing was done")
#eigens$values[eigens$values  < .Machine$double.eps] <- 100 * .Machine$double.eps
eigens$values[eigens$values  < eig.tol] <- 100 * eig.tol
nvar <- dim(x)[1]
tot <- sum(eigens$values)
eigens$values <- eigens$values * nvar/tot
cnames <- colnames(x)
rnames <- rownames(x)
x <- eigens$vectors %*% diag(eigens$values) %*% t(eigens$vectors)
x <- cov2cor(x)
colnames(x) <- cnames
rownames(x) <- rnames}
}
return(x)}
#modified January 9, 2012 to add the try so we don't fail (just complain) if the data are bad.


#identify the most likely candidates for a bad item
"cor.smoother" <- function(x,cut=.01) {
nvar <- ncol(x)
result <- list()
if(nrow(x) != nvar) x <- cor(x,use="pairwise")
bad <- rep(NA,nvar)
good <- rep(TRUE,nvar)
names(good) <- names(bad) <- colnames(x)

for (i in 1:nvar) {
ev <- eigen(x[-i,-i])$values
if(any(ev < 0) ) {bad[i] <- TRUE
good[i] <- FALSE}
bad[i] <- sum((ev < 0),na.rm=TRUE)
}

if(sum(bad+0) > 0 )  {result$good <- bad[(bad > 0)]
result$bad <- good[good]
s <- cor.smooth(x)
possible <- arrayInd(which.max(abs(s-x)),dim(x),.dimnames=colnames(x))
result$likely <- colnames(x)[possible]

result$possible <- arrayInd(which(abs(s-x) > cut),dim(x),.dimnames=colnames(x))
result$possible <- sort(table(colnames(x)[result$possible]),decreasing=TRUE)
} else {result$bad <- c("all ok")}
class(result) <- c("psych","smoother")
return(result)
}

