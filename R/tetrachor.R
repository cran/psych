#adapted from John Fox's Polychor
#this does all the work
"tetrac" <- 
function(x,y=NULL,taux,tauy,correct=TRUE,global=TRUE) {
 binBvn <- function (rho,rc,cc)    #adapted from John Fox's polychor
{ row.cuts <- c(-Inf, rc, Inf)
    col.cuts <- c(-Inf, cc, Inf)
    P <- matrix(0, 2,2)
    R <- matrix(c(1, rho, rho, 1), 2, 2)
    for (i in 1:2) {
        for (j in 1:2) {
            P[i, j] <- pmvnorm(lower = c(row.cuts[i], col.cuts[j]), 
                upper = c(row.cuts[i + 1], col.cuts[j + 1]), 
                corr = R)
        }}
    P   #the estimated 2 x 2 predicted by rho, rc, cc
}
 f <- function(rho,cc,rc) { 
      P <- binBvn(rho, cc, rc) 
       -sum(tab * log(P)) }  #the ML criterion to be minimized
      
 if(is.null(y)) {tab <- x} else {tab <- table(x,y)}
 if((sum(tab) > 1) && (min(tab) == 0) && correct) {
    warning("A cell entry of 0 was replaced with .5.  Check your data!")
    tab[tab==0] <-.5  #correction for continuity

    }
  
  if(global) {cc <- taux
              rc <- tauy } else {
  tot <- sum(tab)
  tab <- tab/tot
  rc <- qnorm(colSums(tab))[1]
  cc <- qnorm(rowSums(tab))[1]
  } 
  rho <- optimize(f,interval=c(-1,1),rc=rc,cc=cc)
  result <- list(rho=rho$minimum,tau=c(cc,rc),objective=rho$objective)
  return(result)
  }
  
 #repeatedly do the analysis to form a matrix of output 
"tetra.mat" <- 
function(x,y=NULL,correct=TRUE,smooth=TRUE,global=TRUE) {nvar <- dim(x)[2]
x <- x -min(x,na.rm=TRUE) #in case the numbers are not 0,1
n.obs <- dim(x)[1]
if(is.null(y)) {
if(max(x,na.rm=TRUE) > 1) {stop("Tetrachoric correlations require dictomous data")}
tau <- -qnorm(colMeans(x,na.rm=TRUE))

mat <- matrix(0,nvar,nvar)
colnames(mat) <- rownames(mat) <- colnames(x)
names(tau) <- colnames(x)
for (i in 2:nvar) {
  for (j in 1:(i-1)) {
  tetra <-  tetrac(x[,i],x[,j],tau[i],tau[j],correct=correct,global=global)
   mat[i,j] <- mat[j,i] <- tetra$rho
   }
   }
   diag(mat) <- 1
  if(smooth) {mat <- cor.smooth(mat) }  #makes it positive semidefinite 
  result <- list(rho = mat,tau = tau,n.obs=n.obs) } else {
  
      # the case of having a y variable
      y <- y -min(y,na.rm=TRUE) #in case the numbers are not 0,1 
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
   
 #the public function
 "tetrachoric" <- 
 function(x,y=NULL,correct=TRUE,smooth=TRUE,global=TRUE) {
 if(!require(mvtnorm)) {stop("I am sorry, you must have mvtnorm installed to use tetrachoric")}
 cl <- match.call() 
 if (!is.matrix(x) && !is.data.frame(x)) {
  if (length(x) ==4) {x <- matrix(x,2,2) } else {
  if(length(x) ==3 ) {x <- pqr(x) } else {
  stop("Data must be either a 1 x 4 vector, a 2 x 2 matrix, or a data.frame/matrix of data")} 
  }}
  nvar <- dim(x)[2]
  n.obs <- dim(x)[1]
 if (n.obs == nvar) {result <- tetrac(x,correct=correct,global=FALSE)} else {
 result <- tetra.mat(x,y=y,correct=correct,smooth=smooth,global=global)}
 
 result$Call <- cl
 class(result) <- c("psych","tetra")
 return(result) 
 }
 
 "tetrachor" <- 
 function(x,correct=TRUE) {
 if(!require(mvtnorm)) {stop("I am sorry, you must have mvtnorm installed to use tetrachor")}
 cl <- match.call() 
 if (!is.matrix(x) && !is.data.frame(x)) {
  if (length(x) ==4) {x <- matrix(x,2,2) } else {
  if(length(x) ==3 ) {x <- pqr(x) } else {
  stop("Data must be either a 1 x 4 vector, a 2 x 2 matrix, or a data.frame/matrix of data")} 
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
function(x,y) {
cc <- complete.cases(x,y)
x <- x[cc]
y <- y[cc]
yf <- as.factor(y)
lev <- levels(yf)
if(length(lev)!=2) {stop("y is not a dichotomous variable")}

ty <- table(y)
 tot <- sum(ty)
 tab <- ty/tot
 zp <- dnorm(qnorm(tab[2]))
 hi <- mean(x[y==lev[2]],na.rm=TRUE)
 lo <- mean(x[y==lev[1]],na.ram=TRUE)
# r <- (hi - lo)*sqrt(prod(tab))/(sd(x,na.rm=TRUE))  #point biserial
 r <- (hi - lo)*(prod(tab))/(zp * sd(x,na.rm=TRUE))
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
mat <- matrix(NA,nrow=ny,ncol=nx)
colnames(mat) <- colnames(x)
rownames(mat) <- colnames(y)
for(i in 1:ny) {
   for (j in 1:nx) {
    mat[i,j] <- biserialc(x[,j],y[,i])
    }}
   return(mat)
}


"polyserial" <-
function(x,y) {
   min.item <- min(y, na.rm = TRUE)
    max.item <- max(y, na.rm = TRUE)
     n.var <- dim(y)[2]
        n.cases <- dim(y)[1]
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



#December 22,2010
#revised July 15, 2011 to work for the various special cases
#meant to combine continuous, polytomous and dichotomous correlations
#revised October 12, 2011 to get around the sd of vectors problem
"mixed.cor" <-
function(x=NULL,p=NULL,d=NULL,smooth=TRUE) {
if(!is.null(x)) {nx <- dim(x)[2]} else {nx <- 0}
if(!is.null(p)) {np <- dim(p)[2]} else {np <- 0}
if(!is.null(d))  {nd <- dim(d)[2]} else {nd <- 0}
if(is.null(nx)) nx <- 1
if(is.null(np)) np <- 1
if(is.null(nd)) nd <- 1
npd  <- nx +np + nd
if(nx > 0) {rx <- cor(x,use="pairwise")} else {rx <- NULL
   rho <- NULL}
if(np > 1) {rp <- polychoric(p,smooth=smooth)}    else {if (np == 1) {
	rho <- 1
	names(rho) <- colnames(p)
	rp <- list(rho=rho,tau=NULL)}  else {rp <- NULL}}
if(nd > 1) {rd <- tetrachoric(d,smooth=smooth)}   else {if (nd == 1) {rd <- list(rho=1,tau=NULL)}  else {rd <- NULL}}

if(nx > 0) {if(np > 0) {rxp <- polyserial(x,p)   #the normal case is for all three to exist
		tmixed <- cbind(rx,t(rxp))
		lmixed <- cbind(rxp,rp$rho)
		rho <- rbind(tmixed,lmixed)} else {rho <- rx}  #we now have x and p
		
		if(nd > 0) { rxd <- biserial(x,d)
		    if(np > 0) {rpd <- biserial(p,d) 
			            topright <- t(cbind(rxd,rpd)) 
			            } else {
			        topright <- t(rxd)}
			tmixed <- cbind(rho,topright) 
			lmixed <- cbind(t(topright),rd$rho)
			rho <- rbind(tmixed,lmixed) }

		} else {  #the case of nx =0
		   if( (np > 0) & (nd >0 )) {rpd <- biserial(p,d)
		       tmixed <- cbind(rp$rho,t(rpd))
		       lmixed <- cbind(rpd,rd$rho)
		       rho <- rbind(tmixed,lmixed)
		       }
		    }
colnames(rho) <- rownames(rho)
mixed <- list(rho=rho,rx=rx,poly=rp,tetra=rd)
return(mixed)
}




"cor.smooth" <- function(x) {
eigens <- eigen(x)
if(min(eigens$values) < .Machine$double.eps)  {warning("Matrix was not positive definite, smoothing was done")
eigens$values[eigens$values  < .Machine$double.eps] <- 100 * .Machine$double.eps
nvar <- dim(x)[1]
tot <- sum(eigens$values)
eigens$values <- eigens$values * nvar/tot
cnames <- colnames(x)
rnames <- rownames(x)
x <- eigens$vectors %*% diag(eigens$values) %*% t(eigens$vectors)
x <- cov2cor(x)
colnames(x) <- cnames
rownames(x) <- rnames}
return(x)}
