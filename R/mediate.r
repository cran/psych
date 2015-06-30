#The basic logic is to find the regression of Y on X, of M on X, and of Y on M
#and then compare the direct (Y on X) also known as c  to c - ab 
#Modified May, 2015 to allow multiple Xs (and multiple Ys)
#Revised June, 2015 for prettier output
"mediate" <-
function(y,x,m,data, mod=NULL,n.obs=NULL,use="pairwise",n.iter=5000,alpha=.05,std=FALSE)  {
 
 #this first part just gets the variables right, depending upon the way we input them
   cl <- match.call()
   #convert names to locations 
   #first, see if they are in formula mode  
   if(class(y) == "formula") { yn <- all.vars(y)
            y <- yn[1]
            x <- yn[2]
            m <- yn[3:length(yn)]
            }  
    if(is.numeric(y )) y <- colnames(data)[y]
    if(is.numeric(x )) x <- colnames(data)[x]
    if(is.numeric(m )) m <- colnames(data)[m]
    if(!is.null(mod) && is.numeric(mod)) mod <- colnames(data)[mod]
    if(ncol(data) == nrow(data)) {raw <- FALSE 
                if(!is.null(mod)) {data <- data[c(y,x,m,mod),c(y,x,m,mod)] } else {data <- data[c(y,x,m),c(y,x,m)]} 
                 } else {
    if(!is.null(mod)) {data <- data[,c(y,x,m,mod)] } else {data <- data[,c(y,x,m)]}
 
  }
   if(!is.matrix(data)) data <- as.matrix(data)
  if((dim(data)[1]!=dim(data)[2]))  {n.obs=dim(data)[1]   #this does not take into account missing data
                    if(!is.null(mod)) data <- scale(data,scale=FALSE)  #0 center 
                    C <- cov(data,use=use)
                    raw <- TRUE
                     }  else {
                    raw <- FALSE
                    C <- data
 
   nvar <- ncol(C)
   if(is.null(n.obs)) {n.obs <- 1000
      message("The data matrix was a correlation matrix and the number of subjects was not specified. \n The data were simulated based upon the observed correlation matrix and  n.obs set to 1000")
   } else { message("The data matrix was simulated based upon the specified number of subjects and the observed correlation matrix.")}
       eX <- eigen(C)   #we use this in the simulation in the case of a correlation matrix
      data <- matrix(rnorm(nvar * n.obs),n.obs)
     data <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(data) )
     colnames(data) <- c(y,x,m)
  }
if(std) { C <- cov2cor(C)}  #use correlations rather than covariances
  
 
  if(!is.null(mod)) {ivXm <- matrix(data[,x] * data[,mod],ncol=1) #add in a moderating variable as a product term
         
       colnames(ivXm) <- paste(abbreviate(x),"X",abbreviate(mod),sep="")
      data <- cbind(data,ivXm)
        
              if(!(mod %in% m)) {m <- c(m,mod,colnames(ivXm))} else {m <- c(m,colnames(ivXm))}
              C <- cov(data,use=use)
              if(std) { C <- cov2cor(C)}
              }
  ######### 
  #########
  

  #now, we are ready to process the data
  #We do the basic regressions as matrix operations
     
        xy <- c(x,y)
        numx <- length(x)
     	numy <- length(y)
     	numm <- length(m)
        nxy <- numx + numy 
        m.matrix <- C[c(x,m),c(x,m),drop=FALSE] #this is the predictor +  moderator matrix
     	df <- n.obs - nxy - 1 

     	xy.matrix <- C[c(x,m),y,drop=FALSE]   #this is the matrix of correlations with the criterion (still just one)

   ##this is the zero order beta -- the total effect
   beta.x <-solve(C[x,x], (C[x,y]))
   
   #There are 3 broad cases that need to be handled somewhat differently in terms of the matrix operations
   # 1 IV, 1 or more mv 
   # multiple IV, 1 MV
   #multiple IV, multiple MV
   #this is the direct path from X to M  
    
    if(numx==1) { a <- solve(C[x,x],t(C[x,m])) } else {a <- solve(C[x,x],(C[x,m])) } 
     
     #R2 for the x predictor(s) by itself (themselves)
     R2.x <- (beta.x *C[x,y])/diag(C[y,y,drop=FALSE])
     	 uniq.x <- C[x,x]*(1- R2.x)

     	if(n.obs > 2) { se.bx <- sqrt(((1-R2.x) * (C[y,y])/ (diag(C[x,x,drop=FALSE]) * (n.obs-2) )  ))} else {se.bx <- NA}
     	 
     	tt <- beta.x/se.bx       #total direct t
     	probt <- 2*(1- pt(abs(tt),n.obs-2)) 
        beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY ~ aX
            
  beta <- beta.xm        #these are the individual predictors, x and then m and then mod
  b <- beta[-c(1:numx),drop=FALSE]   #these are the b paths 
  #  colnames(b) <- colnames(beta)
  #  rownames(b) <- rownames(beta)[-c(1:numx)]
    indirect <-  beta.x - beta[1:numx]  #this is c' = c - ab
    
#now here is where we need to consider the 3 options

    ab <- a * b    #these are the ab products for each a and b path   c' is c - sum of all of these
   
    if((numx>1)&& (numy > 1)) {for (i in 1:numx) {ab[i,] <- a[i,] * b}}
        R2 <- colSums(beta * xy.matrix)/(C[y,y]) 
    
    #the other way is to find the product of the paths
    #this is necessary for the interaction terms
   
    #predict m from x and mod
    if(!is.null(mod)) {
    Cm <- C[c(x,m[-1]),c(x,m[-1])]
    mM <- C[c(x,m[-1]),m[1],drop=FALSE]
	bM <- solve(Cm,mM)
	int.ind <- beta[2] * bM[nrow(bM)] 
	a <- bM[nrow(bM)]
	b.int = beta[nrow(beta)]} else {int.ind <- NA
	                        b.int <- NA}
     	 
   if(n.obs > 1)  {  uniq <- (1-smc(m.matrix))
     	                     df <- n.obs-1-nxy-numm  
     	                     se <- (sqrt((1-R2)/(df))*sqrt(1/uniq))
     	                     se <- t(t(se) * sqrt(diag(C)[y]))/sqrt(diag(m.matrix))
     	                         	               
     	                     tvalue <- beta/se
     	                    
     	                     k <- 1    
     	                     prob <- 2*(1- pt(abs(tvalue),df))
     	                     SE2 <- 4*R2*(1-R2)^2*(df^2)/((n.obs^2-1)*(n.obs+3))
     	                     SE =sqrt(SE2)
     	                     F <- R2*df/(k*(1-R2))
     	                     pF <- 1 - pf(F,k,df)
     	                     shrunkenR2 <- 1-(1-R2)*(n.obs-1)/df }
     ratit <- indirect/beta.x
     ratid <- indirect/beta.xm[1]	
     
     if(n.obs < 1) {message("Bootstrap is not meaningful unless raw data are provided") 
     
          mean.boot <- sd.boot <- ci.quant <- boot <-  se <- tvalue  <- prob <- NA } else {
            if(is.null(mod)) {boot <- boot.mediate(data,x,y,m,n.iter=n.iter,std=std,use=use)
               
               }  else {boot <- boot.moderate( data,x,y,m,mod,n.iter=n.iter,std=std,use=use)}
           
            mean.boot <- colMeans(boot)
            sd.boot <- apply(boot,2,sd)
            ci.quant <- apply(boot,2, function(x) quantile(x,c(alpha/2,1-alpha/2),na.rm=TRUE) )
              }   
      
  
     if(numx==1) {a <- t(a)
                  ab <- t(ab)
                  rownames(a) <- rownames(ab) <- m
                  }  
     b <- as.matrix(b)  #in order to get some names here
       colnames(b) <- colnames(beta)
    rownames(b) <- rownames(beta)[-c(1:numx)]
     if((numx > 1) && (numm > 1) ) {abc <- NA } else {       
    if(n.obs >  1)  {abc <- data.frame(a=a,b=b,ab=ab,mean.ab=mean.boot[-c(1:numx)],ci.ablower = ci.quant[1,-c(1:numx)],ci.abupper= ci.quant[2,-c(1:numx)])}   else {
   # {abc <- data.frame(a=a,b=b,ab=ab,mean.ab=mean.boot[-1],ci.ablower = ci.quant[1,-1],ci.abupper= ci.quant[2,-1])}  else {
   	 abc <- data.frame(a=a,b=b,ab=ab,NA,ci.ablower = NA,ci.abupper= NA)}
   	 }
      direct=beta.xm[1:numx,1,drop=FALSE]         
    result <- list(total=beta.x,se.bx=se.bx,tt=tt,probt=probt,direct=direct, indirect = indirect,se.beta=se,t=tvalue,prob=prob,R2.x=R2.x,R2=R2,
       ratit=ratit,ratid=ratid, mean.boot=mean.boot,sd.boot=sd.boot,ci.quant=ci.quant,boot=boot, mod=mod, int.ind=int.ind,a=a,b=b,ab=ab,b.int=b.int,abc=abc,mod=mod,names=colnames(data),C=C,  Call=cl)
   

    class(result) <- c("psych","mediate")
    if(is.null(mod)) { mediate.diagram(result) } else {moderate.diagram(result)}
    
    return(result)	
     		 }
     		 

boot.mediate <- function(data,x,y,m,n.iter=10,std=FALSE,use="pairwise") {
  n.obs <- nrow(data)
   numx <- length(x)
  numy <- length(y)
  numm <- length(m)
        nxy <- numx + numy 
        result <- matrix(NA,nrow=n.iter,ncol = (numx + numm*numx))
        
       
        
        
  for (iteration in 1:n.iter) {
  
  samp.data <- data[sample.int(n.obs,replace=TRUE),]
  C <- cov(samp.data,use=use)
  if(std) C <- cov2cor(C)
   xy <- c(x,y)
       
        m.matrix <- C[c(x,m),c(x,m)]
     	df <- n.obs - nxy - 1 

     	xy.matrix <- C[c(x,m),y,drop=FALSE]
     	
    beta.x <- solve(C[x,x],(C[x,y]) )  #this is the zero order beta -- the total effect
   if(numx ==1) { a <-  solve(C[x,x],t(C[x,m]) ) } else { a <-  solve(C[x,x],(C[x,m]) )} #the a paths    
   #R2.x <- beta.x *   C[x,y] /(C[y,y])  #R2 for the x predictor by itself 
  # R2.x <- (beta.x *C[x,y])/diag(C[y,y,drop=FALSE])
  #   	 uniq.x <- C[x,x]*(1- R2.x)
   #  	 se.bx <- sqrt(((1-R2.x) * (C[y,y])/ (diag(C[x,x,drop=FALSE]) * (n.obs-2) )  ))
   #  	tt <- beta.x/se.bx       #total direct t
   #  	probt <- 2*(1- pt(abs(tt),n.obs-2))

    beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
   beta <- as.matrix(beta.xm)     #these are the individual predictors, x and then m

 
  b <- beta[-c(1:numx),drop=FALSE]
  ab <- a*b    #each individual path  #this probably is only correct for the numx = 1 model
  if((numx>1)&& (numy > 1)) {for (i in 1:numx) {ab[i,] <- a[i,] * b}}
  indirect <-  beta.x - beta[1:numx]  #this is c' = c - ab

  result[iteration,] <- c(indirect,ab)
  
  }
 return(result) 
}

boot.moderate <- function(data,x,y,m,mod,n.iter=10,std=FALSE,use="pairwise") {
  
  n.obs <- nrow(data)
  numx <- 1
     	numy <- 1
     	numm <- length(m)
        nxy <- numx + numy 
   result <- matrix(NA,nrow=n.iter,ncol = (1+numm))
  for (iteration in 1:n.iter) {
  samp.data <- data[sample.int(n.obs,replace=TRUE),]
  samp.data <- scale(samp.data,scale=FALSE)
  ivXm <- samp.data[,x] * samp.data[,mod]
              samp.data <- cbind(samp.data,ivXm)
  C <- cov(samp.data,use=use)
  if(std) C <- cov2cor(C)
   xy <- c(x,y)
        numx <- 1
     	numy <- 1
     	numm <- length(m)
        nxy <- numx + numy 
        m.matrix <- C[c(x,m),c(x,m)]
     	df <- n.obs - nxy - 1 

     	xy.matrix <- C[c(x,m),y,drop=FALSE]
     	
    beta.x <- C[x,y]/C[x,x]   #this is the zero order beta -- the total effect
    
    beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
    beta <- as.matrix(beta.xm)   
    Cm <- C[c(x,m[-1]),c(x,m[-1])]
    mM <- C[c(x,m[-1]),m[1],drop=FALSE]
	bM <- solve(Cm,mM)
	int.ind <- beta[nrow(beta)] * bM[nrow(bM)]
         
   R2.x <- beta.x *   C[x,y] /(C[y,y])  #R2 for the x predictor by itself 
     
     	 uniq.x <- C[x,x]*(1- R2.x)
     	 se.bx <- sqrt(((1-R2.x) * (C[y,y])/ (C[x,x] * (n.obs-2) )  ))
     	tt <- beta.x/se.bx       #total direct t
     	probt <- 2*(1- pt(abs(tt),n.obs-2))

    beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
   beta <- as.matrix(beta.xm)     #these are the individual predictors, x and then m
      
  indirect <-  beta.x - beta[1]
  a <- bM #I think
 b <- beta[-1]
 ab <- a * b
 
  result[iteration,] <- c(int.ind,ab)
  }
 return(result) 
}


"print.psych.mediate" <- function(x,digits=2) {
 value <- class(x)[2] 

 switch(value, 
mediate  = { cat("\nMediation analysis \n")
#copy from here 
 cat("Call: ")
    print(x$Call)
    dv <- x$names[1]
    iv <- rownames(x$direct)
    mv <- x$names[-c(1:(length(iv)+1))]
    cat("\nThe DV (Y) was ", dv,". The IV (X) was ", iv,". The mediating variable(s) = ", mv,".")
    if(!is.null(x$mod)) cat("  The moderating variable(s) = ", x$names[x$mod])
    
   for(i in 1:length(iv)) { cat("\n\nTotal Direct effect(c) of ",iv[i], " on ", dv," = ",round(x$total[i],digits), "  S.E. = ", round(x$se.bx[i],digits), " t direct = ",round(x$tt[i],digits), "  with probability = ", signif(x$probt[i],digits))
    cat("\nDirect effect (c') of ",iv[i],  " on ", dv," removing ", mv ," = ",round(x$direct[i],digits), "  S.E. = ", round(x$se.beta[i],digits), " t direct = ",round(x$t[i],digits), "  with probability = ", signif(x$prob[i],digits))
     
   if(is.null(x$mod)) { cat("\nIndirect effect (ab) of ",iv[i], " on ", dv," through " ,mv , "  = ", round(x$indirect[i],digits),"\n")
   cat("Mean bootstrapped indirect effect = ",round(x$mean.boot[i],digits), " with standard error = ",round(x$sd.boot[i],digits), " Lower CI = ",round(x$ci.quant[1,i],digits), "   Upper CI = ", round(x$ci.quant[2,i],digits))}
     }
   if(is.null(x$mod)) {
     cat("\nSummary of a, b, and ab estimates and ab confidence intervals\n")
     if(!any(is.na(x$abc))) {print(round(x$abc,digits))} else {
      cat("\n 'a'  paths \n")
      print(round(x$a,digits))
      cat("\n'b' paths \n")
      print(round(x$b,digits))
      cat("\n'ab' paths \n")
      print(round(x$ab,digits))
     }
    
    cat("\nratio of indirect to total effect=  ", round(x$ratit,digits))
    cat("\nratio of indirect to direct effect= ", round(x$ratid,digits))
    }  else {
    cat("\nEffect of interaction of ",iv[1], " with ", iv[2] , "  = ", round(x$direct[3],digits),"  S.E. = ", round(x$se.beta[3,1],digits), " t direct = ",round(x$t[3,1],digits), "  with probability = ", signif(x$prob[3,1],digits))
    cat("\nIndirect effect due to interaction  of ",iv[1], " with ", iv[2] , "  = ", round(x$int.ind,digits))
    cat("\nMean bootstrapped indirect interaction effect = ",round(x$mean.boot[1],digits), " with standard error = ",round(x$sd.boot[1],digits), " Lower CI = ",round(x$ci.quant[1],digits), "   Upper CI = ", round(x$ci.quant[2,i],digits))
    cat("\nSummary of a, b, and ab estimates and ab confidence intervals\n")
   if(!is.na(x$abc)) {print(round(x$abc,digits))} else {
      print(round(x$a,digits))
      print(round(x$b,digits))
      print(round(x$ab,digits))
     }
   
    cat("\nR2 of model = ", round(x$R2,digits)) 
   }}
##end of copy 
)}
   
   
"mediate.diagram" <- function(medi,digits=2,ylim=c(2,8),xlim=c(-1,10),main="Mediation model",...) { 

 dv <- medi$names[1]
    iv <- as.matrix(rownames(medi$direct))
    mv <- medi$names[-c(1:(length(iv)+1))]
    C <- round(medi$C[c(iv,mv,dv),c(iv,mv,dv)],digits)
    

miny <- 5 - max(length(iv),length(mv)) - 1
maxy <- 5 + max(length(iv),length(mv)) + 1
if(missing(ylim)) ylim=c(miny,maxy)
plot(NA,xlim=xlim,ylim=ylim,main=main,axes=FALSE,xlab="",ylab="")
var.names <- c(rownames(medi$direct),colnames(medi$direct),rownames(medi$b))
n.mediate <- length(mv)

m <- list()
 c <- as.matrix(round(medi$total,digits))
 a <- as.matrix(round(medi$a,digits))
 if(ncol(a)==1) a <- t(a)
 b <- as.matrix(round(medi$b,digits))
 cprime <- as.matrix(round(medi$direct,digits))
 numx <- length(cprime)
 x <- list()
 
 if((numx > 1) && (n.mediate > 1) ) {adj <- 3} else {adj <- 2}  #this fixes where to put the labels on the a path

viv <- 1:numx
for(i in 1:numx)  {
if((numx %% 2)== 0) {
viv[i] <- switch(i,7,3,6,4,8,2,9,1,10)  } else { viv[i] <- switch(i,5,7,3,6,4,8,2,9)}
 x[[i]]  <- dia.rect(1,viv[i],iv[i])}
y <- dia.rect(9,5,dv) 

v.loc <- 1:n.mediate
for (mediate in 1:n.mediate) {
  if((n.mediate %% 2) ==0) {v.loc[mediate] <- switch(mediate,7,3,9,1,6,4,7,3,10) } else {
if(numx==1) {v.loc[mediate] <- switch(mediate,7,3,9,1,6,4,8,2)} else {
             v.loc[mediate] <- switch(mediate,5,3,7,2,6,4,8,2)}}
}


v.loc <- sort(v.loc,decreasing=TRUE)
if(n.mediate==1) a <- t(a)
for (mediate in 1:n.mediate) {
m[[mediate]] <- dia.rect(5,v.loc[mediate],mv[mediate] ) 

for(i in 1: numx) {dia.arrow(x[[i]]$right,m[[mediate]]$left,labels=a[i,mediate],adj=adj,...) #a term
      dia.arrow(x[[i]]$right,y$left,labels=paste("c = ",c[i]),pos=3,...)
      dia.arrow(x[[i]]$right,y$left,labels=paste("c' = ",cprime[i]),pos=1,...)}
      dia.arrow(m[[mediate]]$right,y$left,labels=b[mediate],...)  #     
}
rviv <- max(viv)
if(numx >1) {
  for (i in 2:numx) {
  for (k in 1:(i-1)) {dia.curved.arrow(x[[i]]$left,x[[k]]$left,C[i,k],scale=-(numx-1)*(abs(viv[i]-viv[k])/rviv),both=TRUE)} 
  } }
}	


"moderate.diagram" <- function(medi,digits=2,ylim=c(2,8),main="Moderation model",...) {
  
xlim=c(0,10)

plot(NA,xlim=xlim,ylim=ylim,main=main,axes=FALSE,xlab="",ylab="")
var.names <- rownames(medi$direct)
x.names <- rownames(medi$direct)
y.names <- colnames(medi$direct)
beta <- round(medi$direct,digits)
c <- round(medi$total,digits)
cprime <-round(medi$direct[1],digits) 
nx <- length(x.names)
ny <- length(y.names)  #should be 1
top <- max(nx,ny)
xlim=c(-nx/3,10)
ylim=c(0,top)
top <- max(nx,ny)
x <- list()
y <- list()
x.scale <- top/(nx+1)
y.scale <- top/(ny+1)
plot(NA,xlim=xlim,ylim=ylim,main=main,axes=FALSE,xlab="",ylab="")
for(i in 1:nx) {x[[i]] <- dia.rect(3,top-i*x.scale,x.names[i]) }
 
y[[1]] <- dia.rect(7,top-y.scale,y.names[1]) 
 dia.arrow(x[[1]]$right,y[[1]]$left,labels=paste("c = ",c),pos=3,...)
 dia.arrow(x[[1]]$right,y[[1]]$left,labels=paste("c' = ",cprime),pos=1,...)

 
for(i in 2:nx) {
   dia.arrow(x[[i]]$right,y[[1]]$left,labels = beta[i,1],adj=2)
   }
}  
     		 
     	
             
 