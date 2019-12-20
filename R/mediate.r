#The basic logic is to find the regression of Y on X, of M on X, and of Y on M
#and then compare the direct (Y on X) also known as c  to c - ab 
#Modified May, 2015 to allow multiple Xs (and multiple Ys)
#Revised June, 2015 for prettier output
#Revised February 2016 to get moderation to work
#The moderate part is a bit more complicated and needs to be cleaned up
#Three potential cases of moderate
#a) moderation with out mediation  -- not yet treated
#b) moderation of the independent to mediator  
#c) moderation of the mediator to dependent -- not yet treated
# 
#substantially revised May 6, 2016 by using matReg to make simpler to understand
#and November,2017 improved to handle formula input and do more moderations 
#added the zero option to match Hayes Process  
#added the ability to partial out variables
#12/1/18  Added the ability to do quadractic terms (similar to setCor) 
"mediate" <-
function(y,x,m=NULL, data, mod=NULL, z=NULL, n.obs=NULL,use="pairwise",n.iter=5000,alpha=.05,std=FALSE,plot=TRUE,zero=TRUE,main= "Mediation")  {
 
 #this first part just gets the variables right, depending upon the way we input them
   cl <- match.call()
   #convert names to locations 
   #first, see if they are in formula mode  
   if(inherits(y,"formula")) {  ps <- fparse(y)
   y <- ps$y
   x <- ps$x
   m <- ps$m #but, mediation is not done here, so we just add this to x
  # if(!is.null(med)) x <- c(x,med)   #not  necessary, because we automatically put this in
   mod <- ps$prod
   ex <- ps$ex  #the quadratic term
   #but, we want to drop the m variables from x
   x <- x[!ps$x %in% ps$m] 
   z <- ps$z  #are there any variables to partial
} else {ex <- NULL}
   
   all.ab <- NULL  #preset this in case we are just doing regression 
   
    if(is.numeric(y )) y <- colnames(data)[y]
    if(is.numeric(x )) x <- colnames(data)[x]
    if(!is.null(m))  if(is.numeric(m )) m <- colnames(data)[m]
    if(!is.null(mod) ) {if(is.numeric(mod)) {nmod <- length(mod)  #presumably 1 
                                         mod <- colnames(data)[mod] } }
    if(is.null(mod)) {nmod<- 0} else {nmod<- length(mod)}
     var.names <- list(IV=x,DV=y,med=m,mod=mod,z=z,ex=ex)
  
  if(any(!(unlist(var.names) %in% colnames(data)))) {stop ("Variable names not specified correctly")}                                      
    if(ncol(data) == nrow(data)) {raw <- FALSE 
            if(nmod > 0) {stop("Moderation Analysis requires the raw data") } else {data <- data[c(y,x,m,z),c(y,x,m,z)]} 
                 } else { data <- data[,c(y,x,m,z,ex)]
                
                  }
   # if(nmod > 0 ) {data <- data[,c(y,x,mod,m)] } else {data <- data[,c(y,x,m)]} #include the moderation variable
   if(nmod==1) {mod<- c(x,mod)
     nmod <- length(mod) 
  }
  
  
   if(!is.matrix(data)) data <- as.matrix(data)
  if((dim(data)[1]!=dim(data)[2]))  {n.obs=dim(data)[1]   #this does not take into account missing data
                    if(!is.null(mod)) if(zero) data <- scale(data,scale=FALSE)  #0 center 
                    C <- cov(data,use=use)
                    raw <- TRUE
                    if(std) {C <- cov2cor(C)}   #use correlations rather than covariances
                     }  else {
                    raw <- FALSE
                    C <- data
   nvar <- ncol(C)
   
   
   if(is.null(n.obs)) {n.obs <- 1000
      message("The data matrix was a correlation matrix and the number of subjects was not specified. \n n.obs arbitrarily set to 1000")
                       }
   if(!is.null(m)) { # only if we are doing mediation (12/11/18)
      message("The replication data matrices were simulated based upon the specified number of subjects and the observed correlation matrix.")
       eX <- eigen(C)   #we use this in the bootstrap replications in the case of a correlation matrix
       data <- matrix(rnorm(nvar * n.obs),n.obs)
       
     
      data <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(data) )
      colnames(data) <- c(y,x,m,z)}
     }
 

  

  if ((nmod > 0 ) | (!is.null(ex))) {if(!raw) {stop("Moderation  analysis requires the raw data")
   } else { if(zero) {data <- scale(data,scale=FALSE)}
   }
   }
   
if (nmod > 0 ) {
              prods <- matrix(NA,ncol=length(ps$prod),nrow=nrow(data))
                     colnames(prods) <- paste0("V",1:length(ps$prod))
                     for(i in 1:length(ps$prod)) {
                       prods[,i] <- apply(data[,ps$prod[[i]]],1,prod) 
                       colnames(prods)[i] <- paste0(ps$prod[[i]],collapse="*")
                       }
                      
                      data <- cbind(data,prods)
                      x <- c(x,colnames(prods))
                }
             
     if(!is.null(ex)) {
                 
                    quads <- matrix(NA,ncol=length(ex),nrow=nrow(data))  #find the quadratric terms
                    colnames(quads) <- ex
                     for(i in 1:length(ex)) {
                      quads[,i] <- data[,ex[i]] * data[,ex[i]]
                      colnames(quads)[i] <- paste0(ex[i],"^2")
                     }
                     data <- cbind(data,quads)
                     x <- c(x,colnames(quads))
                    }
   #We have now added in products and quadratics, if desired
                    
             if(raw)  {C <- cov(data,use=use) } #else {C <- data}
              if(std) {  C <- cov2cor(C)} 
             
              
  ######### 
  #########
  #now, we are ready to process the data
  #We do the basic regressions as matrix operations using matReg

        xy <- c(x,y)
        numx <- length(x) 
     	numy <- length(y)
     	if(!is.null(m)) {numm <- length(m)
        nxy <- numx + numy 
        m.matrix <- C[c(x,m),c(x,m),drop=FALSE] #this is the predictor +  moderator matrix
        } else {numm <- 0
             nxy <- numx}  #the case of moderation without mediation 
     	df <- n.obs - nxy - 1 

     	xy.matrix <- C[c(x,m),y,drop=FALSE]   #this is the matrix of correlations with the criterion

#this next section adds the intercept to the regressions 
 C.int <- matrix(NA,nrow=NROW(C)+1,ncol=ncol(C)+1)
  C.int[1,] <- C.int[,1] <- 0
  C.int[-1,-1] <-C
 if(!raw) { std <- TRUE
     C.int[1,1] <- 1
     means <- rep(0,NCOL(data))} else {
     C.int[1,1] <- 0 
    if(!std){ means <-colMeans(data,na.rm=TRUE)} else { means <- rep(0,ncol(C))}
     means <- c(1,means)
     C.int <- C.int *  (n.obs-1) + means %*% t(means) * n.obs
     } 
 rownames(C.int) <- colnames(C.int) <- c("Intercept",colnames(C))
if(std) { C.int <- cov2cor(C.int) }  

#first, find the complete regression model 
cprime.reg <- matReg(c("Intercept",x,m),y,C=C.int,n.obs=n.obs,z=z,means=means,std=std,raw=raw) #we have added  the intercept term here (11/25/19)

 ##this is the zero order beta -- the total effect
   total.reg <- matReg(x,y,m=m,z=z,C=C,n.obs=n.obs,std=std)   #include m for correct df, add in the z here to do partial correlations
   direct <- total.reg$beta
   if(!is.null(z)) {colnames(direct) <- paste0(colnames(direct),"*")
                   rownames(direct) <- paste0(rownames(direct),"*") }
   #There are 3 broad cases that need to be handled somewhat differently in terms of the matrix operations
   # 1 IV, 1 or more mv 
   # multiple IV, 1 MV
   #multiple IV, multiple MV
   #this is the direct path from X to M  
   #For the purposes of moderation, at least for now, think of this as just 2 or 3 IVs

 #get a, b and cprime effects and their se
  
 if(numm > 0) {a.reg <- matReg(x=c("Intercept",x),y=m,C=C.int,z=z,n.obs=n.obs,means=means,std=std)   #the default case is to have at least one mediator
 b.reg <- matReg(c(x,m),y,C=C,z=z,n.obs=n.obs)

 
# cprime.reg <- matReg(c("Intercept",x,m),y,C=C.int,n.obs=n.obs,z=z,means=means,std=std) #we have added  the intercept term here (11/25/19)

 a <- a.reg$beta[-1,,drop=FALSE]
 b <- b.reg$beta[-(1:numx),,drop=FALSE]
 c <- total.reg$beta
 cprime <- cprime.reg$beta
 

   # ab <- a * b    #these are the ab products for each a and b path   c' is c - sum of all of these
#  all.ab <- matrix(NA,ncol=numm*numy,nrow=numx)
#We currently only get the all.ab  terms if there is just 1 dv



 all.ab <- matrix(NA,ncol=numm*numx,nrow=numy)
 # if((numx == 1) & (numy==1)) {ab <- a * t(b)}  else {ab <- matrix(NA,ncol=numm*numx,nrow=numy)
#                                  for (i in 1:numx) { ab <- t(a[i,] *  b[,1:numy])}}
#fixed ? November 18, 2019 to handle more than 1 dv

  for(i in 1:numx) {
  if((numx == 1) & (numy==1)) {all.ab  <- a * t(b)}  else {
       all.ab <- matrix(NA,ncol=numm*numx,nrow=numy)
     for (i in 1:numx) { all.ab <- t(a[i,] *  b[,1:numy])}}               
  
 # all.ab[i,] <- a[i,] * t(b[,1])
 }   #just do the first column of b   (this is problematic, perhaps and doesn't work for two dependent variables)
 # colnames(all.ab) <- rep(m, numy)

 # colnames(all.ab) <- m
 # rownames(all.ab) <- x
   ab <- a %*% b      #are we sure that we want to do the matrix sum? 
   indirect <- c - ab
   

     
 if(is.null(n.obs)) {message("Bootstrap is not meaningful unless raw data are provided or the number of subjects is specified.") 
     
          mean.boot <- sd.boot <- ci.quant <- boot <-  se <- tvalue  <- prob <- NA } else {
         #  if(is.null(mod)) {
           boot <- boot.mediate(data,x,y,m,z,n.iter=n.iter,std=std,use=use)   #this returns a list of vectors
  #the first values are the indirect (c') (directly found), the later values are c-ab from  the products
            mean.boot <- colMeans(boot)
            sd.boot <- apply(boot,2,sd)
            ci.quant <- apply(boot,2, function(x) quantile(x,c(alpha/2,1-alpha/2),na.rm=TRUE)) 

           # mean.boot <- matrix(mean.boot[1:(numx*numy)],nrow=numx)
           mean.boot <- matrix(mean.boot,nrow=numx)
           # sd.boot <- matrix(sd.boot[1:(numx*numy)],nrow=numx)
           sd.boot <- matrix(sd.boot,nrow=numx)
          #  ci.ab <- matrix(ci.quant,nrow=2*numx)
          ci.ab <- ci.quant
           # colnames(mean.boot) <- colnames(sd.boot)  <- c(y,m)
            rownames(mean.boot) <- rownames(sd.boot) <- x
            
            boots <- list(mean=mean.boot,sd=sd.boot,ci=ci.quant,ci.ab=ci.ab)
            
           
              } 
      } else {    #the case of just an interaction term 
    a.reg <- b.reg <- reg <- NA
    a <- b <- c <- ab <- cprime <- boot<- boots <- indirect <- NA}  
              
    #beta.x is the effect without the mediators
    #direct is the effect with the mediators
    #indirect is ab from the difference
    #ab is ab from the product  of a and b paths

if(!is.null(z)) {var.names$IV <- paste0(var.names$IV,"*")
                var.names$DV <- paste0(var.names$DV,"*")
                var.names$med <- paste0(var.names$med,"*")
                colnames(C) <- rownames(C) <- paste0(colnames(C),"*")
               
                
}

         
     result <- list(var.names=var.names,a=a,b=b,ab=ab,all.ab = all.ab,c=c,direct=direct,indirect=indirect,cprime = cprime, total.reg=total.reg,a.reg=a.reg,b.reg=b.reg,cprime.reg=cprime.reg,boot=boots,boot.values = boot,sdnames=colnames(data),data=data,C=C,  Call=cl)
   
    class(result) <- c("psych","mediate")
       
   if(plot) {if(is.null(m)) {moderate.diagram(result) } else {  mediate.diagram(result,main=main) }
   }
  
    return(result)	
     		 }
   

#a helper function to find regressions from covariances
#May 6, 2016
#Fixed November 29, 2018 to handle se of partialed variables correctly
#modified September 25, 2019 to find intercepts and standard errors using momements
#modified even more November 25, 2019 to return the intercepts and R2
matReg <- function(x,y,C,m=NULL,z=NULL,n.obs=0,means=NULL,std=FALSE,raw=TRUE) {
 if(is.null(n.obs)) n.obs <- 0
   numx <- length(x)   #this is the number of predictors (but we should adjust by the number of covariates)   
   numz <- length(z)
   numy <- length(y)
   #df <- n.obs -1 - numx - length(z) - length(m)    #but this does not take into account the mediating variables
   #note that the x variable includes the intercept and thus uses up one extra df
   df <- n.obs  - numx -numz #We have partialed out z, should we use the df from it?  This is changed 11/26/19 to reduce df for z  
   Cr <- cov2cor(C)  
        	if(!is.null(z)){numz <- length(z)      #partial out the z variables
     	                zm <- C[z,z,drop=FALSE]
     	                za <- C[x,z,drop=FALSE]
     	                zb <- C[y,z,drop=FALSE]
     	                zmi <- solve(zm)
     	                 x.matrix <- C[x,x,drop=FALSE] - za %*% zmi %*% t(za)
     	                 y.matrix <- C[y,y,drop=FALSE] - zb %*% zmi %*% t(zb)
     	                xy.matrix <- C[x,y,drop=FALSE] - za  %*% zmi %*% t(zb)
     	                 C <- cbind(rbind(y.matrix,xy.matrix),rbind(t(xy.matrix),x.matrix))
     	                
     	                 }
   
    if(numx==1) { beta <- solve(C[x,x,drop=FALSE],(C[x,y,drop=FALSE])) 
                 colnames(beta) <- y
        } else {
        beta <- solve(C[x,x],(C[x,y])) }    #this is the same as setCor and is a x * x matrix
        
    if(!is.matrix(beta)) {beta <- matrix(beta,nrow=length(beta))}   #beta is a matrix of beta weights 
    if(is.character(x)) {rownames(beta) <- x}  else {rownames(beta) <- colnames(C)[x]}
    if(is.character(y)) { colnames(beta) <- y} else { colnames(beta) <- colnames(C)[y]}
      
      x.inv <- solve(C[x,x]) #solve x.matrix    #taken from setCor
      yhat <- t(C[x,y,drop=FALSE]) %*% x.inv %*% C[x,y,drop=FALSE]
      resid <- C[y,y]- yhat
      
     if(!std ) {
        df <- n.obs - numx - numz

         Residual.se <- sqrt(diag(resid /df))  #this is the df  n.obs - length(x))
         se <- MSE <- diag(resid )/(df) 
             
               if(length(y) > 1) {SST <- diag(C[y,y] - means[y]^2 * n.obs)} else {SST <- ( C [y,y] - means[y]^2 * n.obs)}
     	       R2 <- (SST - diag(resid)) /SST 
     	       se.beta <- list()
     	        for (i in 1:length(y)) {
     	        se.beta[[i]] <- sqrt(MSE[i] * diag(x.inv))
                 }
                se <- matrix(unlist(se.beta),ncol=numy)
      if(length(y) > 1) {SST <- diag(C [y,y] - means[y]^2 * n.obs)} else {SST <- ( C [y,y] - means[y]^2 * n.obs)}
	       R2 <- (SST - diag(resid)) /SST 
	       
	   
      } else {      R2 <- colSums(beta * C[x,y])/diag(C[y,y,drop=FALSE])   #the standardized case
        uniq <- 1-(1-1/diag(solve(Cr[x,x,drop=FALSE])))  #1- smc
     

        if(n.obs > 2) { # se <- (sqrt((1-R2)/(n.obs-1 - numx-numz)) %*% t(sqrt(1/uniq)))  #these are the standardized se
                        se <- (sqrt((1-R2)/(df)) %*% t(sqrt(1/uniq)))   #setCor uses df = n.obs - numx - 1
                        se <- t( se * sqrt(diag(C[y,y,drop=FALSE])) %*% t(sqrt(1/diag(C[x,x,drop=FALSE]))) )  #But does this work in the general case?

        
                    colnames(se) <- colnames(beta) } else {se <- NA}
                   if(raw) {   #used to compare models  -- we need to adjust this for dfs
                           Residual.se <-   sqrt((1-R2)* df/(df-1)) } else {  #this is a kludge and is necessary to treat the SSR correctly
                       Residual.se <- sqrt((1-R2)/df * (n.obs-1))}
                   }
                   
                
                    if(!any(is.na(se))) { tvalue <- beta/se
                                        # prob <- 2*(1- pt(abs(tvalue),df))
                                          prob <- -2 *  expm1(pt(abs(tvalue),df,log.p=TRUE))
                                         } else {tvalue <- prob <- df <- NA}
  result <- list(beta=beta,se=se, t=tvalue,df=df,prob=prob,R2=R2,SE.resid=Residual.se)
  return(result)       }  		 
#######




partialReg <- function(C,x,y,m,z) {
   x <- c(x,m)
   y <- c(y,m)
                        numz <- length(z)      #partial out the z variables
     	                zm <- C[z,z,drop=FALSE]
     	                za <- C[x,z,drop=FALSE]
     	                zb <- C[y,z,drop=FALSE]
     	                zmi <- solve(zm)
     	                 x.matrix <- C[x,x,drop=FALSE] - za %*% zmi %*% t(za)
     	                 y.matrix <- C[y,y,drop=FALSE] - zb %*% zmi %*% t(zb)
     	                xy.matrix <- C[x,y,drop=FALSE] - za  %*%zmi %*% t(zb)
     	                 C <- cbind(rbind(y.matrix,xy.matrix),rbind(t(xy.matrix),x.matrix))
     	                 }
#finally fixed November 3, 2019
boot.mediate <- function(data,x,y,m,z,n.iter=10,std=FALSE,use="pairwise") {

  n.obs <- nrow(data)
  numx <- length(x)
  numy <- length(y)
  numm <- length(m)
  nxy <- numx + numy 
 result <- matrix(NA,nrow=n.iter,ncol = (numx*numy+ numm*numy*numx ))
 if((numm > 1) &  (numx > 1)) ab <- matrix(0,nrow=numx,ncol=numy)     
  for (iteration in 1:n.iter) {
  
  samp.data <- data[sample.int(n.obs,replace=TRUE),]
  C <- cov(samp.data,use=use)
  if(!is.null(z)) C <- partialReg(C,x,y,m,z)  #partial out z
  if(std) C <- cov2cor(C)
   xy <- c(x,y)     
        m.matrix <- C[c(x,m),c(x,m)]
     #	df <- n.obs - nxy - 1 
     	xy.matrix <- C[c(x,m),y,drop=FALSE]
     	
   if(numx ==1) { beta.x <- solve(C[x,x],t(C[x,y]) ) } else  {beta.x <- solve(C[x,x],C[x,y]) }  #this is the zero order beta -- the total effect
  # if(numx ==0) { a <-  solve(C[x,x,drop=FALSE],t(C[x,m,drop=FALSE]) ) } else { a <-  solve(C[x,x,drop=FALSE],(C[x,m,drop=FALSE]) )} #the a paths    
   a <-  solve(C[x,x,drop=FALSE],(C[x,m,drop=FALSE]))
    beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
    beta <- as.matrix(beta.xm)     #these are the individual predictors, x and then m

  b <- beta[-c(1:numx),,drop=FALSE]

 if((numx == 1) & (numy==1)) {ab <- a[-1] * t(b)}  else {ab <- array(NA,dim=c(numx,numm,numy))
                                 for(j in 1:numy) {
                                 for(k in 1:numm) {
                                  for (i in 1:numx) { ab[i,k,j] <- t(a[i,k] *  b[k,j])}}   #this needs to be fixed for two numm>1
                                  }}
                                  
 # ab <- a %*% b    #each individual path  #this probably is only correct for the numx = 1 model
#  if((numx > 1) & (numy > 1)) {for (i in 1:numx) {ab[i,] <- a[i,] * b}}
  # ab <- a * b  #we don't really need this
#all.ab <- matrix(NA,ncol=numm*numy,nrow=numx*numm)   
   #The number of all.ab terms is ( numx * numm) * (numm * numy)
 # for(j in 1:numm*numy) {
 # for(i in 1:numx*numm) {all.ab[i,j] <- a[i,] * t(b[j,i])}  #this just does one column of b
  #
  
  #consider muliple ivs 
 # all.ab <-outer(a,t(b))    #this is cute,but actually not correct
  #all.ab <- matrix(all.ab,nrow=numx*numm,ncol=numm*numy)
  all.ab <- ab
  indirect <-  beta.x - beta[1:numx,1:numy]  #this is c' = c - ab
 # result[iteration,] <- c(indirect,ab)  #this is a list of vectors

  result[iteration,] <- c(indirect, all.ab)  #this is a list of vectors -- do we really need all all.ab since it is identical to indirect?
   }
 return(result)  }

 
   
"mediate.diagram" <- function(medi,digits=2,ylim=c(3,7),xlim=c(-1,10),show.c=TRUE, main="Mediation model",cex=1,l.cex=1,...) { 
    if(missing(l.cex)) l.cex <- cex
     dv <- medi$var.names[["DV"]]
   # iv <- medi$var.names[["IV"]]
     iv <- as.matrix(rownames(medi$direct))
    mv <- medi$var.names[["med"]]
    mod <- medi$var.names[["mod"]]
   # numx <- length(medi$var.names[["IV"]])
   numx <- NROW(iv)
    numy <- length(dv)
    direct <- round(medi$direct,digits)
    C <- round(medi$C[c(iv,mv,dv),c(iv,mv,dv)],digits)
#if have moderated effects, this is the same as more xs    

miny <- 5 - max(length(iv)/2,length(mv),2) - .5
maxy <- 5 + max(length(iv)/2,length(mv),2) + .5

if(missing(xlim)) xlim=c(-numx * .67,10)
if(missing(ylim)) ylim=c(miny,maxy)
plot(NA,xlim=xlim,ylim=ylim,main=main,axes=FALSE,xlab="",ylab="")
var.names <- c(rownames(medi$direct),colnames(medi$direct),rownames(medi$b))
if(is.null(mv)) {n.mediate <- 0} else {n.mediate <- length(mv)}

m <- list()
 #c <- as.matrix(round(medi$total,digits))
 c <- as.matrix(round(medi$c,digits))
 a <- as.matrix(round(medi$a,digits))
 if(ncol(a)==1) a <- t(a)
 b <- as.matrix(round(medi$b,digits))
 cprime <- as.matrix(round(medi$cprime,digits))
 
 x <- list()
 
 if((numx > 1) && (n.mediate > 1) ) {adj <- 3} else {adj <- 2}  #this fixes where to put the labels on the a path

viv <- 1:numx
for(i in 1:numx)  {
if((numx %% 2)== 0) {
viv[i] <- switch(i,7,3,6,4,8,2,9,1,10)  } else { viv[i] <- switch(i,5,7,3,6,4,8,2,9)}
 x[[i]]  <- dia.rect(1,viv[i],iv[i],cex=cex,...)}
 
 vdv <- 1:numy
 y <- list()
 for (i in 1:numy) {
 if((numy %% 2)== 0) {
vdv[i] <- switch(i,6,4,7,3,8,2,9,1,10)  } else { vdv[i] <- switch(i,5,7,3,6,4,8,2,9)}
 y[[i]] <- dia.rect(9,vdv[i],dv[i],cex=cex,...) 
 }

#y <- dia.rect(9,5,dv) 

v.loc <- 1:n.mediate
if(n.mediate > 0) {
for (mediate in 1:n.mediate) {
  if((n.mediate %% 2) ==0) {v.loc[mediate] <- switch(mediate,7,3,9,1,6,4,7,3,10) } else {
    switch(numx,
   1:  {v.loc[mediate] <- switch(mediate,7,3,8,1,6,4,9,2)}, 
   2 : {v.loc[mediate] <- switch(mediate,5,3,7,2,6,4,8,2)},
   3:  {v.loc[mediate] <- switch(mediate,5.5,3,7,2,5,4,8,2)},
   4:  {v.loc[mediate] <- switch(mediate,5,3,7,2,6,4,8,2)},
   5:  {v.loc[mediate] <- switch(mediate,6,3,7,2,5,4,8,2)},
   6:  {v.loc[mediate] <- switch(mediate,5,3,7,2,6,4,8,2)},
   7:   {v.loc[mediate] <- switch(mediate,6,3,7,2,5,4,8,2)})
}
}
}


v.loc <- sort(v.loc,decreasing=TRUE)
if(n.mediate ==0) { for(j in 1:numy) {
			for(i in 1: numx) {
     			 dia.arrow(x[[i]]$right,y[[j]]$left,labels=paste("c = ",direct[i,j]),pos=0,cex=l.cex,...)}
	}  
} else {
if(n.mediate==1) a <- t(a)
for (mediate in 1:n.mediate) {
	m[[mediate]] <- dia.rect(5,v.loc[mediate],mv[mediate],cex=cex,... ) 
		for(j in 1:numy) {
			for(i in 1: numx) {dia.arrow(x[[i]]$right,m[[mediate]]$left,labels=a[i,mediate],adj=adj,cex=l.cex,...) #a term
     			if(show.c) {dia.arrow(x[[i]]$right,y[[j]]$left,labels=paste("c = ",c[i,j]),pos=3,cex=l.cex,...)}
     			 dia.arrow(x[[i]]$right,y[[j]]$left,labels=paste("c' = ",cprime[i+1,j]),pos=1,cex=l.cex,...)}   #we have an intercept
      			dia.arrow(m[[mediate]]$right,y[[j]]$left,labels=b[mediate,j],cex=l.cex,...)  #     
			}
		} 
		}
	rviv <- max(viv)
if(numx >1) {
  for (i in 2:numx) {
  for (k in 1:(i-1)) {dia.curved.arrow(x[[i]]$left,x[[k]]$left,C[i,k],scale=-(numx-1)*(abs(viv[i]-viv[k])/rviv),both=TRUE,dir="u",cex=l.cex,...)} 
  } }
	}


"moderate.diagram" <- function(medi,digits=2,ylim=c(2,8),main="Moderation model",cex=1,l.cex=1,...) {
 
if(missing(l.cex)) l.cex <- cex  
xlim=c(0,10)

#plot(NA,xlim=xlim,ylim=ylim,main=main,axes=FALSE,xlab="",ylab="")
var.names <- rownames(medi$direct)
x.names <- rownames(medi$direct)
y.names <- colnames(medi$direct)
beta <- round(medi$direct,digits)

nx <- length(x.names)
ny <- length(y.names)
top <- max(nx,ny)
xlim=c(-nx/3,10)
ylim=c(0,top)
top <- max(nx,ny)
x <- list()
y <- list()
x.scale <- top/(nx+1)
y.scale <- top/(ny+1)

plot(NA,xlim=xlim,ylim=ylim,main=main,axes=FALSE,xlab="",ylab="",...)



for(i in 1:nx) {x[[i]] <- dia.rect(2,top-i*x.scale,x.names[i],cex=cex,...) }
for(j in 1: ny) {y[[j]] <- dia.rect(7,top-j*y.scale,y.names[j],cex=cex,...) }
y[[1]] <- dia.rect(7,top-y.scale,y.names[1],cex=cex,...) 
# dia.arrow(x[[1]]$right,y[[j]]$left,labels=paste("c = ",c),pos=3,...)
#dia.arrow(x[[1]]$right,y[[j]]$left,labels=paste("c' = ",cprime),pos=1,...)

for(j in 1:ny){
for(i in 1:nx) {
   dia.arrow(x[[i]]$right,y[[j]]$left,labels = beta[i,j],adj=2,cex=l.cex,...)
   }
   }
if(nx >1) {
  for (i in 2:nx) {
  for (k in 1:(i-1)) {dia.curved.arrow(x[[i]]$left,x[[k]]$left,round(medi$C[i+1,k+1],2),scale= -(abs(k-i)),both=TRUE,dir="u",cex=l.cex,...)} 
  } }

}  
     		 
  #finally got the print to work on multiple dvs    11/24/19 	
  "summary.psych.mediate" <- function(x,digits=2,short=FALSE) {
 cat("Call: ")
    print(x$Call)
    dv <- x$var.names[["DV"]]
   # iv <- x$var.names[["IV"]]
    mv <- x$var.names[["med"]]
    mod <- x$var.names[["mod"]]
   # dv <- x$names[1]
    iv <- rownames(x$direct)
    niv <- length(iv)
    nmed <- length(mv)
    ndv <- length(dv)
     nz <- length(x$var.names[["z"]])
     
 if(nmed < 1) {
 cat("\nNo mediator specified leads to traditional regression \n") } else {    
 cat("\nDirect effect estimates (traditional regression)    (c') \n")}
 
     for(j in 1:ndv) {
     
    if (niv==1) { dfd <- round(data.frame(direct=x$cprime.reg$beta[,j],se = x$cprime.reg$se[,j],t=x$cprime.reg$t[,j],df=x$cprime.reg$df),digits)
     dfdp <- cbind(dfd,p=signif(x$cprime.reg$prob[,j],digits=digits+1)) } else {
     dfd <- round(data.frame(direct=x$cprime.reg$beta[1:(niv+1+nmed),j],se = x$cprime.reg$se[1:(niv+1+nmed),j],t=x$cprime.reg$t[1:(niv+1+nmed),j],df=x$cprime.reg$df),digits)
     dfdp <- cbind(dfd,p=signif(x$cprime.reg$prob[1:(niv+1+nmed),j],digits=digits+1))
     }
      colnames(dfdp) <- c(dv[j],"se","t","df","Prob")
     
   print(dfdp)
     F <-  x$cprime.reg$df * x$cprime.reg$R2[j]/(((nrow(x$cprime.reg$beta)-1) * (1-x$cprime.reg$R2[j])))
      pF <-  -expm1(pf(F,nrow(x$cprime.reg$beta)-1,x$cprime.reg$df,log.p=TRUE)) 
      cat("\nR =", round(sqrt(x$cprime.reg$R2[j]),digits),"R2 =", round(x$cprime.reg$R2[j],digits),  "  F =", round(F,digits), "on",nrow(x$cprime.reg$beta)-1, "and", x$cprime.reg$df,"DF   p-value: ",signif(pF,digits+1), "\n") 
   }
 
 if(nmed > 0) {
     cat("\n Total effect estimates (c) \n")
      
        for(j in 1:ndv) {

    dft <- round(data.frame(direct=x$total.reg$beta[,j],se = x$total.reg$se[,j],t=x$total.reg$t[,j],df=x$total.reg$df),digits)
    dftp <- cbind(dft,p = signif(x$total.reg$prob[,j],digits=digits+1))
    colnames(dftp) <- c(dv[j],"se","t","df","Prob")
    rownames(dftp) <- rownames(x$total.reg$beta)
     print(dftp)
    }
   

   
    
     
    cat("\n 'a'  effect estimates \n")

  if(niv==1) {
    for(i in 1:nmed) {
    	dfa <- round(data.frame(a = x$a.reg$beta[,i],se = x$a.reg$se[,i],t = x$a.reg$t[,i],df= x$a.reg$df),digits)
    	dfa <- cbind(dfa,p=signif(x$a.reg$prob[,i],digits=digits+1))
    	if(NROW(dfa) ==1) {rownames(dfa) <- rownames(x$a.reg$beta)
    	colnames(dfa) <-  c(colnames(x$a.reg$beta),"se","t","df", "Prob")} else {
    	rownames(dfa) <- rownames(x$a.reg$beta)
    	colnames(dfa) <-  c(colnames(x$a.reg$beta)[i],"se","t","df", "Prob")}
    	
    	print(dfa)} 
    	} else {
    	
     	for (i in 1:nmed) {
     	dfa <- round(data.frame(a = x$a.reg$beta[,i],se = x$a.reg$se[,i],t = x$a.reg$t[,i],df= x$a.reg$df),digits)
    	dfa <- cbind(dfa,p=signif(x$a.reg$prob[,i],digits=digits+1))
     	rownames(dfa) <-rownames(x$a.reg$beta)
     	colnames(dfa) <-  c(colnames(x$a.reg$beta)[i],"se","t","df","Prob") 
     	print(dfa) }
     	
     	}
     	        
      cat("\n 'b'  effect estimates \n")
      for (j in 1:ndv) {
      if(niv==1) {
     dfb <- round(data.frame(direct=x$b.reg$beta[-(1:niv),j],se = x$b.reg$se[-(1:niv),j],t=x$b.reg$t[-(1:niv),j], df=x$b.reg$df),digits)
     dfb <- cbind(dfb,p=signif(x$b.reg$prob[-(1:niv),j],digits=digits+1))} else {
      dfb <- round(data.frame(direct=x$b.reg$beta[-(1:niv),j],se = x$b.reg$se[-(1:niv),j],t=x$b.reg$t[-(1:niv),j],df=x$b.reg$df),digits)
     dfb <- cbind(dfb,p=signif(x$b.reg$prob[-(1:niv),j],digits=digits+1))}
     rownames(dfb) <- rownames(x$b.reg$beta)[-(1:niv)]
     colnames(dfb) <-  c(dv[j],"se","t","df", "Prob")
      print(dfb)
      }
 
      cat("\n 'ab'  effect estimates (through mediators)\n")
     
 for (j in 1:ndv) {
     
      dfab  <-round(data.frame(indirect = x$ab[,j],boot = x$boot$mean[,j],sd=x$boot$sd[,j],
                        #  lower=x$boot$ci[1,1:niv],   #was niv perhaps should be ndv?
                         #  upper=x$boot$ci[2,1:niv]),digits)
                           lower=x$boot$ci[1,(1:niv + niv*(j-1))],
                           upper=x$boot$ci[2,(1:niv + niv*(j-1))]),digits)
                         #  lower=x$boot$ci[1,(j*niv )],
                          # upper=x$boot$ci[2,(j*niv )]),digits)
      rownames(dfab) <- rownames(x$ab)
      colnames(dfab)[1] <- dv[j]
      print(round(dfab,digits))
      }
      
    #now show the individual ab effects (just works for 1 dv)
for(k in 1: ndv) {
    if(nmed > 1) {
    cat("\n 'ab' effects estimates for each mediator for",colnames(x$ab)[k], "\n")
    for (j in 1:nmed) {
       dfab  <-round(data.frame(#indirect = x$all.ab[,j],
             boot = x$boot$mean[,j+ndv*k],sd=x$boot$sd[,j+ndv*k],
                           lower=x$boot$ci[1,(j*niv*k +1):(j*niv*k +niv)],
                           upper=x$boot$ci[2,(j*niv*k +1):(j*niv*k +niv)]),digits)
      rownames(dfab) <- rownames(x$ab)
      colnames(dfab)[1] <- mv[j]
      print(round(dfab,digits))
      }
     }
    }
    }
    } 
          
 