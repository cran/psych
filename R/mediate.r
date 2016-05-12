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
"mediate" <-
function(y,x,m=NULL,data, mod=NULL,n.obs=NULL,use="pairwise",n.iter=5000,alpha=.05,std=FALSE,plot=TRUE)  {
 
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
   if(!is.null(m))  if(is.numeric(m )) m <- colnames(data)[m]
    if(!is.null(mod) ) {if(is.numeric(mod)) {nmod <- length(mod)  #presumably 1 
                                         mod <- colnames(data)[mod] } }
    if(is.null(mod)) {nmod<- 0} else {nmod<- length(mod)}
                                         
    if(ncol(data) == nrow(data)) {raw <- FALSE 
            if(nmod > 0) {stop("Moderation Analysis requires the raw data") } else {data <- data[c(y,x,m),c(y,x,m)]} 
                 } else {
    if(nmod > 0 ) {data <- data[,c(y,x,mod,m)] } else {data <- data[,c(y,x,m)]} #include the moderation variable
   
  }
  
  var.names <- list(IV=x,DV=y,med=m,mod=mod)
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
      message("The data matrix was a correlation matrix and the number of subjects was not specified. \n The replication data matrices were simulated based upon the observed correlation matrix and  n.obs set to 1000")
   } else { message("The replication data matrices were simulated based upon the specified number of subjects and the observed correlation matrix.")}
       eX <- eigen(C)   #we use this in the simulation in the case of a correlation matrix
      data <- matrix(rnorm(nvar * n.obs),n.obs)
     data <- t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(data) )
     colnames(data) <- c(y,x,m)
  }
if(std) {C <- cov2cor(C)}   #use correlations rather than covariances
  
 
  if(nmod > 0 ) {if(!raw) {stop("Moderation analysis requires the raw data")
   } else {ivXm <- matrix(data[,x] * data[,mod],ncol=length(mod)) #add in a moderating variable as a product term
       
       colnames(ivXm) <- paste0(abbreviate(x),"X",abbreviate(mod))
      data <- cbind(data,ivXm)
        
             # if(!(mod %in% m)) {m <- c(m,mod,colnames(ivXm))} else {m <- c(m,colnames(ivXm))}  #this is questionable C
              if(nmod > 0)   {ivX <- c(x,mod,colnames(ivXm))
              x <- ivX  #the x variables now include the moderator and the x * mod products
              var.names$IV <- x  
              } 
              C <- cov(data,use=use)
              if(std) {  C <- cov2cor(C)}
              }
              }
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

     	xy.matrix <- C[c(x,m),y,drop=FALSE]   #this is the matrix of correlations with the criterion (still just one)

   ##this is the zero order beta -- the total effect
   total.reg <- matReg(x,y,C,n.obs)
   direct <- total.reg$beta

   #There are 3 broad cases that need to be handled somewhat differently in terms of the matrix operations
   # 1 IV, 1 or more mv 
   # multiple IV, 1 MV
   #multiple IV, multiple MV
   #this is the direct path from X to M  
   #For the purposes of moderation, at least for now, think of this as just 2 or 3 IVs

 #get a, b and cprime effects and their se
  
 if(numm > 0) {a.reg <- matReg(x,m,C,n.obs)   #the default case is to have at least one mediator
 b.reg <- matReg(c(x,m),y,C,n.obs)
 cprime.reg <- matReg(c(x,m),y,C,n.obs) 

 a <- a.reg$beta
 b <- b.reg$beta[-(1:numx),,drop=FALSE]
 c <- total.reg$beta
 cprime <- cprime.reg$beta
 

   # ab <- a * b    #these are the ab products for each a and b path   c' is c - sum of all of these
   ab <- a %*% b
   indirect <- c - ab
   
#    if((numx > 1)&& (numy > 1)) {for (i in 1:numx) {ab[i,] <- a[i,] * b}}

  

     
     if(is.null(n.obs)) {message("Bootstrap is not meaningful unless raw data are provided or the number of subjects is specified.") 
     
          mean.boot <- sd.boot <- ci.quant <- boot <-  se <- tvalue  <- prob <- NA } else {
         #  if(is.null(mod)) {
           boot <- boot.mediate(data,x,y,m,n.iter=n.iter,std=std,use=use)   #this returns a list of vectors
                                                                            #the first x elements are indirect, the last x * m are ab effects
    #        }  else {
      #         boot <- boot.moderate( data,x,y,m,mod,n.iter=n.iter,std=std,use=use)
     #       }
      
            mean.boot <- colMeans(boot)
            sd.boot <- apply(boot,2,sd)
            ci.quant <- apply(boot,2, function(x) quantile(x,c(alpha/2,1-alpha/2),na.rm=TRUE)) 
            mean.boot <- matrix(mean.boot)
            ci.ab <- matrix(ci.quant,nrow=2*numx * numy)
            boots <- list(mean=mean.boot,sd=sd.boot,ci=ci.quant,ci.ab=ci.ab)
              } 
      } else {    #the case of just an interaction term 
    a.reg <- NA
    b.reg <- NA
    c.reg <- NA
    a <- b <- c <- ab <- cprime <- boot<- boots <- indirect <- cprime.reg <- NA}  
              
    #beta.x is the effect without the mediators
    #direct is the effect with the mediators
    #indirect is ab from the difference
    #ab is ab from the product  of a and b paths

 #    if(n.obs > 1) {if(ncol(mean.boot) == ncol(beta)) {colnames(mean.boot) <- colnames(beta)} else {colnames(mean.boot)[1:ncol(beta)] <- colnames(mean.boot)}  #temp
 #   rownames(b) <- rownames(beta)[-c(1:numx)]  #the ivs (including the moderators) are the first elements in beta
  #             }
     result <- list(var.names=var.names,a=a,b=b,ab=ab,c=c,direct=direct,indirect=indirect,cprime = cprime, total.reg=total.reg,a.reg=a.reg,b.reg=b.reg,cprime.reg=cprime.reg,boot=boots,boot.values = boot,sdnames=colnames(data),C=C,  Call=cl)
   
    class(result) <- c("psych","mediate")
   if(plot) { if(is.null(mod)) { mediate.diagram(result) } else {mediate.diagram(result,main="Moderation model")}
   }
  
    return(result)	
     		 }
   

#a helper function to find regressions from covariances
#May 6, 2016
matReg <- function(x,y,C,n.obs=0) {
   numx <- length(x)
   df <- n.obs -1 - numx
   Cr <- cov2cor(C)
    if(numx==1) { beta <- solve(C[x,x,drop=FALSE],(C[x,y,drop=FALSE])) 
                 colnames(beta) <- y
        } else {
        beta <- solve(C[x,x],(C[x,y])) } 
     if(!is.matrix(beta)) {beta <- matrix(beta,nrow=length(beta))}   #beta is a matrix of beta weights 
    if(is.character(x)) {rownames(beta) <- x}  else {rownames(beta) <- colnames(C)[x]}
      if(is.character(y)) { colnames(beta) <- y} else { colnames(beta) <- colnames(C)[y]}
         R2 <- colSums(beta * C[x,y])/diag(C[y,y,drop=FALSE])
        uniq <- 1-(1-1/diag(solve(Cr[x,x,drop=FALSE])))  #1- smc
        if(n.obs > 2) {  se <- (sqrt((1-R2)/(n.obs-1 - numx)) %*% t(sqrt(1/uniq)))  #these are the standardized se
                         se <- t( se * sqrt(diag(C[y,y,drop=FALSE])) %*% t(sqrt(1/diag(C[x,x,drop=FALSE]))) )  #But does this work in the general case?
             
                    colnames(se) <- colnames(beta) } else {se <- NA}
                    if(!any(is.na(se))) { tvalue <- beta/se
                                         prob <- 2*(1- pt(abs(tvalue),df))
                                         } else {tvalue <- prob  <- NA}
  result <- list(beta=beta,se=se, t=tvalue,prob=prob,R2=R2)
  return(result)       }  		 


boot.mediate <- function(data,x,y,m,n.iter=10,std=FALSE,use="pairwise") {
  n.obs <- nrow(data)
  numx <- length(x)
  numy <- length(y)
  numm <- length(m)
  nxy <- numx + numy 
 result <- matrix(NA,nrow=n.iter,ncol = (numx*numy))
 if((numm > 1) &  (numx > 1)) ab <- matrix(0,nrow=numx,ncol=numy)     
  for (iteration in 1:n.iter) {
  
  samp.data <- data[sample.int(n.obs,replace=TRUE),]
  C <- cov(samp.data,use=use)
  if(std) C <- cov2cor(C)
   xy <- c(x,y)     
        m.matrix <- C[c(x,m),c(x,m)]
     #	df <- n.obs - nxy - 1 
     	xy.matrix <- C[c(x,m),y,drop=FALSE]
     	
   if(numx ==1) { beta.x <- solve(C[x,x],t(C[x,y]) ) } else  {beta.x <- solve(C[x,x],C[x,y]) }  #this is the zero order beta -- the total effect
   if(numx ==0) { a <-  solve(C[x,x,drop=FALSE],t(C[x,m,drop=FALSE]) ) } else { a <-  solve(C[x,x,drop=FALSE],(C[x,m,drop=FALSE]) )} #the a paths    
  
    beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
    beta <- as.matrix(beta.xm)     #these are the individual predictors, x and then m

  b <- beta[-c(1:numx),drop=FALSE]
# if(numx ==1) {ab <- a * t(b)}  else { for (i in 1:numx) {ab[i,] <- a[i,] * b}}
 # ab <- a %*% b    #each individual path  #this probably is only correct for the numx = 1 model
#  if((numx > 1) & (numy > 1)) {for (i in 1:numx) {ab[i,] <- a[i,] * b}}
  # ab <- a * b  #we don't really need this
  indirect <-  beta.x - beta[1:numx]  #this is c' = c - ab
 # result[iteration,] <- c(indirect,ab)  #this is a list of vectors
  result[iteration,] <- c(indirect)  #this is a list of vectors
   }
 return(result)  }

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
  beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
   beta <- as.matrix(beta.xm)     #these are the individual predictors, x and then m
      
  indirect <-  beta.x - beta[1]
  a <- bM #I think
 b <- beta[-1]
 ab <- a * b
  result[iteration,] <- c(int.ind,ab)
  }
 return(result) }


"print.psych.mediate" <- function(x,digits=2,short=FALSE) {
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
    
   # if(dim(x$a)) {mv <- names(x$a)} else {mv <- colnames(x$a)
    cat("\nThe DV (Y) was ", dv,". The IV (X) was ", iv,". The mediating variable(s) = ", mv,".")
   if(!is.null(x$mod)) cat("  The moderating variable(s) = ",mod)
   for(j in 1:ndv) { 
   for(i in 1:niv) { cat("\n\nTotal Direct effect(c) of ",iv[i], " on ", dv[j]," = ",round(x$direct[i,j],digits), "  S.E. = ", round(x$total.reg$se[i,j],digits), " t direct = ",round(x$total.reg$t[i,j],digits), "  with probability = ", signif(x$total.reg$prob[i,j],digits))
    cat("\nDirect effect (c') of ",iv[i],  " on ", dv[i]," removing ", mv ," = ",round(x$indirect[i,j],digits), "  S.E. = ", round(x$cprime.reg$se[i,j],digits), " t direct = ",round(x$cprime.reg$t[i,j],digits), "  with probability = ", signif(x$cprime.reg$prob[i,j],digits))
     
   if(is.null(x$mod)) { cat("\nIndirect effect (ab) of ",iv[i], " on ", dv[j]," through " ,mv , "  = ", round(x$ab[i,j],digits),"\n")
   cat("Mean bootstrapped indirect effect = ",round(x$boot$mean[j],digits), " with standard error = ",round(x$boot$sd[j],digits), " Lower CI = ",round(x$boot$ci[1,i],digits), "   Upper CI = ", round(x$boot$ci[2,i],digits))}
     }
     cat("\nR2 of model = ", round(x$cprime.reg$R2[j],digits)) 
     
     }
     if(short) cat("\n To see the longer output, specify short = FALSE in the print statement")
   
 if(is.null(x$mod)) {

 #   cat("\nratio of indirect to total effect=  ", round(x$ratit,digits))
 #   cat("\nratio of indirect to direct effect= ", round(x$ratid,digits))
    
    
    cat("\n\n Full output  \n")
     cat("\n Total effect estimates (c) \n")
        for(j in 1:ndv) {
    dft <- round(data.frame(direct=x$total.reg$beta[,j],se = x$total.reg$se[,j],t=x$total.reg$t[,j]),digits)
    dftp <- cbind(dft,p = signif(x$total.reg$prob[,j],digits=digits+1))
    colnames(dftp) <- c(dv[j],"se","t","Prob")
    rownames(dftp) <- rownames(x$total.reg$beta)
    }
   
   
    print(dftp)
    
     cat("\nDirect effect estimates     (c') \n")
     for(j in 1:ndv) {
    if (niv==1) { dfd <- round(data.frame(direct=x$cprime.reg$beta[,j],se = x$cprime.reg$se[,j],t=x$cprime.reg$t[,j]),digits)
     dfdp <- cbind(dfd,p=signif(x$cprime.reg$prob[,j],digits=digits+1)) } else {
     dfd <- round(data.frame(direct=x$cprime.reg$beta[1:niv,j],se = x$cprime.reg$se[1:niv,j],t=x$cprime.reg$t[1:niv,j]),digits)
     dfdp <- cbind(dfd,p=signif(x$cprime.reg$prob[1:niv,j],digits=digits+1))
     }
      colnames(dfdp) <- c(dv[j],"se","t","Prob")
     }
   print(dfdp)
    
     
    cat("\n 'a'  effect estimates \n")

  if(niv==1) {
    	dfa <- round(data.frame(a = x$a.reg$beta[1,1:nmed],se = x$a.reg$se[1,1:nmed],t = x$a.reg$t[1,1:nmed]),digits)
    	dfa <- cbind(dfa,p=signif(x$a.reg$prob[1,nmed],digits=digits+1))
    	rownames(dfa) <- colnames(x$a.reg$beta)
    	colnames(dfa) <-  c(rownames(x$a.reg$beta),"se","t","Prob")
    	print(dfa)}  else {
    	
     	for (i in 1:nmed) {
     	dfa <- round(data.frame(a = x$a.reg$beta[,i],se = x$a.reg$se[,i],t = x$a.reg$t[,i]),digits)
    	dfa <- cbind(dfa,p=signif(x$a.reg$prob[,i],digits=digits+1))
     	rownames(dfa) <-rownames(x$a.reg$beta)
     	colnames(dfa) <-  c(colnames(x$a.reg$beta)[i],"se","t","Prob") 
     	print(dfa) }
     	
     	}
     	        
      cat("\n 'b'  effect estimates \n")
      for (j in 1:ndv) {
      if(niv==1) {
     dfb <- round(data.frame(direct=x$b.reg$beta[-(1:niv),j],se = x$b.reg$se[-(1:niv),j],t=x$b.reg$t[-(1:niv),j]),digits)
     dfb <- cbind(dfb,p=signif(x$b.reg$prob[-(1:niv),j],digits=digits+1))} else {
      dfb <- round(data.frame(direct=x$b.reg$beta[-(1:niv),j],se = x$b.reg$se[-(1:niv),j],t=x$b.reg$t[-(1:niv),j]),digits)
     dfb <- cbind(dfb,p=signif(x$b.reg$prob[-(1:niv),j],digits=digits+1))}
     rownames(dfb) <- rownames(x$b.reg$beta)[-(1:niv)]
     colnames(dfb) <-  c(dv[j],"se","t", "Prob")
      print(dfb)
      }
 
      cat("\n 'ab'  effect estimates \n")
      #does not work yet for ndv > 1  for boot
 for (j in 1:ndv) {
      dfab  <- round(data.frame(indirect=x$ab[,j],boot = x$boot$mean[,j],sd=x$boot$sd,lower = x$boot$ci[1,],upper= x$boot$ci[2,]),digits)
      rownames(dfab) <- rownames(x$ab)
      colnames(dfab)[1] <- dv[j]
      print(round(dfab,digits))
      }

    }  else {
    cat("\n\nEffect of interaction of ",iv[1], " with ", iv[2] , "  = ", round(x$direct[3],digits),"  S.E. = ", round(x$direct.reg$se[3,1],digits), " t direct = ",round(x$direct.reg$t[3,1],digits), "  with probability = ", signif(x$direct.reg$prob[3,1],digits))
    cat("\nIndirect effect due to interaction  of ",iv[1], " with ", iv[2] , "  = ", round(x$indirect,digits))
    cat("\nMean bootstrapped indirect interaction effect = ",round(x$boot$mean[1],digits), " with standard error = ",round(x$boot$sd[1],digits), " Lower CI = ",round(x$boot$ci.ab[1],digits), "   Upper CI = ", round(x$boot$ci.ab[2,i],digits))
    cat("\nSummary of a, b, and ab estimates and ab confidence intervals\n")
    }     
}

   
   
"mediate.diagram" <- function(medi,digits=2,ylim=c(3,7),xlim=c(-1,10),show.c=TRUE, main="Mediation model",...) { 

    
     dv <- medi$var.names[["DV"]]
   # iv <- medi$var.names[["IV"]]
     iv <- as.matrix(rownames(medi$direct))
    mv <- medi$var.names[["med"]]
    mod <- medi$var.names[["mod"]]
    numx <- length(medi$var.names[["IV"]])
    numy <- length(dv)
    direct <- round(medi$direct,digits)
    C <- round(medi$C[c(iv,mv,dv),c(iv,mv,dv)],digits)
    

miny <- 5 - max(length(iv)/2,length(mv),2) - .5
maxy <- 5 + max(length(iv)/2,length(mv),2) + .5



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
 x[[i]]  <- dia.rect(1,viv[i],iv[i])}
 
 vdv <- 1:numy
 y <- list()
 for (i in 1:numy) {
 if((numy %% 2)== 0) {
vdv[i] <- switch(i,6,4,7,3,8,2,9,1,10)  } else { vdv[i] <- switch(i,5,7,3,6,4,8,2,9)}
 y[[i]] <- dia.rect(9,vdv[i],dv[i]) 
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
     			 dia.arrow(x[[i]]$right,y[[j]]$left,labels=paste("c = ",direct[i,j]),pos=0,...)}
	}  
} else {
if(n.mediate==1) a <- t(a)
for (mediate in 1:n.mediate) {
	m[[mediate]] <- dia.rect(5,v.loc[mediate],mv[mediate] ) 
		for(j in 1:numy) {
			for(i in 1: numx) {dia.arrow(x[[i]]$right,m[[mediate]]$left,labels=a[i,mediate],adj=adj,...) #a term
     			if(show.c) {dia.arrow(x[[i]]$right,y[[j]]$left,labels=paste("c = ",c[i,j]),pos=3,...)}
     			 dia.arrow(x[[i]]$right,y[[j]]$left,labels=paste("c' = ",cprime[i,j]),pos=1,...)}
      			dia.arrow(m[[mediate]]$right,y[[j]]$left,labels=b[mediate,j],...)  #     
			}
		} 
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
     		 
     	
             
 