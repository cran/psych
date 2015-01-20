"mediate" <-
function(y,x,m,data, mod=NULL,n.obs=NULL,use="pairwise",n.iter=5000,alpha=.05,std=FALSE)  {
 
   cl <- match.call()
  
  if(!is.matrix(data)) data <- as.matrix(data)
  if((dim(data)[1]!=dim(data)[2]))  {n.obs=dim(data)[1]   #this does not take into account missing data
                    if(!is.null(mod)) data <- scale(data,scale=FALSE)
                    C <- cov(data,use=use)
                    raw <- TRUE
                     }  else {
                    raw <- FALSE
                    C <-data
                    n.obs=0

                  }
if(std) { C <- cov2cor(C)}  #use correlations rather than covariances
   #convert names to locations 
   #first, see if they are in formula mode  
   if(class(y) == "formula") { yn <- all.vars(y)
            y <- yn[1]
            x <- yn[2]
            m <- yn[3:length(yn)]
            }  
                   
   if(is.character(x)) x <- which(colnames(data) == x)
   if(is.character(y)) y <- which(colnames(data) %in% y) 
   if(is.character(m)) m <- which( colnames(data) %in%  m)
  
  if(!is.null(mod)) {ivXm <- matrix(data[,x] * data[,mod],ncol=1)
              colnames(ivXm) <- paste(abbreviate(colnames(data)[x]),"X",abbreviate(colnames(data)[mod]),sep="")
              data <- cbind(data,ivXm)
              if(!(mod %in% m)) {m <- c(m,mod,ncol(data))} else {m <- c(m,ncol(data))}
              C <- cov(data,use=use)
              if(std) { C <- cov2cor(C)}
              }
  #now, we are ready to process the data
  #some of this is taken from set.cor for the more general case
     # nm <- dim(data)[1]
        xy <- c(x,y)
        numx <- 1
     	numy <- 1
     	numm <- length(m)
        nxy <- numx + numy 
        m.matrix <- C[c(x,m),c(x,m)]
     	df <- n.obs - nxy - 1 

     	xy.matrix <- C[c(x,m),y,drop=FALSE]
     	
    beta.x <- C[x,y]/C[x,x]   #this is the zero order beta -- the total effect
    a <- C[x,m]/C[x,x]  #this is the direct path from X to M  
    R2.x <- beta.x *   C[x,y] /(C[y,y])  #R2 for the x predictor by itself 
     
     	 uniq.x <- C[x,x]*(1- R2.x)
     	 se.bx <- sqrt(((1-R2.x) * (C[y,y])/ (C[x,x] * (n.obs-2) )  ))
     	 
     	tt <- beta.x/se.bx       #total direct t
     	probt <- 2*(1- pt(abs(tt),n.obs-2))

    beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
    beta <- as.matrix(beta.xm)     #these are the individual predictors, x and then m and then mod
    b <- beta[-1]  #these are the b paths 
    indirect <-  beta.x - beta[1]  #this is c' = c - ab
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
          mean.boot <- sd.boot <- ci.quant <- boot <-  NA } else {
            if(is.null(mod)) {boot <- boot.mediate(data,x,y,m,n.iter=n.iter,std=std) }  else {boot <- boot.moderate( data,x,y,m,mod,n.iter=n.iter,std=std)}
            mean.boot <- mean(boot)
            sd.boot <- sd(boot) 
            ci.quant <- quantile(boot,c(alpha/2,1-alpha/2),na.rm=TRUE)   
              }    
    result <- list(total=beta.x,se.bx=se.bx,tt=tt,probt=probt,direct=beta.xm, indirect = indirect,se.beta=se,t=tvalue,prob=prob,R2.x=R2.x,R2=R2,
       ratit=ratit,ratid=ratid, mean.boot=mean.boot,sd.boot=sd.boot,ci.quant=ci.quant,boot=boot, mod=mod, int.ind=int.ind,a=a,b=b,b.int=b.int,mod=mod,names=colnames(data),  Call=cl)
    class(result) <- c("psych","mediate")
    if(is.null(mod)) { mediate.diagram(result) } else {moderate.diagram(result)}
    
    return(result)	
     		 }
     		 

boot.mediate <- function(data,x,y,m,n.iter=10,std=FALSE) {
  result <- matrix(NA,n.iter)
  n.obs <- nrow(data)
  for (iteration in 1:n.iter) {
  samp.data <- data[sample.int(n.obs,replace=TRUE),]
  C <- cov(samp.data)
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
         
   R2.x <- beta.x *   C[x,y] /(C[y,y])  #R2 for the x predictor by itself 
     
     	 uniq.x <- C[x,x]*(1- R2.x)
     	 se.bx <- sqrt(((1-R2.x) * (C[y,y])/ (C[x,x] * (n.obs-2) )  ))
     	tt <- beta.x/se.bx       #total direct t
     	probt <- 2*(1- pt(abs(tt),n.obs-2))

    beta.xm <- solve(m.matrix,xy.matrix)       #solve the equation bY~aX
   beta <- as.matrix(beta.xm)     #these are the individual predictors, x and then m
      
  indirect <-  beta.x - beta[1]
  result[iteration] <- indirect
  }
 return(result) 
}

boot.moderate <- function(data,x,y,m,mod,n.iter=10,std=FALSE) {
  result <- matrix(NA,n.iter)
  n.obs <- nrow(data)
  for (iteration in 1:n.iter) {
  samp.data <- data[sample.int(n.obs,replace=TRUE),]
  samp.data <- scale(samp.data,scale=FALSE)
  ivXm <- samp.data[,x] * samp.data[,mod]
              samp.data <- cbind(samp.data,ivXm)
  C <- cov(samp.data)
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
  
  result[iteration] <- int.ind
  }
 return(result) 
}


"print.psych.mediate" <- function(x,digits=2) {
 value <- class(x)[2] 

 switch(value, 
mediate = { cat("\nMediation analysis \n")
 cat("Call: ")
    print(x$Call)
    dv <- colnames(x$direct)
    iv <- rownames(x$direct)
    cat("\nThe DV (Y) was ", dv,". The IV (X) was ", iv[1],". The mediating variable(s) = ", iv[-1],".")
    if(!is.null(x$mod)) cat("  The moderating variable(s) = ", x$names[x$mod])
    
    cat("\nTotal Direct effect(c) of ",iv[1], " on ", dv," = ",round(x$total,digits), "  S.E. = ", round(x$se.bx,digits), " t direct = ",round(x$tt,digits), "  with probability = ", signif(x$probt,digits))
    cat("\nDirect effect of ",iv[1],  " on ", dv," removing ", iv[-1] ," = ",round(x$direct[1],digits), "  S.E. = ", round(x$se.beta[1],digits), " t direct = ",round(x$t[1],digits), "  with probability = ", signif(x$prob[1],digits))
    
   if(is.null(x$mod)) { cat("\nIndirect effect (c') of ",iv[1], " on ", dv," through " ,iv[-1] , "  = ", round(x$indirect,digits),"\n")
     
    cat("\nMean bootstrapped indirect effect = ",round(x$mean.boot,digits), " with standard error = ",round(x$sd.boot,digits), " Lower CI = ",round(x$ci.quant[1],digits), "   Upper CI = ", round(x$ci.quant[2],digits))
    
    
    cat("\nratio of indirect to total effect=  ", round(x$ratit,digits))
    cat("\nratio of indirect to direct effect= ", round(x$ratid,digits))
    }  else {
    cat("\nEffect of interaction of ",iv[1], " with ", iv[2] , "  = ", round(x$direct[3],digits),"  S.E. = ", round(x$se.beta[3,1],digits), " t direct = ",round(x$t[3,1],digits), "  with probability = ", signif(x$prob[3,1],digits))
    cat("\nIndirect effect due to interaction  of ",iv[1], " with ", iv[2] , "  = ", round(x$int.ind,digits))
    cat("\nMean bootstrapped indirect interaction effect = ",round(x$mean.boot,digits), " with standard error = ",round(x$sd.boot,digits), " Lower CI = ",round(x$ci.quant[1],digits), "   Upper CI = ", round(x$ci.quant[2],digits))
    }

    cat("\nR2 of model = ", round(x$R2,digits))
})
}

   
   
"mediate.diagram" <- function(medi,digits=2,...) {
  
xlim=c(0,10)
ylim=c(0,10)
plot(NA,xlim=xlim,ylim=ylim,main="Mediation model",axes=FALSE,xlab="",ylab="")
var.names <- c(rownames(medi$direct)[1],colnames(medi$direct),rownames(medi$direct)[-1])
n.mediate <- length(var.names) - 2

m <- list()
c <- round(medi$total,digits)
 a <- round(medi$a,digits)
 b <- round(medi$b,digits)
 cprime <- round(medi$direct[1],digits)
x <- dia.rect(1,5,var.names[1])
y <- dia.rect(9,5,var.names[2]) 
for (mediate in 1:n.mediate) {
oddeven <- ((mediate %% 2)-.5)*2
m[[mediate]] <- dia.rect(5,5+(floor((mediate+1)/2)*2 )*oddeven,var.names[mediate+2] ) 
#dia.arrow(x$right,m[[mediate]]$left,labels=paste("a = ",a[mediate]),...) 
#dia.arrow(m[[mediate]]$right,y$left,labels=paste("b = ",b[mediate]),...)
dia.arrow(x$right,m[[mediate]]$left,labels=a[mediate],...) 
dia.arrow(m[[mediate]]$right,y$left,labels=b[mediate],...)

#m <- dia.rect(5,8,var.names[3])
dia.arrow(x$right,y$left,labels=paste("c = ",c),pos=3,...)
dia.arrow(x$right,y$left,labels=paste("c' = ",cprime),pos=1,...)

}
}	


"moderate.diagram" <- function(medi,digits=2,main="Moderation",...) {
  
xlim=c(0,10)
ylim=c(0,10)
plot(NA,xlim=xlim,ylim=ylim,main="Mediation model",axes=FALSE,xlab="",ylab="")
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
     		 
     	
             
 