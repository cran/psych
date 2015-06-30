"set.cor" <-
function(y,x,data,z=NULL,n.obs=NULL,use="pairwise",std=TRUE,square=FALSE,main="Regression Models")  {
setCor(y=y,x=x,data=data,z=z,n.obs=n.obs,use=use,std=std,square=square,main=main)}


"setCor" <-
function(y,x,data,z=NULL,n.obs=NULL,use="pairwise",std=TRUE,square=FALSE,main="Regression Models")  {

 #a function to extract subsets of variables (a and b) from a correlation matrix m or data set m
  #and find the multiple correlation beta weights + R2 of the a set predicting the b set
  #seriously rewritten, March 24, 2009 to make much simpler
  #minor additons, October, 20, 2009 to allow for print and summary function
  #major addition in April, 2011 to allow for set correlation
  #added calculation of residual matrix December 30, 2011
  #added option to allow square data matrices
  #modified December, 2014  to allow for covariances as well as to fix a bug with single x variable
  #modified April, 2015 to handle data with non-numeric values in the data, which are not part of the analysis
  
   cl <- match.call()
    if(is.numeric(y )) y <- colnames(data)[y]
    if(is.numeric(x )) x <- colnames(data)[x]
    if(is.numeric(z )) z <- colnames(data)[z]
  
  
  if((dim(data)[1]!=dim(data)[2]) | square)  {n.obs=dim(data)[1]   #this does not take into account missing data
                  if(!is.null(z))   {data <- data[,c(y,x,z)]} else {data <- data[,c(y,x)]}
                   if(!is.matrix(data)) data <- as.matrix(data)  
                   if(!is.numeric(data)) stop("The data must be numeric to proceed")
                    C <- cov(data,use=use)
                    if(std) {m <- cov2cor(C)} else {m <- C}
                    raw <- TRUE}  else {
                    raw <- FALSE
                    if(!is.matrix(data)) data <- as.matrix(data)  
                    C <- data
                    if(std) {m <- cov2cor(C)} else {m <- C}}
   #convert names to locations                 
   #if(is.character(x)) x <- which(colnames(data) == x)
   #if(is.character(y)) y <- which(colnames(data) == y) 
   #if(!is.null(z) && is.character(z)) z <-  which(colnames(data) == z)               
 
        nm <- dim(data)[1]
        xy <- c(x,y)
        numx <- length(x)
     	numy <- length(y)
     	numz <- 0
        nxy <- numx+numy
        m.matrix <- m[c(x,y),c(x,y)]
     	x.matrix <- m[x,x,drop=FALSE]
     	xc.matrix <- C[x,x,drop=FALSE]
     	xy.matrix <- m[x,y,drop=FALSE]
     	xyc.matrix <- C[x,y,drop=FALSE]
     	y.matrix <- m[y,y,drop=FALSE]
     	

        	
     	if(!is.null(z)){numz <- length(z)      #partial out the z variables
     	                zm <- m[z,z,drop=FALSE]
     	                za <- m[x,z,drop=FALSE]
     	                 zb <- m[y,z,drop=FALSE]
     	                 x.matrix <- x.matrix - za %*% solve(zm) %*% t(za)
     	                 y.matrix <- y.matrix - zb %*% solve(zm) %*% t(zb)
     	                xy.matrix <- xy.matrix - za  %*% solve(zm) %*% t(zb)
     	                m.matrix <- cbind(rbind(y.matrix,xy.matrix),rbind(t(xy.matrix),x.matrix))
     	                
     	                 }
     	if(numx == 1 ) {beta <- matrix(xy.matrix,nrow=1)/x.matrix[1,1]
     	                } else   #this is the case of a single x 
      				 { beta <- solve(x.matrix,xy.matrix)       #solve the equation bY~aX
       					 beta <- as.matrix(beta) 
       					 }


       	yhat <- t(xy.matrix) %*% solve(x.matrix) %*% (xy.matrix)
       	resid <- y.matrix - yhat
    	if (numy > 1 ) { 
    	               if(is.null(rownames(beta))) {rownames(beta) <- x}
    	                if(is.null(colnames(beta))) {colnames(beta) <- y}
     	 
     	 R2 <- colSums(beta * xy.matrix)/diag(y.matrix) } else {  
     	 colnames(beta) <- y
     		
     		 R2 <- sum(beta * xy.matrix)/y.matrix
     		 R2 <- matrix(R2)
    		      	 	 rownames(beta) <- x
    		 rownames(R2) <- colnames(R2) <- y
     		 }
     	

     	
     	#now find the unit weighted correlations  
     	#reverse items in X and Y so that they are all positive signed
     	 px <- principal(x.matrix)
            keys.x <- diag(as.vector(1- 2* (px$loadings < 0 )) )
         py <- principal(y.matrix)
            keys.y <- diag(as.vector(1- 2* (py$loadings < 0 ) ))

        Vx <- sum( keys.x %*% x.matrix %*% t(keys.x))
        Vy <- sum( keys.y %*% y.matrix %*% t(keys.y))
     		 
     	ruw <- colSums(abs(xy.matrix))/sqrt(Vx)	 
     	Ruw <- sum(diag(keys.x) %*% xy.matrix %*% t(keys.y))/sqrt(Vx * Vy)
     	
   
     	if(numy < 2) {Rset <- 1 - det(m.matrix)/(det(x.matrix) )
     	             Myx <- solve(x.matrix) %*% xy.matrix  %*% t(xy.matrix)
     	             cc2 <- cc <- T <- NULL} else {if (numx < 2) {Rset <- 1 - det(m.matrix)/(det(y.matrix) )
     	            Myx <-  xy.matrix %*% solve(y.matrix) %*% t(xy.matrix)
     	            cc2 <- cc <- T <- NULL} else {Rset <- 1 - det(m.matrix)/(det(x.matrix) * det(y.matrix))
     	            if(numy > numx) {
     	            Myx <- solve(x.matrix) %*% xy.matrix %*% solve(y.matrix) %*% t(xy.matrix)} else { Myx <- solve(y.matrix) %*% t(xy.matrix )%*% solve(x.matrix) %*% (xy.matrix)}
     	           }
     	            cc2 <- eigen(Myx)$values
     	            cc <- sqrt(cc2)
     	            T <- sum(cc2)/length(cc2)             
     	            }
     	       
     	if(!is.null(n.obs)) {k<- length(x)
     	                     
     	                    # uniq <- (1-smc(x.matrix,covar=!std))
     	                     uniq <- (1-smc(x.matrix))
     	                     se.beta <- list() 
     	                     for (i in 1:length(y)) {
     	                     df <- n.obs-k-1
     	                     se.beta[[i]] <- (sqrt((1-R2[i])/(df))*sqrt(1/uniq))}    #why is uniq not indexed if R2 is?
     	                     se <- matrix(unlist(se.beta),ncol=length(y))
     	                     colnames(se) <- colnames(beta)
     	                     rownames(se) <- rownames(beta)
     	                     
     	                     se <- t(t(se) * sqrt(diag(C)[y]))/sqrt(diag(xc.matrix))
     	                     tvalue <- beta/se
     	                prob <- 2*(1- pt(abs(tvalue),df))
     	                     SE2 <- 4*R2*(1-R2)^2*(df^2)/((n.obs^2-1)*(n.obs+3))
     	                     SE =sqrt(SE2)
     	                     F <- R2*df/(k*(1-R2))
     	                     pF <- 1 - pf(F,k,df)
     	                     shrunkenR2 <- 1-(1-R2)*(n.obs-1)/df 
     	                     
     	               #find the shrunken R2 for set cor  (taken from CCAW p 615)
     	                     u <- numx * numy
     	                     m1 <- n.obs - max(numy ,(numx+numz)) - (numx + numy +3)/2 
     	                    
     	                     s <- sqrt((numx ^2 * numy^2 -4)/(numx^2 + numy^2-5)) 
     	                     if(numx*numy ==4) s <- 1
     	                   
     	                     v <- m1 * s + 1 - u/2
     	                     R2set.shrunk <- 1 - (1-Rset) * ((v+u)/v)^s
     	                     L <- 1-Rset
     	                     L1s <- L^(-1/s)
     	                     Rset.F <- (L1s-1)*(v/u)
     	                     df.m <- n.obs -  max(numy ,(numx+numz))  -(numx + numy +3)/2
     	                      s1 <- sqrt((numx ^2 * numy^2 -4)/(numx^2 + numy^2-5)) #see cohen p 321
     	                     if(numx^2*numy^2 < 5) s1 <- 1
     	                     df.v <- df.m * s1 + 1 - numx * numy/2 #which is just v
     	                    # df.v <- (u+v)  #adjusted for bias to match the CCAW results
     	                    #Rset.F <- Rset.F * (u+v)/v 
     	                   
     	                    Chisq <- -(n.obs - 1 -(numx + numy +1)/2)*log((1-cc2))
     	                   
     	                              }
     	
     #	if(numx == 1) {beta <-  beta * sqrt(diag(C)[y])
     #	   } else {beta <-  t(t(beta) * sqrt(diag(C)[y]))/sqrt(diag(xc.matrix))} #this puts the betas into the raw units
        
     	if(is.null(n.obs)) {set.cor <- list(beta=beta,R=sqrt(R2),R2=R2,Rset=Rset,T=T,cancor = cc, cancor2=cc2,raw=raw,residual=resid,ruw=ruw,Ruw=Ruw,x.matrix=x.matrix,y.matrix=y.matrix,Call = cl)} else {
     	              set.cor <- list(beta=beta,se=se,t=tvalue,Probability = prob,R=sqrt(R2),R2=R2,shrunkenR2 = shrunkenR2,seR2 = SE,F=F,probF=pF,df=c(k,df),Rset=Rset,Rset.shrunk=R2set.shrunk,Rset.F=Rset.F,Rsetu=u,Rsetv=df.v,T=T,cancor=cc,cancor2 = cc2,Chisq = Chisq,raw=raw,residual=resid,ruw=ruw,Ruw=Ruw,x.matrix=x.matrix,y.matrix=y.matrix,Call = cl)}
     	class(set.cor) <- c("psych","setCor")
     	setCor.diagram(set.cor,main=main)
     	return(set.cor)
     	
     	}
#modified July 12,2007 to allow for NA in the overall matrix
#modified July 9, 2008 to give statistical tests
#modified yet again August 15 , 2008 to convert covariances to correlations
#modified January 3, 2011 to work in the case of a single predictor 
#modified April 25, 2011 to add the set correlation (from Cohen)
#modified April 21, 2014 to allow for mixed names and locations in call
#modified February 19, 2015 to just find the covariances of the data that are used in the regression
#this gets around the problem that some users have large data sets, but only want a few variables in the regression


setCor.diagram <- function(sc,main="Regression model",digits=2,show=TRUE,...) { 

beta <- round(sc$beta,digits)
x.matrix <- round(sc$x.matrix,digits)
y.matrix <- round(sc$y.matrix,digits)
x.names <- rownames(sc$beta)
y.names <- colnames(sc$beta)
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
plot(NA,xlim=xlim,ylim=ylim,main=main,axes=FALSE,xlab="",ylab="")
for(i in 1:nx) {x[[i]] <- dia.rect(3,top-i*x.scale,x.names[i]) }
 
for (j in 1:ny) {y[[j]] <- dia.rect(7,top-j*y.scale,y.names[j]) }
for(i in 1:nx) {
  for (j in 1:ny) {
   dia.arrow(x[[i]]$right,y[[j]]$left,labels = beta[i,j],adj=4-j)
   }
} 
if(nx >1) {
  for (i in 2:nx) {
  for (k in 1:(i-1)) {dia.curved.arrow(x[[i]]$left,x[[k]]$left,x.matrix[i,k],scale=-(abs(i-k)),both=TRUE)} 
  } }
  
  if(ny>1) {for (i in 2:ny) {
  for (k in 1:(i-1)) {dia.curved.arrow(y[[i]]$right,y[[k]]$right,y.matrix[i,k],scale=(abs(i-k)))} 
  }}
  for(i in 1:ny) {dia.self(y[[i]],side=3,scale=.2) }
 if(show) {text((10-nx/3)/2,0,paste("unweighted matrix correlation  = ",round(sc$Ruw,digits)))}
}
		 
     		 
     	