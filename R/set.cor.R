"set.cor" <-
function(y,x,data,z=NULL,n.obs=NULL,use="pairwise",square=FALSE)  {
 #a function to extract subsets of variables (a and b) from a correlation matrix m or data set m
  #and find the multiple correlation beta weights + R2 of the a set predicting the b set
  #seriously rewritten, March 24, 2009 to make much simpler
  #minor additons, October, 20, 2009 to allow for print and summary function
  #major addition in April, 2011 to allow for set correlation
  #added calculation of residual matrix December 30, 2011
  #added option to allow square data matrices
   cl <- match.call()
   cor <- TRUE  #fix this to allow this to a parameter
  if(!is.matrix(data)) data <- as.matrix(data)
  if((dim(data)[1]!=dim(data)[2]) | square)  {n.obs=dim(data)[1]   #this does not take into account missing data
                    C <- cov(data,use=use)
                    if(cor) {m <- cov2cor(C)} else {m <- C}
                    raw <- TRUE}  else {
                    raw <- FALSE
                    C <-data
                    if(cor) {m <- cov2cor(C)} else {m <- C}}
   #convert names to locations                 
   if(is.character(x)) x <- which(colnames(data) == x)
   if(is.character(y)) y <- which(colnames(data) == y) 
   if(!is.null(z) && is.character(z)) z <-  which(colnames(data) == z)               
 
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
     	if(!is.null(z)){numz <- length(z)
     	                zm <- m[z,z,drop=FALSE]
     	                 za <- m[x,z,drop=FALSE]
     	                 zb <- m[y,z,drop=FALSE]
     	                 x.matrix <- x.matrix - za %*% solve(zm) %*% t(za)
     	                 y.matrix <- y.matrix - zb %*% solve(zm) %*% t(zb)
     	                xy.matrix <- xy.matrix - za  %*% solve(zm) %*% t(zb)
     	                m.matrix <- cbind(rbind(y.matrix,xy.matrix),rbind(t(xy.matrix),x.matrix))
     	                
     	                 }
     	if(numx == 1 ) {beta <- matrix(xy.matrix,nrow=1)
     	                } else   #this is the case of a single x 
      				 { beta <- solve(x.matrix,xy.matrix)       #solve the equation bY~aX
       					 beta <- as.matrix(beta) 
       					 }
       				 
       	yhat <- t(xy.matrix) %*% solve(x.matrix) %*% (xy.matrix)
       	resid <- y.matrix - yhat
     	if (numy >1 ) { if(is.null(rownames(beta))) {rownames(beta) <- colnames(m)[x]}
     	                if(is.null(colnames(beta))) {colnames(beta) <- colnames(m)[y]}
     	 
     	 R2 <- colSums(beta * xy.matrix) } else { colnames(beta) <- colnames(data)[y]
     		 R2 <- sum(beta * xy.matrix)
     	 	 rownames(beta) <- colnames(data)[x]
     		 names(R2) <- colnames(data)[y]}
     		 
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
     	                     
     	                     uniq <- (1-smc(x.matrix))
     	                     se.beta <- list() 
     	                     for (i in 1:length(y)) {
     	                     df <- n.obs-k-1
     	                     se.beta[[i]] <- (sqrt((1-R2[i])/(df))*sqrt(1/uniq))}
     	                     se <- matrix(unlist(se.beta),ncol=length(y))
     	                     colnames(se) <- colnames(beta)
     	                     rownames(se) <- rownames(beta)
     	                     tvalue <- beta/se
     	                     se <- t(t(se) * sqrt(diag(C)[y]))/sqrt(diag(xc.matrix))
     	                     
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
     	
     	if(numx == 1) {beta <-  beta * sqrt(diag(C)[y])
     	   } else {beta <-  t(t(beta) * sqrt(diag(C)[y]))/sqrt(diag(xc.matrix))} #this puts the betas into the raw units
        
     	if(is.null(n.obs)) {set.cor <- list(beta=beta,R=sqrt(R2),R2=R2,Rset=Rset,T=T,cancor = cc, cancor2=cc2,raw=raw,residual=resid,ruw=ruw,Ruw=Ruw,Call = cl)} else {
     	              set.cor <- list(beta=beta,se=se,t=tvalue,Probability = prob,R=sqrt(R2),R2=R2,shrunkenR2 = shrunkenR2,seR2 = SE,F=F,probF=pF,df=c(k,df),Rset=Rset,Rset.shrunk=R2set.shrunk,Rset.F=Rset.F,Rsetu=u,Rsetv=df.v,T=T,cancor=cc,cancor2 = cc2,Chisq = Chisq,raw=raw,residual=resid,ruw=ruw,Ruw=Ruw,Call = cl)}
     	class(set.cor) <- c("psych","set.cor")
     	return(set.cor)
     	}
#modified July 12,2007 to allow for NA in the overall matrix
#modified July 9, 2008 to give statistical tests
#modified yet again August 15 , 2008 to convert covariances to correlations
#modified January 3, 2011 to work in the case of a single predictor 
#modified April 25, 2011 to add the set correlation (from Cohen)
#modified April 21, 2014 to allow for mixed names and locations in call