"mat.regress" <-
function(y,x,data,z=NULL,n.obs=NULL,use="pairwise",square=FALSE)  {
 #a function to extract subsets of variables (a and b) from a correlation matrix m or data set m
  #and find the multiple correlation beta weights + R2 of the a set predicting the b set
  #seriously rewritten, March 24, 2009 to make much simpler
  #minor additons, October, 20, 2009 to allow for print and summary function
  #major addition in April, 2011 to allow for set correlation
  
  message("mat.regress has been replaced by setCor, please change your call") 
  setCor(y,x,data,z=NULL,n.obs=NULL,use="pairwise",square=FALSE)} 
 #  
#    cl <- match.call()
#   if(!is.matrix(data)) data <- as.matrix(data)
#   if((dim(data)[1]!=dim(data)[2]) |square) {n.obs=dim(data)[1]
#                     C <- cov(data,use=use)
#                     m <- cov2cor(C)
#                      raw <- TRUE}  else {
#                     raw <- FALSE  
#                     C <-data
#                     m <- cov2cor(C)}
#                    
#  
#         nm <- dim(data)[1]
#         xy <- c(x,y)
#         numx <- length(x)
#      	numy <- length(y)
#         nxy <- numx+numy
#         m.matrix <- m[c(x,y),c(x,y)]
#      	a.matrix <- m[x,x]
#      	ac.matrix <- C[x,x]
#      	b.matrix <- m[x,y]
#      	bc.matrix <- C[x,y]
#      	y.matrix <- m[y,y]
#      	if(numx == 1 ) {beta <- matrix(b.matrix,nrow=1)
#      	                } else   #this is the case of a single x 
#        { beta <- solve(a.matrix,b.matrix)       #solve the equation bY~aX
#         beta <- as.matrix(beta) }
#      	if (numy >1 ) { if(is.null(rownames(beta))) {rownames(beta) <- colnames(m)[x]}
#      	                if(is.null(colnames(beta))) {colnames(beta) <- colnames(m)[y]}
#      	 
#      	 R2 <- colSums(beta * b.matrix) } else { colnames(beta) <- colnames(data)[1]
#      		 R2 <- sum(beta * b.matrix)
#      	 	names(beta) <- colnames(data)[x]
#      		 names(R2) <- colnames(data)[y]}
#      	if(numy < 2) {Rset <- 1 - det(m.matrix)/(det(a.matrix) )
#      	            } else {if (numx < 2) {Rset <- 1 - det(m.matrix)/(det(y.matrix) )
#      	            } else {Rset <- 1 - det(m.matrix)/(det(a.matrix) * det(y.matrix))}
#      	            }
#      	if(!is.null(n.obs)) {k<- length(x)
#      	                     
#      	                     uniq <- (1-smc(a.matrix))
#      	                     se.beta <- list() 
#      	                     for (i in 1:length(y)) {
#      	                     df <- n.obs-k-1
#      	                     se.beta[[i]] <- (sqrt((1-R2[i])/(df))*sqrt(1/uniq))}
#      	                     se <- matrix(unlist(se.beta),ncol=length(y))
#      	                     colnames(se) <- colnames(beta)
#      	                     rownames(se) <- rownames(beta)
#      	                     tvalue <- beta/se
#      	                     se <- t(t(se) * sqrt(diag(C)[y]))/sqrt(diag(ac.matrix))
#      	                     
#      	                prob <- 2*(1- pt(abs(tvalue),df))
#      	                     SE2 <- 4*R2*(1-R2)^2*(df^2)/((n.obs^2-1)*(n.obs+3))
#      	                     SE =sqrt(SE2)
#      	                     F <- R2*df/(k*(1-R2))
#      	                     pF <- 1 - pf(F,k,df)
#      	                     shrunkenR2 <- 1-(1-R2)*(n.obs-1)/df 
#      	                     
#      	               #find the shrunken R2 for set cor  (taken from CCAW p 615)
#      	                     u <- numx * numy
#      	                     m1 <- n.obs - numx * numy - (numx + numy +3)/2 
#      	                     s <- sqrt((numx ^2 * numy^2 -4)/(numx^2 + numy^2))
#      	                     v <- m1 * s + 1 - u/2
#      	                     R2set.shrunk <- 1 - (1-Rset) * ((v+u)/v)^s
#      	                              }
#      	
#      	if(numx == 1) {beta <-  beta * sqrt(diag(C)[y])
#      	   } else {beta <-  t(t(beta) * sqrt(diag(C)[y]))/sqrt(diag(ac.matrix))} #this puts the betas into the raw units
#         
#      	if(is.null(n.obs)) {mat.regress <- list(beta=beta,R=sqrt(R2),R2=R2,Rset=Rset,raw=raw,Call = cl)} else {
#      	              mat.regress <- list(beta=beta,se=se,t=tvalue,Probability = prob,R=sqrt(R2),R2=R2,shrunkenR2 = shrunkenR2,seR2 = SE,F=F,probF=pF,df=c(k,df),Rset=Rset,Rset.shrunk=R2set.shrunk,raw=raw,Call = cl)}
#      	class(mat.regress) <- c("psych","set.cor")
#      	return(mat.regress)
#      	}
#modified July 12,2007 to allow for NA in the overall matrix
#modified July 9, 2008 to give statistical tests
#modified yet again August 15 , 2008 to convert covariances to correlations
#modified January 3, 2011 to work in the case of a single predictor 
#modified April 25, 2011 to add the set correlation (from Cohen)
