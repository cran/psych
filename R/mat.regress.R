"mat.regress" <-
function(m,x,y,n.obs=NULL,digits=2)  {
 #a function to extract subsets of variables (a and b) from a correlation matrix m or data set m
  #and find the multiple correlation beta weights + R2 of the a set predicting the b set
  #seriously rewritten, March 24, 2009 to make much simpler
  
  if(!is.matrix(m)) m <- as.matrix(m)
  if(dim(m)[1]!=dim(m)[2]) {n.obs=dim(m)[1]
                    C <- cov(m,use="pairwise")
                    m <- cov2cor(C)}  else {
                    C <- m
                    m <- cov2cor(m)}
                   
 
        nm <- dim(m)[1]
        xy <- c(x,y)
        numx <- length(x)
     	numy <- length(y)
        nxy <- numx+numy
     	a.matrix <- m[x,x]
     	ac.matrix <- C[x,x]
     	b.matrix <- m[x,y]
     	bc.matrix <- C[x,y]
     	
     beta<- solve(a.matrix,b.matrix)       #solve the equation bY~aX
     	if (length(y) >1 ) { rownames(beta) <- colnames(m)[x]
     	 colnames(beta) <- colnames(m)[y]
     	 
     	R2 <- colSums(beta * b.matrix) } else { R2 <- sum(beta * b.matrix)
     	 names(beta) <- colnames(m)[x]
     	 names(R2) <- colnames(m)[y]}
     	if(!is.null(n.obs)) {k<- length(x)
     	                     
     	                     uniq <- (1-smc(a.matrix))
     	                     se.beta <- list() 
     	                     for (i in 1:length(y)) {
     	                     df <- n.obs-k-1
     	                     se.beta[[i]] <- (sqrt((1-R2[i])/(df))*sqrt(1/uniq))}
     	                     se <- matrix(unlist(se.beta),ncol=length(y))
     	                     colnames(se) <- colnames(beta)
     	                     rownames(se) <- rownames(beta)
     	                     tvalue <- beta/se
     	                     se <- t(t(se) * sqrt(diag(C)[y]))/sqrt(diag(ac.matrix))
     	                     
     	                prob <- 2*(1- pt(abs(tvalue),df))
     	                     SE2 <- 4*R2*(1-R2)^2*(df^2)/((n.obs^2-1)*(n.obs+3))
     	                     SE =sqrt(SE2)
     	                     F <- R2*df/(k*(1-R2))
     	                     pF <- 1 - pf(F,k,df)
     	                     shrunkenR2 <- 1-(1-R2)*(n.obs-1)/df    	                     }
     	
     	beta <-  t(t(beta) * sqrt(diag(C)[y]))/sqrt(diag(ac.matrix)) #this puts the betas into the raw units
        
     	if(is.null(n.obs)) {mat.regress <- list(beta=round(beta,digits),R=round(sqrt(R2),digits),R2=round(R2,digits))} else {
     	              mat.regress <- list(beta=round(beta,digits),se=round(se,digits),t=round(tvalue,digits),Probability = round(prob,digits),R=round(sqrt(R2),digits),R2=round(R2,digits),shrunkenR2 = round(shrunkenR2,digits),seR2 = round(SE,digits),F=round(F,digits),probF= round(pF,digits+1),df=c(k,df))}
     	return(mat.regress)
     	}
#modified July 12,2007 to allow for NA in the overall matrix
#modified July 9, 2008 to give statistical tests
#modified yet again August 15 , 2008to convert covariances to correlations
