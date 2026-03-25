	#Developed 12/1/2025 - 1/18/2026 
	#revised March 2, 2026 to do iterations
	
	# Closed form estimation of Confirmatory factoring
	#code taken from equations for Spearman method by
	#CFA from Dhaene and Rosseel (SEM, 2024 (31) 1-13)
	
	# D and R discuss multiple approaches but recommend the Spearman algorithm as best
	# follows a paper by Guttman 
	#still have problems in estimating factor indeterminancy
	
	CFA <- function(model=NULL,x=NULL,all=FALSE,cor="cor", use="pairwise", n.obs=NA, orthog=FALSE, weight=NULL, correct=0, method="regression", 
	missing=FALSE,impute="none" ,Grice=FALSE, n.iter =100) {
	cl <- match.call()
	# 3 ways of defining the model: a keys list ala scoreIems, matrix, or lavaan style
	keys <- model
	if(is.null(x)) {x <- model
		 model <- NULL
		  X <- matrix(rep(1,ncol(x)),ncol=1)
		rownames(X) <- colnames(x)} else { 
	 if(is.null(model) ) {
	 X <- matrix(rep(1,ncol(x)),ncol=1)
	 rownames(X) <- rownames(x)} else {  
			if(is.list(model)) {X <- make.keys(x,model)} else {if(is.matrix(model)) {X <- model} else {X <- lavParse(model)}
			}
	
	  X[X >0] <- 1
	  X[X < 0] <- -1  #just a matrix of 1, 0, -1
	 }  #converts the keys variables to 1, 0 matrix (if not already in that form)
	} 
	nvar <- nrow(X)
	if(all) {xx <- matrix(0,ncol=ncol(X),ncol(x))
			 rownames(xx) <- rownames(x)
			 xx[rownames(X),] <- X
			 colnames(xx) <- colnames(X)
			 X <- xx }
	
	if( is.list(model)) {select <- sub("-","",unlist(model))} else {
		select<- rownames(X)}  # to speed up scoring one set of scales from many items
	if(any( !(select %in% colnames(x)) )) {
	 cat("\nVariable names in keys are incorrectly specified. Offending items are ", select[which(!(select %in% colnames(x)))],"\n")
	 stop("I am stopping because of improper input.   See above for a list of bad item(s). ")}
			if(!all) {X <-X[select,,drop=FALSE]}   
	
	if(isCovariance(x)) { if(!all) {x <-x[select,select,drop=FALSE]}     #use isCovariance because some datasets fail correlation test
	               R <- S <- x 
	               nvar <- ncol(x)
	                if(is.na(n.obs)){ n.obs <- 100
			 cat("\n Number of observations not specified. Arbitrarily set to 100.\n")}
		 original.data  <- x
		 } else {
	 n.obs <- nrow(x)
	 if(!all) {x <-x[,select,drop=FALSE]} 
	 original.data <- x   #save these for scores  -- use for iterations 
	 nvar <- ncol(x)}
	 
	 #first do one pass to get the estimates
	 #then do this n.iter times to get the errors
	 
	 cf.original <- cfa.main(X,x,cor,use,correct,weight,orthog,method,n.obs,model,original.data,Grice=Grice,cl=cl) 
	 
	#replicateslist <- parallel::mclapply(1:n.iter,function(xx) { 
   replicateslist <- lapply(1:n.iter,function(xx) {



 if(isCovariance(x)) {#create data sampled from multivariate normal with observed correlation
                       nvar <- ncol(x)
                                      mu <- rep(0, nvar)
                                      #X <- mvrnorm(n = n.obs, mu, Sigma = r, tol = 1e-06, empirical = FALSE)
                                      #the next 3 lines replaces mvrnorm (taken from mvrnorm, but without the checks)
                                      eX <- eigen(x)
                                      XX <- matrix(rnorm(nvar * n.obs),n.obs)
                                      x <-  t(eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(XX))
                            } else {x <- original.data[sample(n.obs,n.obs,replace=TRUE),]}
                            
  fs <- cfa.main(X,x,cor,use,correct,weight,orthog,method,n.obs,model=model,original.data=original.data,Grice=Grice,cl=cl) #call CFA with the appropriate parameters
                  
                  if(!is.null(fs$Phi)) {
                  replicates <- list(loadings=fs$loadings,Phi=fs$Phi[lower.tri(fs$Phi)])} else {replicates <- list(loadings=fs$loadings)} 
                  }
                  )
                        #end iterative function             
 
 #ci <- list( replicates=replicateslist)

 ci <- summarize.replicates(cf.original, replicateslist, X=X)
 rep.load <- ci$median
 #now,do the stats and scores based upon the replicated loadings
 nfact <- ncol(X)
 dof <- nvar * (nvar-1)/2 - nvar - (!orthog) * nfact*(nfact-1)/2 #number of correlations - number of factor loadings estimated - correlations
 stats <- fa.stats(r=cf.original$r,f=rep.load ,phi=cf.original$Phi,n.obs=n.obs, dof=dof, Grice=Grice)
 cf.original$stats <- stats
 results <- list(cfa=cf.original,ci= ci)
 class(results) <- c("psych", "CFA")
 return(results)
}


#the function to iterate does the actual CFA analysis
	  
cfa.main <- function(X, x, cor=cor,  use=use,correct=correct,weight=weight,orthog=orthog,method=method,n.obs=n.obs,model=model,original.data,Grice, 
	missing=FALSE,impute="none",cl ){  
	 if(!isCovariance(x)) {#find the correlations
#we have cleaned and checked the data, ready to do correlations	 
	 
	 
	  switch(cor, 
		   cor = {if(!is.null(weight))  {r <- cor.wt(x,w=weight)$r} else  {
										 r <- cor(x,use=use)}
										 },
		   cov = {if(!is.null(weight))  {r <- cor.wt(x,w=weight,cor=FALSE)$r} else  {
										 r <- cov(x,use=use)}
										   covar <- TRUE},  #fixed 10/30/25
		   wtd = { r <- cor.wt(x,w=weight)$r},
		   spearman = {r <- cor(x,use=use,method="spearman")},
		   kendall = {r <- cor(x,use=use,method="kendall")},
		   tet = {r <- tetrachoric(x,correct=correct,weight=weight)$rho},
		   poly = {r <- polychoric(x,correct=correct,weight=weight)$rho},
		   tetrachoric = {r <- tetrachoric(x,correct=correct,weight=weight)$rho},
		   polychoric = {r <- polychoric(x,correct=correct,weight=weight)$rho},
		   mixed = {r <- mixedCor(x,use=use,correct=correct)$rho},
		   Yuleb = {r <- YuleCor(x,,bonett=TRUE)$rho},
		   YuleQ = {r <- YuleCor(x,1)$rho},
		   YuleY = {r <- YuleCor(x,.5)$rho } 
		   )
			
		 R <- S <- r
		 
		} else {R <- S<- x}	 
			 
		
		 
		
		diagS <- diag(S)  #save this for later
	
	#flip correlations to go along with directions of keys
	
	flip  <- X %*% rep(1,ncol(X))    #combine all keys into one   to allow us to select the over all model
	if(is.null(model)) {
	p1 <- principal(x,scores=FALSE)} else {p1 <- list()
	   p1$loadings <- flip }
	  #  flip <- 1- 2* (p1$loadings < 0)}
		
						
						
	
	flipd <- diag(nrow=ncol(R), ncol=ncol(R))
	rownames(flipd) <- colnames(flipd)<- rownames(R)
	diag(flipd)[rownames(flipd) %in% rownames(flip)] <- flip
	S  <- flipd %*% S %*% t(flipd)   #negatively keyed variables are flipped
	
	mask <- abs(X %*% t(X))
	#S <- mask * S     #select a block diagonal form
	
	
	#need to find the communalities by each factor using the Spearman method
	nfact <- ncol(X)
	nvar <- nrow(X)
	
	dof <- nvar * (nvar-1)/2 - nvar - (!orthog) * nfact*(nfact-1)/2 #number of correlations - number of factor loadings estimated - correlations
	keys <- X
	X <- abs(X)  
	H2 <- NULL
	for (i in 1:nfact) {    #this is the brute force way of doing it, look more carefully at D and R
	  H2 <- c(H2,Spearman.h2(S[X[,i]>0,X[,i]>0]))}   #H2 is the vector of 1 factor communalities
	  
	 if(ncol(S)  == length(H2)) { 
	  diag(S) <- H2   #put the communalities on the diagonal
	} else {for(i in 1:length(H2)) { #pad out the communalities for the groups we have
	   S[names(H2)[i],names(H2)[i]] <- H2[i]}
	} 
	
	if(nrow(X) != nrow(S))  {
	Xplus <- matrix(0,nrow(S)-nrow(X),ncol(X))
	  rownames(Xplus) <- rownames(S)[!(rownames(S) %in% rownames(X))]
	 X <- rbind(X,Xplus) }
	 
	 
	L0 <- t(X)%*% S   #DR equation 4
	P0 <- t(X) %*% S %*% X  #equation 5
	
	if(ncol(P0)>1) {
	   D <- diag(diag(P0))  # variances
	   L <- t(L0) %*% invMatSqrt(D)  #divide by sd eq 6     factor structure
	   P <- invMatSqrt(D) %*% P0 %*% invMatSqrt(D)   #converts P0 to correlations eq 7
	if(orthog) {P <- diag(1,nfact)} 
	Lambda <- L %*% solve(P)  #not simple structure  --- better to use L *X  %*%  P   eq 8
	#Lambda <- (L * X) %*% (P)    #not equation 8, but produces Pattern   
	flipf <- matrix(flip,ncol=nfact,nrow=nrow(Lambda))
	
	Lambda <- Lambda*flipf   #get the signs right 
	L <- L *flipf
	
	colnames(Lambda) <-   colnames(X)# paste0("CF",1:nfact)
	rownames(Lambda) <- colnames(R)
	
	colnames(P) <- rownames(P) <- colnames(Lambda)
	#Ls <- Lambda  * abs(keys)   #forced simple structure
	L <- L * X       #force a simple structure
	
	} else { D= diag(1)
		   L <- Lambda <- as.matrix(sqrt(diag(S)))
		   flip <- 1- 2* (p1$loadings < 0)
		   L <- L * flip     #for the case of negative general factor loadings
		   rownames(L) <- colnames(R)
		   if(is.null(colnames(L))) colnames(L) <- colnames(X) <- "CF1"
		   P <-  NULL
		   Ls <- NULL
		   Phi  <-NULL
		   }
	
	Lambda <- as.matrix(Lambda)
	
	
	
	#scores <- factor.scores(x=original.data,f=L,Phi=P,method=method, missing=missing,impute=impute, Grice=Grice)

	#residual <- R-L %*% t(Lambda)
	if(!is.null(P)) {residual <- R-L %*% P %*% t(L)} else {residual <- R-L %*% t(L)} 
	rownames(X) <- rownames(L) 
	colnames(L) <- colnames(X)
	#stats <- fa.stats(r=R,f=L ,phi=P,n.obs=n.obs, dof=dof, Grice=Grice)
	result <- list(loadings=L,Phi=P,Pattern=Lambda, communality=H2, dof = dof,
					#stats=stats,
					 r=R, residual=residual,
	  #rms = stats$rms,
	 # RMSEA =stats$RMSEA[1],
	  #chi = stats$chi,
	 # BIC  = stats$BIC,
	  #STATISTIC = stats$STATISTIC,
	 # null.chisq = stats$null.chisq,
	 # null.dof = stats$null.dof,
	  #scores= scores$scores,
	# weights=scores$weights,
	 model = X,   #the original model as parsed
	  Call=cl
	  )
 
	  
	class(result) <- c("psych", "CFA")
	return(result)
	}



summarize.replicates <- function(f.original,replicates,X=NULL, p =.05,digits=2) {
n.iter <- length(replicates)
replicates <- matrix(unlist(replicates),nrow=n.iter,byrow=TRUE)
#replicates <- replicates[X>0]
iqr <- apply(replicates,2,quantile, p=c(p/2,.5,1-p/2))
means <- colMeans(replicates,na.rm=TRUE)
sds <- apply(replicates,2,sd,na.rm=TRUE)
nfactors <- ncol(X)
nvar <- nrow(X)

median.mat <- matrix(iqr[2,1:(nvar*nfactors)],ncol=nfactors)
low  <- matrix(iqr[1,1:(nvar*nfactors)],ncol=nfactors)
high <- matrix(iqr[3,1:(nvar*nfactors)],ncol=nfactors)

 median<- median.mat[abs(X)>0]
low <- low[abs(X)>0] 
high <-high[abs(X)>0]

 
   
#ci.lower <-  M + qnorm(p/2) * SD
#ci.upper <- M + qnorm(1-p/2) * SD

ci <- data.frame(estimate=median,lower = low,upper=high)
#class(means) <- "loadings"


fnames  <- colnames(f.original$loadings)
rownames(ci) <-  rownames(f.original$loadings)
colnames(median.mat) <-fnames
rownames(median.mat) <- rownames(ci)

#for(i in 1:nfactors){
#cat("\n Factor ",i, fnames[i], "\n")
#print(round(ci[X[,i]>0,],digits))
#}
result <- list(ci=ci,X=X,median=median.mat)
return(result)
} 




##########	
	sumLower <- function(R) {
	return(sum(R[lower.tri(R)]))}
	
	
#############	
#added the positive option to treat non positive definite mCXatrices
	Spearman.h2 <- function(R,positive=TRUE){   #from Harman chapter 7 
	  diag(R) <- 0
	  if(positive) {sumr <- colSums(abs(R))} else {
		  sumr <- colSums(R)}   #Harman doesn't consider negatives
	  sumr2 <- colSums(R^2)
	  if(positive) { h2 <- (sumr^2 - sumr2)/(2*(sumLower(abs(R))-sumr))}  else { 
		   h2 <- (sumr^2 - sumr2)/(2*(sumLower(R)-sumr))}
		   h2[is.na(h2)] <- 0  #to handle case of all correlations exactly zero
	return(h2)}
	  
	  
	  


