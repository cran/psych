	#Developed 12/1/2025 - 1/18/2026 
	# Closed form estimation of Confirmatory factoring
	#code taken from equations for Spearman method by
	#CFA from Dhaene and Rosseel (SEM, 2024 (31) 1-13)
	
	# D and R discuss multiple approaches but recommend the Spearman algorithm as best
	# follows a paper by Guttman 
	#still have problems in estimating factor indeterminancy
	
	CFA <- function(model=NULL,r=NULL,all=FALSE,cor="cor", use="pairwise", n.obs=NA, orthog=FALSE, weight=NULL, correct=0, method="regression", 
	missing=FALSE,impute="none" ,Grice=FALSE) {
	cl <- match.call()
	# 3 ways of defining the model: a keys list ala scoreIems, matrix, or lavaan style
	keys <- model
	if(is.null(r)) {r <- model
		 model <- NULL
		  X <- matrix(rep(1,ncol(r)),ncol=1)
		rownames(X) <- colnames(r)} else { 
	 if(is.null(model) ) {
	 X <- matrix(rep(1,ncol(r)),ncol=1)
	 rownames(X) <- rownames(r)} else {  
			if(is.list(model)) {X <- make.keys(r,model)} else {if(is.matrix(model)) {X <- model} else {X <- lavParse(model)}
			}
	
	  X[X >0] <- 1
	  X[X < 0] <- -1  #just a matrix of 1, 0, -1
	 }  #converts the keys variables to 1, 0 matrix (if not already in that form)
	} 
	
	if(all) {xx <- matrix(0,ncol=ncol(X),ncol(r))
			 rownames(xx) <- rownames(r)
			 xx[rownames(X),] <- X
			 colnames(xx) <- colnames(X)
			 X <- xx }
	
	if( is.list(model)) {select <- sub("-","",unlist(model))} else {
		select<- rownames(X)}  # to speed up scoring one set of scales from many items
	if(any( !(select %in% colnames(r)) )) {
	 cat("\nVariable names in keys are incorrectly specified. Offending items are ", select[which(!(select %in% colnames(r)))],"\n")
	 stop("I am stopping because of improper input.   See above for a list of bad item(s). ")}
			if(!all) {X <-X[select,,drop=FALSE]}   
	
	 if(!isCovariance(r)) {n.obs <- nrow(r)   #use isCpvariance because some datasets fail correlation test
	 if(!all) {r <-r[,select,drop=FALSE]} 
	 original.data <- r   #save these for scores 
	  switch(cor, 
		   cor = {if(!is.null(weight))  {r <- cor.wt(r,w=weight)$r} else  {
										 r <- cor(r,use=use)}
										 },
		   cov = {if(!is.null(weight))  {r <- cor.wt(r,w=weight,cor=FALSE)$r} else  {
										 r <- cov(r,use=use)}
										   covar <- TRUE},  #fixed 10/30/25
		   wtd = { r <- cor.wt(r,w=weight)$r},
		   spearman = {r <- cor(r,use=use,method="spearman")},
		   kendall = {r <- cor(r,use=use,method="kendall")},
		   tet = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
		   poly = {r <- polychoric(r,correct=correct,weight=weight)$rho},
		   tetrachoric = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
		   polychoric = {r <- polychoric(r,correct=correct,weight=weight)$rho},
		   mixed = {r <- mixedCor(r,use=use,correct=correct)$rho},
		   Yuleb = {r <- YuleCor(r,,bonett=TRUE)$rho},
		   YuleQ = {r <- YuleCor(r,1)$rho},
		   YuleY = {r <- YuleCor(r,.5)$rho } 
		   )
			
		 R <- S <- r
		 
			 } else { if(!all) {r <-r[select,select,drop=FALSE]}  
			 R <- S <- r
		 if(is.na(n.obs)){ n.obs <- 100
			 cat("\n Number of observations not specified. Arbitrarily set to 100.\n")}
		 original.data  <- r}
		 
		 
		
		diagS <- diag(S)  #save this for later
	
	#flip correlations to go along with directions of keys
	
	flip  <- X %*% rep(1,ncol(X))    #combine all keys into one   to allow us to select the over all model
	if(is.null(model)) {
	p1 <- principal(r,scores=FALSE)}
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
	
	
	scores <- factor.scores(x=original.data,f=L,Phi=P,method=method, missing=missing,impute=impute, Grice=Grice)

	#residual <- R-L %*% t(Lambda)
	if(!is.null(P)) {residual <- R-L %*% P %*% t(L)} else {residual <- R-L %*% t(L)} 
	rownames(X) <- rownames(L) 
	colnames(L) <- colnames(X)
	stats <- fa.stats(r=R,f=L ,phi=P,n.obs=n.obs, dof=dof, Grice=Grice)
	result <- list(loadings=L,Phi=P,Pattern=Lambda, communality=H2, dof = dof,
					stats=stats, r=r, residual=residual,
	  #rms = stats$rms,
	 # RMSEA =stats$RMSEA[1],
	  #chi = stats$chi,
	 # BIC  = stats$BIC,
	  #STATISTIC = stats$STATISTIC,
	 # null.chisq = stats$null.chisq,
	 # null.dof = stats$null.dof,
	  scores= scores$scores,
	 weights=scores$weights,
	 model = X,   #the original model as parsed
	  Call=cl
	  )
	class(result) <- c("psych", "CFA")
	return(result)
	}



##########	
	sumLower <- function(R,diag=FALSE){
	   return( sum(R[lower.tri(R,diag)]))}
	   
	
	
#############	
#added the positive option to treat non positive definite matrices
	Spearman.h2 <- function(R,positive=TRUE){   #from Harman chapter 7 
	  diag(R) <- 0
	  if(positive) {sumr <- colSums(abs(R))} else {
		  sumr <- colSums(R)}   #Harman doesn't consider negatives
	  sumr2 <- colSums(R^2)
	  if(positive) { h2 <- (sumr^2 - sumr2)/(2*(sumLower(abs(R))-sumr))}  else { 
		   h2 <- (sumr^2 - sumr2)/(2*(sumLower(R)-sumr))}
	return(h2)}
	  
	  
	  
################

#CFA.bifactor

#########	
	CFA.bifactor <- function(
	   model=NULL,
	   r,
	   all=FALSE,
	   g=FALSE,
	   cor="cor", 
	   use="pairwise",
	    n.obs=NA,
	    orthog=FALSE, 
	    weight=NULL, 
	    correct=0, 
	method="regression", 
	missing=FALSE,
	impute="none",
	Grice=FALSE ) {
	
	cl <- match.call()  #save for result
	keys <- model
	#first some housekeeping
	 if(is.list(model)) {X <- make.keys(r,model)} else {if(is.matrix(model)) {X <- model} else {X <- lavParse(model)}  #just a matrix of 1, 0, -1
	  X[X >0] <- 1
	  X[X < 0] <- -1  #just a matrix of 1, 0, -1
	 }  #converts the keys variables to 1, 0 matrix (if not already in that form)
	
	if( is.list(model)) {select <- sub("-","",unlist(keys))} else {
		select<- rownames(X)}  # to speed up scoring one set of scales from many items
	if(any( !(select %in% colnames(r)) )) {
	 cat("\nVariable names in keys are incorrectly specified. Offending items are ", select[which(!(select %in% colnames(r)))],"\n")
	 stop("I am stopping because of improper input.   See above for a list of bad item(s). ")}
			if(!all) {r <-r[,select] 
			X <-X[select,,drop=FALSE]}
			nfact <- ncol(X)
	#now, we use the keys to reverse code some items
	
	 
	 #do we need to find correlations, and if so, to flip them?
	  original.data <- r   #save these for second pass 
	 if(!isCovariance(r)) {n.obs <- nrow(r)
	
	  switch(cor, 
		   cor = {if(!is.null(weight))  {r <- cor.wt(r,w=weight)$r} else  {
										 r <- cor(r,use=use)}
										 },
		   cov = {if(!is.null(weight))  {r <- cor.wt(r,w=weight,cor=FALSE)$r} else  {
										 r <- cov(r,use=use)}
										   covar <- TRUE},  #fixed 10/30/25
		   wtd = { r <- cor.wt(r,w=weight)$r},
		   spearman = {r <- cor(r,use=use,method="spearman")},
		   kendall = {r <- cor(r,use=use,method="kendall")},
		   tet = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
		   poly = {r <- polychoric(r,correct=correct,weight=weight)$rho},
		   tetrachoric = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
		   polychoric = {r <- polychoric(r,correct=correct,weight=weight)$rho},
		   mixed = {r <- mixedCor(r,use=use,correct=correct)$rho},
		   Yuleb = {r <- YuleCor(r,,bonett=TRUE)$rho},
		   YuleQ = {r <- YuleCor(r,1)$rho},
		   YuleY = {r <- YuleCor(r,.5)$rho } 
		   )
		   }
	  #flip correlations to go along with directions of keys
	nvar <- ncol(r)
	flip  <- X %*% rep(1,ncol(X))    #combine all keys into one   to allow us to select the over all model
	
	flipd <- diag(nrow=ncol(r), ncol=ncol(r))
	rownames(flipd) <- colnames(flipd)<- colnames(r)
	diag(flipd)[rownames(flipd) %in% rownames(flip)] <- flip
	#in the interesting case of no general factor we need to flip some variables
	
	#We do not need to flip them here, because they are flipped in CFA
	#r  <- flipd %*% r %*% t(flipd)   #negatively keyed variables are flipped   
		
	if(!g) {     	
			gf <- CFA(r=r,all=TRUE, n.obs=n.obs)
			null.chisq <- gf$stats$null.chisq   #this is null of original matrix
	n.obs <- gf$stats$n.obs
	resid.g <- gf$residual
	diag(resid.g) <- 1
	group <- CFA(keys,resid.g,all=all,n.obs=n.obs,Grice=Grice)
	bifactor <- cbind(gf$loadings,group$loadings)
	h2 <- rowSums(bifactor^2)
	colnames(bifactor)[1] <- "g"
	om.stats <- omegaStats(r, gload=gf$loading,group= group$loadings,h2=h2,Phi=group$Phi)    #to calculate omega statistics
	group$tats$dof <- dof <- nvar*(nvar-1)/2 - nvar - sum(abs(X)) 
	R2 <- fsi(bifactor,r=r,Grice=Grice)
	scores <- factor.scores(x=original.data,f=cbind(gf$loading,group$loadings),Phi=NULL,method=method, missing=missing,impute=impute, Grice=Grice)
	stats=group$stats
	group$stats$dof <- dof
	result <- list(loadings=bifactor,communality=h2,n.obs=n.obs,
	    stats=group$stats,
	     dof=dof,
	    Phi=group$Phi, R2=R2,
	    omegaStats = om.stats, 
	    r = r, 
	    scores= scores$scores,
	 weights=scores$weights,
	  rms = stats$rms,
	  RMSEA =stats$RMSEA[1],
	  chi = stats$chi,
	  BIC  = stats$BIC,
	  STATISTIC = stats$STATISTIC,
	   null.chisq =  null.chisq,
	  null.dof = stats$null.dof, 
	    Call=cl)
	    } else {   #the hierarachical case
		group <- CFA(model,r=r,all=all,n.obs=n.obs,Grice=Grice)
		gload <- CFA(r=group$Phi,all=all,n.obs=group$stats$n.obs,Grice=Grice)
		loadings <- group$loadings
		h2 <- rowSums(loadings^2)
		n.obs<- group$n.obs
		null.chisq <- group$stats$null.chisq   #this is null of original matrix
		dof <- nvar*(nvar-1)/2 - sum(abs(X)) - nfact   # 
		omegaStats <- NULL  # fix omegaStats for this case
		
		scores <- factor.scores(x=original.data,f=cbind(group$loadings),
			Phi=NULL,method=method, 
			missing=missing,impute=impute, Grice=Grice)
		stats=group$stats
		stats$dof <- dof
		stats$STATISTIC <- stats$STATISTIC + gload$stats$STATISTIC
		stats$BIC <- stats$BIC + gload$stats$BIC
	result <- list(loadings=loadings,communality=h2 ,gload=gload, 
	   Phi = group$Phi, n.obs=n.obs,
	   stats=stats,
	   dof=dof,
	   omegaStats = omegaStats,
	   r=r,	  
	   scores= scores$scores,
	   weights=scores$weights,
	    rms = stats$rms,
	  RMSEA =stats$RMSEA[1],
	  chi = stats$chi,
	  BIC  = stats$BIC,
	  STATISTIC = stats$STATISTIC, 
	  null.chisq = stats$null.chisq,
	  null.dof = stats$null.dof,	 
	     Call=cl)}

	class(result) <- c("psych", "CFA.bifactor")
	
	return(result) 
	}
	
