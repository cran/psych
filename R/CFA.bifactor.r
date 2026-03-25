##Developed 12/1/2025 - 1/18/2026 
	# Closed form estimation of Confirmatory factoring
	#code taken from equations for Spearman method by
	#CFA from Dhaene and Rosseel (SEM, 2024 (31) 1-13)
	
	# D and R discuss multiple approaches but recommend the Spearman algorithm as best
	# follows a paper by Guttman 
#modified 3/6/26 to work with CFA using iterations	
	

#########	
	CFA.bifactor <- function(
	   model=NULL,
	   x=NULL,
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
	Grice=FALSE ,
	n.iter=1) {
	
	cl <- match.call()  #save for result
	keys <- model
	#first some housekeeping
	r <- x
	 if(is.list(model)) {X <- make.keys(r,model)} else {if(is.matrix(model)) {X <- model} else {X <- lavParse(model)}  #just a matrix of 1, 0, -1
	  X[X >0] <- 1
	  X[X < 0] <- -1  #just a matrix of 1, 0, -1
	 }  #converts the keys variables to 1, 0 matrix (if not already in that form)
	
	if( is.list(model)) {select <- sub("-","",unlist(keys))} else {
		select<- rownames(X)}  # to speed up scoring one set of scales from many items
	if(any( !(select %in% colnames(x)) )) {
	 cat("\nVariable names in keys are incorrectly specified. Offending items are ", select[which(!(select %in% colnames(x)))],"\n")
	 stop("I am stopping because of improper input.   See above for a list of bad item(s). ")}
			if(!all) {x <-x[,select] 
			X <-X[select,,drop=FALSE]}
			nfact <- ncol(X)
	#now, we use the keys to reverse code some items
	
	 
	 #do we need to find correlations, and if so, to flip them?
	  original.data <- x   #save these for second pass 
	 if(!isCovariance(x)) {n.obs <- nrow(x)
	
	  switch(cor, 
		   cor = {if(!is.null(weight))  {r <- cor.wt(x,w=weight)$r} else  {
										 r <- cor(x,use=use)}
										 },
		   cov = {if(!is.null(weight))  {r <- cor.wt(r,w=weight,cor=FALSE)$r} else  {
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
			gf <- CFA(x=r,all=TRUE, n.obs=n.obs,n.iter=n.iter)
		
			gf.ci <- gf$ci
			gf<- gf$cfa
			null.chisq <- gf$stats$null.chisq   #this is null of original matrix
	n.obs <- gf$stats$n.obs
	resid.g <- gf$residual
	diag(resid.g) <- 1
	group <- CFA(keys,resid.g,all=all,n.obs=n.obs,Grice=Grice,n.iter=n.iter)

	gr.ci <- group$ci
	group <- group$cfa  
	bifactor <- cbind(gf$loadings,gr.ci$median)
	h2 <- rowSums(bifactor^2)
	colnames(bifactor)[1] <- "g"
	om.stats <- omegaStats(r, gload=gf$loading,group=gr.ci$median,h2=h2,Phi=group$Phi)    #to calculate omega statistics
	group$tats$dof <- dof <- nvar*(nvar-1)/2 - nvar - sum(abs(X)) 
	R2 <- fsi(bifactor,r=r,Grice=Grice)
	#loadings replaced with median loadings from the iterated CFA solution
	 
	scores <- factor.scores(x=original.data,f=cbind(gf.ci$median,gr.ci$median), Phi=NULL,method=method, missing=missing,impute=impute, Grice=Grice)
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
	 gf.ci = gf.ci$ci ,
	 gr.ci= gr.ci$ci,
	  X=X,	
	    Call=cl)
	    } else {   #the hierarachical case
	    
		group <- CFA(model,x=r,all=all,n.obs=n.obs,Grice=Grice,n.iter=n.iter)
		gr.ci <- group$ci
	    group <- group$cfa  
		gload <- CFA(x=group$Phi,all=all,n.obs=group$stats$n.obs,Grice=Grice,n.iter=n.iter)
		gf.ci <- gload$ci
		gload <- gload$cfa
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
    gf.ci = gf.ci$ci ,
	 gr.ci= gr.ci$ci,
	 X=X,	
	  
	     Call=cl)}

	class(result) <- c("psych", "CFA.bifactor")
	
	return(result) 
	}
	
