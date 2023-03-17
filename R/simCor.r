#heavily modified August 15, 2020 to allow for skewing of data, etc.
"simCor" <- 
"sim.correlation" <- function(R,n=1000,data=FALSE,scale=TRUE, skew=c("none","log","lognormal","sqrt","abs"),vars=NULL,latent=FALSE,quant=NULL) {
     eX <- eigen(R)
     nvar <- ncol(R) 
     if(missing(skew)) skew <- "none"
     if(missing(vars)) vars <- 1:nvar
     if(skew=="lognormal") {observed <- matrix(rlnorm(nvar * n),n)} else {
     observed <- matrix(rnorm(nvar * n),n) }  #create the data
    
    if(latent) {  #apply the transformations first   -- probably not very helpful
       min.obs <- apply(observed,2,min)
        if(skew == "abs") { observed[,vars]  <- abs(observed[,vars]) } else {
       observed <- t(t(observed)-min.obs + .1)
       } 
     switch(skew,
     	none = observed <- observed,
     	log = observed[,vars] <- log(observed[,vars]),
     	lognormal = observed[,vars] <- observed[,vars],
    	 sqrt = observed[,vars]  <-sqrt(observed[,vars]),
     	abs = observed <- observed  #they are already done
     )
      if(!missing(quant)) {
    cut <- quantile(observed,quant)
    observed[observed[,vars] < cut] <- 0
    observed[observed[,vars] >= cut] <- 1}
    } 
    
    observed <- -t( eX$vectors %*% diag(sqrt(pmax(eX$values, 0)), nvar) %*%  t(observed)) #apply the eigen solution
    if(scale) observed <- scale(observed)
   
    if(!latent) {
     min.obs <- apply(observed,2,min)
        if(skew == "abs") { observed[,vars] <- abs(observed[,vars]) } else {
       observed <- t(t(observed)-min.obs + .1)
       } 
     switch(skew,
    	 none = observed <- observed,
    	 log = observed[,vars] <- log(observed[,vars]),
    	 lognormal = observed[,vars] <- observed[,vars],
     	sqrt = observed[,vars]  <-sqrt(observed[,vars]),
     	abs = observed[,vars] <- observed[,vars]
     )
     if(scale) observed <- scale(observed)
     
     
       if(!missing(quant)) {
       cut <- rep(0,nvar)
       if(length(quant) == 1) quant <- rep(quant,nvar)
       for(i in 1:length(vars)) {cut[i] <- quantile(observed[,vars[i]],quant[i])
    observed[observed[,vars[i]] < cut[i] ,vars[i]] <- -999
    observed[observed[,vars[i]] >= cut[i],vars[i]] <- 1}
    }
    observed[observed== -999] <- 0
    }
    if(!is.null(colnames(R)) ) {colnames(observed) <- colnames(R)} else {colnames(R) <- paste0("V",1:nvar)}
    if(data) {result <- observed} else {
  	 result <- cor(observed)}
  	 return(result)}
