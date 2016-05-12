#vss is just an alias to VSS to be consistent with naming conventions
#added the RMSEA, BIC, SABIC and complexity criteria 1/27/14
"VSS" <-
function (x,n=8,rotate="varimax",diagonal=FALSE,fm="minres",n.obs=NULL,plot=TRUE,title="Very Simple Structure",use="pairwise",cor="cor",...)     #apply the Very Simple Structure Criterion for up to n factors on data set x
 {vss(x=x,n=n,rotate=rotate,diagonal=diagonal,fm=fm,n.obs=n.obs,plot=plot,title=title,use=use,cor=cor,...)  }

"vss" <-
function (x,n=8,rotate="varimax",diagonal=FALSE,fm="minres",n.obs=NULL,plot=TRUE,title="Very Simple Structure",use="pairwise",cor="cor",...)     #apply the Very Simple Structure Criterion for up to n factors on data set x
 
 #x is a data matrix
  #n is the maximum number of factors to extract  (default is 8)
  #rotate is a string "none" or "varimax" for type of rotation (default is "varimax"
  #diagonal is a boolean value for whether or not we should count the diagonal  (default=FALSE)
  # ... other parameters for factanal may be passed as well  
  #e.g., to do VSS on a covariance/correlation matrix with up to 8 factors and 3000 cases:
  #VSS(covmat=msqcovar,n=8,rotate="none",n.obs=3000)
   { 
  cl <- match.call()
  if (rotate=="oblimin") {if(!requireNamespace('GPArotation')) {stop("You must have GPArotation installed to use oblimin rotation")}}
  old_rotate=rotate  #used to remember which rotation to use
  
            #start Function definition
  #first some preliminary functions
  #complexrow sweeps out all except the c largest loadings
  #complexmat applies complexrow to the loading matrix
 

complexrow <- function(x,c)     #sweep out all except c loadings
    {  n=length(x)          	#how many columns in this row?
       temp <- x                #make a temporary copy of the row
       x <- rep(0,n)            #zero out x
       for (j in 1:c) 
       {
       	locmax <- which.max(abs(temp))                     #where is the maximum (absolute) value
      	 x[locmax] <- sign(temp[locmax])*max(abs(temp))    #store it in x
       	temp[locmax] <- 0                                  #remove this value from the temp copy
       }
     return(x)                                             #return the simplified (of complexity c) row 
    }
    
 complexmat <- function(x,c)           #do it for every row   (could tapply somehow?)
	{
	nrows <- dim(x)[1]
	ncols <- dim(x)[2]
	for (i in 1:nrows)
   		{x[i,] <- complexrow(x[i,],c)}   #simplify each row of the loading matrix
 	return(x)
	 }  
	 
 map <- function(x,n) { 
 nvar <- dim(x)[2]
 min.partial <- rep(NA,n)
 e <- eigen(x)
 evect <- e$vectors
 comp <- evect %*% diag(sqrt(e$values))
 if( n >= nvar) {n1 <- nvar -1} else {n1 <- n}
 for (i in 1:n1) {
 c11.star <- x - comp[,1:i] %*% t(comp[,1:i]) 
 d <- diag(1/sqrt(diag(c11.star)))
 rstar <- d %*% c11.star %*% d
 diag(rstar) <- 0
 min.partial[i] <- sum(rstar * rstar)  /(nvar*(nvar-1))
 }
 return(min.partial)
 }
 
  if(dim(x)[2] < n) n <- dim(x)[2]  
  #now do the main Very Simple Structure  routine

  complexfit <- array(0,dim=c(n,n))        #store these separately for complex fits
  complexresid <-  array(0,dim=c(n,n))
  
  
  vss.df <- data.frame(dof=rep(0,n),chisq=NA,prob=NA,sqresid=NA,fit=NA,RMSEA=NA,BIC=NA,SABIC=NA,complex=NA,eChisq=NA,SRMR=NA,eCRMS=NA,eBIC=NA) #keep the basic results here 
 
  if (dim(x)[1]!=dim(x)[2]) { n.obs <- dim(x)[1]
     switch(cor, 
       cor = {x <- cor(x,use=use)},
       cov = {x <- cov(x,use=use) 
              covar <- TRUE},
       tet = {x <- tetrachoric(x)$rho},
       poly = {x <- polychoric(x)$rho},
       mixed = {x <- mixed.cor(x,use=use)$rho},
       Yuleb = {x <- YuleCor(x,,bonett=TRUE)$rho},
       YuleQ = {x <- YuleCor(x,1)$rho},
       YuleY = {x <- YuleCor(x,.5)$rho } 
       )
              # x <- cor(x,use="pairwise")  The case statement allows many types of correlations
                }  else {if(!is.matrix(x)) x <- as.matrix(x)}
              # if given a rectangular 
  if(is.null(n.obs)) {message("n.obs was not specified and was arbitrarily set to 1000.  This only affects the chi square values.")
        n.obs <- 1000}
        
        map.values <- map(x,n)
 if (n >  dim(x)[2]) {n <- dim(x)[2]}       #in cases where there are very few variables

 for (i in 1:n)                            #loop through 1 to the number of factors requested
 { PHI <- diag(i) 
 if(i<2) {(rotate="none")} else {rotate=old_rotate} 
   if(!(fm=="pc")) { 
   		f <- fa(x,i,rotate=rotate,n.obs=n.obs,warnings=FALSE,fm=fm,scores="none",cor=cor,...)   #do a factor analysis with i factors and the rotations specified in the VSS call
 	 if (i==1)
  		 {original <- x         #just find this stuff once
		 sqoriginal <- original*original    #squared correlations
		 totaloriginal <- sum(sqoriginal) - diagonal*sum(diag(sqoriginal) )   #sum of squared correlations - the diagonal
		}} else {f <- principal(x,i)
	    if (i==1)
  			 {original <- x       #the input to pc is a correlation matrix, so we don't need to find it again
			 sqoriginal <- original*original    #squared correlations
		 	totaloriginal <- sum(sqoriginal) - diagonal*sum(diag(sqoriginal) )   #sum of squared correlations - the diagonal
			 }
		if((rotate=="varimax") & (i>1)) {f <- varimax(f$loadings) 
		                                 PHI <- diag(i)} else {
		if(((rotate=="promax")| (rotate=="Promax") )& (i>1))  {f <- Promax(f$loadings) 
          								 PHI <- f$Phi} else {
		if((rotate=="oblimin")& (i>1))  {f <- GPArotation::oblimin(f$loadings)
		    U <- f$Th
           phi <- t(U) %*% U
           PHI <- cov2cor(phi)
		}}
	     }}
		
 	load <- as.matrix(f$loadings )                    #the loading matrix
 	
   	model <- load %*% PHI %*%  t(load)       #reproduce the correlation matrix by the factor law R=  FF'
 	residual <- original-model              #find the residual  R* = R - FF'
 	sqresid <- residual*residual            #square the residuals
 	totalresid <- sum(sqresid)- diagonal * sum(diag(sqresid) )      #sum squared residuals - the main diagonal
 	fit <- 1-totalresid/totaloriginal       #fit is 1-sumsquared residuals/sumsquared original     (of off diagonal elements

 	if ((fm !="pc")) {
 			vss.df[i,"dof"] <- f$dof                   #degrees of freedom from the factor analysis
 			vss.df[i,"chisq"] <- f$STATISTIC             #chi square from the factor analysis
 			vss.df[i,"prob"] <- f$PVAL                  #probability value of this complete solution\
 			vss.df[i,"eChisq"] <- f$chi
 			vss.df[i,"SRMR"] <- f$rms
 			vss.df[i,"eRMS"] <- f$rms
 			vss.df[i,"eCRMS"] <- f$crms
 			vss.df[i,"eBIC"] <- f$EBIC
 			
 			if(!is.null(f$RMSEA)) {vss.df[i,"RMSEA"] <- f$RMSEA[1]} else {vss.df[i,"RMSEA"] <- NA}
 			if(!is.null(f$BIC)) {vss.df[i,"BIC"] <- f$BIC} else {vss.df[i,"BIC"] <- NA}
 			if(!is.null(f$SABIC)) {vss.df[i,"SABIC"] <- f$SABIC} else {vss.df[i,"SABIC"] <- NA}
 			if(!is.null(f$complexity)) {vss.df[i,"complex"] <- mean(f$complexity)} else {vss.df[i,"complex"] <- NA}
 			 }
  	vss.df[i,"sqresid"] <- totalresid              #residual given complete model
  	vss.df[i,"fit"] <- fit                     #fit of complete model
  	
  	
  	
  	
     #now  do complexities -- how many factors account for each item
 
  for (c in 1:i)   
  	
  	 { 
  	 	simpleload <- complexmat(load,c)             #find the simple structure version of the loadings for complexity c
  		model <- simpleload %*% PHI %*% t(simpleload)           #the model is now a simple structure version  R ? SS'
  		residual <- original- model                   #R* = R - SS'       
  		sqresid <- residual*residual
  		totalsimple <- sum(sqresid) -diagonal * sum(diag(sqresid))    #default is to not count the diagonal 
  		simplefit <- 1-totalsimple/totaloriginal
  		complexresid[i,c] <-totalsimple
  		complexfit[i,c] <- simplefit
  	 }
  	
}     #end of i loop for number of factors



vss.stats <- data.frame(vss.df,cfit=complexfit,cresidual=complexresid)
if (plot) VSS.plot(vss.stats,title=title)
vss.results <- list(title=title,map=map.values,cfit.1=complexfit[,1],cfit.2= complexfit[,2],vss.stats=vss.stats,call=cl)
class(vss.results) <- c("psych" ,"vss")
return(vss.results)
   
    }     #end of VSS function


"nfactors" <-
function(x,n=20,rotate="varimax",diagonal=FALSE,fm="minres",n.obs=NULL,title="Number of Factors",pch=16,use="pairwise",cor="cor",...) {
vs <- vss(x=x,n=n,rotate=rotate,diagonal=diagonal,fm=fm,n.obs=n.obs,plot=FALSE,title=title,use=use,cor=cor,...) 
old.par <- par(no.readonly = TRUE) # save default, for resetting... 
on.exit(par(old.par))     #and when we quit the function, restore to original values
op <- par(mfrow=c(2,2))
x <- vs$vss.stats
n = dim(x)
    plot(x$cfit.1, ylim = c(0, 1), typ = "b", ylab = "Very Simple Structure Fit", 
        xlab = "Number of Factors",main="Very Simple Structure",pch=49)
          lines(x$cfit.1)
    x$cfit.2[1] <- NA
    x$cfit.3[1] <- NA
    x$cfit.3[2] <- NA
    lines(x$cfit.2)
    points(x$cfit.2,pch=50)
     lines(x$cfit.3)
    points(x$cfit.3,pch=51)
   
plot(vs$vss.stats[,"complex"],xlab="Number of factors",ylab="Complexity",typ="b",main="Complexity",pch=pch,...)
plot(vs$vss.stats[,"eBIC"],xlab="Number of factors",ylab="Empirical BIC",typ="b",main="Empirical BIC",pch=pch,...)
plot(vs$vss.stats[,"SRMR"],xlab="Number of factors",ylab="SRMR",typ="b",main="Root Mean Residual",pch=pch,...)
results <- list(title=title,map=vs$map,vss.stats=vs$vss.stats[,1:16],call=vs$call)
class(results) <- c("psych","vss")
return(results)
}


