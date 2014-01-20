"VSSem" <-
function (x,n=8,rotate="varimax",diagonal=FALSE,pc="pa",n.obs=NULL,...)     #apply the Very Simple Structure Criterion for up to n factors on data set x
#find the maximum likelihood goodness of fit criterion 
  #x is a data matrix
  #n is the maximum number of factors to extract  (default is 8)
  #rotate is a string "none" or "varimax" for type of rotation (default is "none"
  #diagonal is a boolean value for whether or not we should count the diagonal  (default=FALSE)
  # ... other parameters for factanal may be passed as well  
  #e.g., to do VSS on a covariance/correlation matrix with up to 8 factors and 3000 cases:
  #VSS(covmat=msqcovar,n=8,rotate="none",n.obs=3000)
  
  
 {             #start Function definition
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
    
  #now do the main Very Simple Structure  routine

  complexfit <- array(0,dim=c(n,n))        #store these separately for complex fits
  complexchi <- array(0,dim=c(n,n))
  complexchi2 <- array(0,dim=c(n,n))
   complexdof <- array(0,dim=c(n,n))
  complexresid <-  array(0,dim=c(n,n))
  
  vss.df <- data.frame(dof=rep(0,n),chisq=0,prob=0,sqresid=0,fit=0) #keep the basic results here 
 
  if (dim(x)[1]!=dim(x)[2]) { n.obs <- dim(x)[1]
               x <- cor(x,use="pairwise") }  else {if(!is.matrix(x)) x <- as.matrix(x)}
              # if given a rectangular 
  if(is.null(n.obs)) {message("n.obs was not specified and was arbitrarily set to 1000.  This only affects the chi square values.")
        n.obs <- 1000}
 if (n >  dim(x)[2]) {n <- dim(x)[2]}         #in cases where there are very few variables
 n.variables <- dim(x)[2]     
 for (i in 1:n)                            #loop through 1 to the number of factors requested
 { 
   if(!(pc=="pc")) { if ( pc=="pa") {
   		f <- fa(x,i,fm="pa",rotate=rotate,n.obs=n.obs,...)   #do a factor analysis with i factors and the rotations specified in the VSS call
 	 if (i==1)
  		 {original <- x         #just find this stuff once
		 sqoriginal <- original*original    #squared correlations
		 totaloriginal <- sum(sqoriginal) - diagonal*sum(diag(sqoriginal) )   #sum of squared correlations - the diagonal
		}}  else { 
   	f <- fa(x,i,fm=pc,rotate=rotate,covmat=x,n.obs=n.obs,...)  #do a factor analysis with i factors and the rotations specified in the VSS call
 	 if (i==1)
  		 {original <- x         #just find this stuff once
		 sqoriginal <- original*original    #squared correlations
		 totaloriginal <- sum(sqoriginal) - diagonal*sum(diag(sqoriginal) )   #sum of squared correlations - the diagonal
		}}
	  } else {f <- principal(x,i)
	    if (i==1)
  			 {original <- x       #the input to pc is a correlation matrix, so we don't need to find it again
			 sqoriginal <- original*original    #squared correlations
		 	totaloriginal <- sum(sqoriginal) - diagonal*sum(diag(sqoriginal) )   #sum of squared correlations - the diagonal
			 }
		if((rotate=="varimax") & (i>1)) {f <- varimax(f$loadings)} else {
		if((rotate=="promax") & (i>1))  {f <- promax(f$loadings)}
	     }}
		
 	load <- as.matrix(f$loadings )          #the loading matrix
   	model <- load %*% t(load)               #reproduce the correlation matrix by the factor law R=  FF'
 	residual <- original-model              #find the residual  R* = R - FF'
 	sqresid <- residual*residual            #square the residuals
 	totalresid <- sum(sqresid)- diagonal * sum(diag(sqresid) )      #sum squared residuals - the main diagonal
 	fit <- 1-totalresid/totaloriginal       #fit is 1-sumsquared residuals/sumsquared original     (of off diagonal elements
 	
 	if ((pc!="pc")) {                       #factor.pa reports the same statistics as mle, although the fits are not as good
 			vss.df[i,1] <- f$dof                   #degrees of freedom from the factor analysis
 			vss.df[i,2] <- f$STATISTIC             #chi square from the factor analysis
 			vss.df[i,3] <- f$PVAL                  #probability value of this complete solution
 			 }
  	vss.df[i,4] <- totalresid              #residual given complete model
  	vss.df[i,5] <- fit                     #fit of complete model
  	
  	
  	
  	
     #now  do complexities -- how many factors account for each item
 
  for (c in 1:i)   
  	
  	 { 
  	 	simpleload <- complexmat(load,c)             #find the simple structure version of the loadings for complexity c
  		model <- simpleload%*%t(simpleload)           #the model is now a simple structure version  R ? SS'
  		residual <- original- model                   #R* = R - SS'       
  		sqresid <- residual*residual
  		totalsimple <- sum(sqresid) -diagonal * sum(diag(sqresid))    #default is to not count the diagonal 
  		simplefit <- 1-totalsimple/totaloriginal
  		complexresid[i,c] <-totalsimple
  		complexfit[i,c] <- simplefit
  		
  	#find the chi square value for this level of complexity  (see factor.pa for more details on code)	
  		 diag(model) <- 1   
    model.inv <- solve(model)
    nfactors <- i 
    m.inv.r <- model.inv %*% original
    dof <-  n.variables * (n.variables-1)/2 - n.variables * c + (nfactors *(nfactors-1)/2)
    objective <- sum(diag((m.inv.r))) - log(det(m.inv.r)) -n.variables 
    if (!is.null(n.obs)) {STATISTIC <-  objective * (n.obs-1) -(2 * n.variables + 5)/6 -(2*nfactors)/3
    	if (dof > 0) {PVAL <- pchisq(STATISTIC, dof, lower.tail = FALSE)} else PVAL <- NA}
    complexchi[i,c]  <- STATISTIC
    complexdof[i,c] <- dof
    res1 <- residual
    diag(res1) <- 1
    complexchi2[i,c] <- -(n.obs - n.variables/3 -1.8) *log(det(res1))
  	 }
  	
}     #end of i loop for number of factors


vss.stats <- data.frame(vss.df,cfit=complexfit,chisq=complexchi,complexchi2,complexdof,cresidual=complexresid)
return(vss.stats)
   
    }     #end of VSS function

