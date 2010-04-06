"sim.item" <-
function (nvar = 72 ,nsub = 500, 
    circum = FALSE, xloading =.6, yloading = .6, gloading=0, xbias=0,  ybias = 0,categorical=FALSE, low=-3,high=3,truncate=FALSE,cutpoint=0) 
	{ 
	avloading <- (xloading+yloading)/2
	
	errorweight <- sqrt(1-(avloading^2  + gloading^2))  #squared errors and true score weights add to 1
    g <- rnorm(nsub) 
	truex <- rnorm(nsub)* xloading  +xbias #generate normal true scores for x + xbias
	truey <- rnorm(nsub) * yloading + ybias #generate normal true scores for y + ybias

	if (circum)  #make a vector of radians (the whole way around the circle) if circumplex
	{radia <- seq(0,2*pi,len=nvar+1)  
      rad <- radia[which(radia<2*pi)]        #get rid of the last one
     } else   rad <- c(rep(0,nvar/4),rep(pi/2,nvar/4),rep(pi,nvar/4),rep(3*pi/2,nvar/4))  #simple structure 
        
	error<- matrix(rnorm(nsub*(nvar)),nsub)    #create normal error scores

	#true score matrix for each item reflects structure in radians
	trueitem <- outer(truex, cos(rad)) + outer(truey,sin(rad)) 

	item<- gloading * g +  trueitem  + errorweight*error   #observed item = true score + error score 
    if (categorical) {
        
    	item = round(item)       #round all items to nearest integer value
		item[(item<= low)] <- low     
		item[(item>high) ] <- high   
		}
	if (truncate) {item[item < cutpoint] <- 0  }
	return (item) 
	}  


"sim.rasch" <-
function (nvar = 5 , n = 500, low=-3,high=3,d=NULL, a=1) 
	{ 
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))}
	theta <- rnorm(n)
	item <- matrix(1/(1+exp(a*(-theta %+% t( d)))),n,nvar)
	#now convert these probabilities to outcomes
    item[] <- rbinom(n*nvar, 1, item)                  
   result <- list(items=item,tau=d,theta=theta)
	return (result) 
	}  
	
	
"sim.irt" <- 
function (nvar = 5 ,n = 500, low=-3,high=3,a=NULL,c=0,z=1,d=NULL,mod="logistic") 
	{ 
	if(mod=="logistic") {result <- sim.npl(nvar,n,low,high,a,c,z,d)} else {result <- sim.npn(nvar,n,low,high,a,c,z,d)}
	return (result) 
	}  

"sim.npn" <- 
function (nvar = 5 ,n = 500, low=-3,high=3,a=NULL,c=0,z=1,d=NULL) 
	{ 
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} 
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n) # the latent variable
	
	item <- matrix(c+(z-c)*pnorm(a*(theta %+% t( -d))),n,nvar)
	#now convert these probabilities to outcomes
	
    item[] <- rbinom(n*nvar, 1, item)
	              
   result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=d,theta=theta)
	return (result) 
	}  
	
	
"sim.npl" <- 
function (nvar = 5 ,n = 500, low=-3,high=3,a=NULL,c=0,z=1,d=NULL) 
	{ 
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} 
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n)
	item <- matrix(c+(z-c)/(1+exp(a*(-theta %+% t( d)))),n,nvar)
    item[] <- rbinom(n*nvar, 1, item) #now convert these probabilities to outcomes	          
    result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=d,theta=theta)
	return (result) 
	} 