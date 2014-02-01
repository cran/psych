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
	colnames(item) <- paste("V",1:nvar,sep="")
	return (item) 
	}  

"sim.spherical" <-
function (simple=FALSE, nx=7,ny=12 ,nsub = 500,  xloading =.55, yloading = .55, zloading=.55, gloading=0, xbias=0,  ybias = 0, zbias=0,categorical=FALSE, low=-3,high=3,truncate=FALSE,cutpoint=0) 
	{ 

	nvar <- nx * (ny-1)
	errorweight <- sqrt(1-(xloading^2  + yloading^2 + zloading^2+  gloading^2))  #squared errors and true score weights add to 1
    g <- rnorm(nsub) 
	truex <- rnorm(nsub)  +xbias #generate normal true scores for x + xbias
	truey <- rnorm(nsub)  + ybias #generate normal true scores for y + ybias
	truez <- rnorm(nsub)  + zbias #generate normal true scores for y + ybias
    true <- matrix(c(g,truex,truey,truez),nsub)
    
	 #make a vector of radians (the whole way around the circle) if circumplex
	if(!simple) {f1  <- rep(cos(seq(0,pi,length.out = nx)),each=ny-1)*xloading
	 f2 <- rep(cos(seq(0,2*pi*(ny-1)/ny,length.out=ny-1)),nx)*yloading
    f3 <- rep(sin(seq(0,2*pi*(ny-1)/ny,length.out=ny-1)),nx)*zloading
    f2 <- f2 * (1-f1^2)  #added to make the cylinder a sphere
    f3 <- f3 * (1-f1^2)
    f <- matrix(c(rep(gloading,(ny-1)*nx),f1,f2,f3),(ny-1)*nx)
      
     } else {    #simple structure -- 
     nvar <- 4*ny
     f1  <- rep(c(1,0,-1,0),each=ny)
	 f2 <- rep(c(1,0,-1,0),ny)
     f3 <- rep(c(0,1,0,-1),ny)
    f <- matrix(c(rep(gloading,ny*4),f1,f2,f3),4*ny)
     }
        
	error<- matrix(rnorm(nsub*(nvar)),nsub)    #create normal error scores

	#true score matrix for each item reflects structure in radians
	

	item <- true %*% t(f) + errorweight*error   #observed item = true score + error score 
    if (categorical) {
        
    	item = round(item)       #round all items to nearest integer value
		item[(item<= low)] <- low     
		item[(item>high) ] <- high   
		}
	if (truncate) {item[item < cutpoint] <- 0  }
	colnames(item) <- paste("V",1:nvar,sep="")
	return (item) 
	}  



"sim.rasch" <-
function (nvar = 5 , n = 500, low=-3,high=3,d=NULL, a=1,mu=0,sd=1) 
	{ 
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))}
	theta <- rnorm(n,mu,sd)
	item <- matrix(t(1/(1+exp(a*t(-theta %+% t( d))))),n,nvar)
	#now convert these probabilities to outcomes
    item[] <- rbinom(n*nvar, 1, item)  
   colnames(item) <- paste("V",1:nvar,sep="")
   result <- list(items=item,tau=d,theta=theta)
	return (result) 
	}  
	
	
"sim.irt" <- 
function (nvar = 5 ,n = 500,low=-3,high=3,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,mod="logistic") 
	{ 
	if(mod=="logistic") {result <- sim.npl(nvar,n,low,high,a,c,z,d,mu,sd)} else {result <- sim.npn(nvar,n,low,high,a,c,z,d,mu,sd)}
	return (result) 
	}  

"sim.npn" <- 
function (nvar = 5 ,n = 500, low=-3,high=3,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1) 
	{ 
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} else {if(length(d)==1) d <- rep(d,nvar)}
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n,mu,sd) # the latent variable
	
	item <- matrix(t(c+(z-c)*pnorm(a*t(theta %+% t(- d)))),n,nvar)  #need to transpose and retranpose to get it right
	#now convert these probabilities to outcomes
	
    item[] <- rbinom(n*nvar, 1, item)
	colnames(item) <- paste("V",1:nvar,sep="")        
   result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=z,theta=theta)
	return (result) 
	}  
	
	
"sim.npl" <- 
function (nvar = 5 ,n = 500, low=-3,high=3,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1) 
	{ 
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} else {if(length(d)==1) d <- rep(d,nvar)}
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n,mu,sd)
	item <- matrix(t(c+(z-c)/(1+exp(a*t((-theta %+% t( d)))))),n,nvar)
    item[] <- rbinom(n*nvar, 1, item) #now convert these probabilities to outcomes	  
    colnames(item) <- paste("V",1:nvar,sep="")
    result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=z,theta=theta)
	return (result) 
	} 


"sim.poly" <- 
function (nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5,mod="logistic") 
	{ 
	if(mod=="normal") {result <- sim.poly.npn(nvar,n,low,high,a,c,z,d,mu,sd,cat)} else {result <- sim.poly.npl(nvar,n,low,high,a,c,z,d,mu,sd,cat)}
	return (result) 
	}  	
	
"sim.poly.npn" <- 
function (nvar = 5 ,n = 500, low=-2,high=2,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1,cat=5) 
	{ cat <- cat - 1
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} else {if(length(d)==1) d <- rep(d,nvar)}
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n,mu,sd) # the latent variable
	
	item <- matrix(t(c+(z-c)*pnorm(a*t(theta %+% t(- d)))),n,nvar)  #need to transpose and retranpose to get it right
	#now convert these probabilities to outcomes
	item[] <- rbinom(n*nvar, cat, item)
   
	colnames(item) <- paste("V",1:nvar,sep="")        
   result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=z,theta=theta)
	return (result) 
	}  	
	
"sim.poly.npl" <- 
function (nvar = 5 ,n = 500, low=-2,high=2,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1,cat=5) 
	{ cat <- cat - 1
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} else {if(length(d)==1) d <- rep(d,nvar)}
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n,mu,sd)
	item <- matrix(t(c+(z-c)/(1+exp(a*t((-theta %+% t( d)))))),n,nvar)
	item[] <- rbinom(n*nvar, cat, item)

      
    colnames(item) <- paste("V",1:nvar,sep="")
     
    result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=z,theta=theta)
	return (result) 
	} 
	
"sim.poly.ideal" <- 
function (nvar = 5 ,n = 500,low=-2,high=2,a=NULL,c=0,z=1,d=NULL, mu=0,sd=1,cat=5,mod="logistic") 
	{ 
	if(mod=="normal") {result <- sim.poly.ideal.npn(nvar,n,low,high,a,c,z,d,mu,sd,cat)} else {result <- sim.poly.ideal.npl(nvar,n,low,high,a,c,z,d,mu,sd,cat)}
	return (result) 
	}  	

	
"sim.poly.ideal.npl.absolute" <- 
function (nvar = 5 ,n = 500, low=-2,high=2,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1,cat=5) 
	{ cat <- cat -1  #binomial is based upon one fewer than categories
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} else {if(length(d)==1) d <- rep(d,nvar)}
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n,mu,sd)
	item <- matrix(t(c+(z-c)/(1+exp(a*t(abs(-theta %+% t( d)))))),n,nvar)
	item[] <- rbinom(n*nvar, cat, item)

      
    colnames(item) <- paste("V",1:nvar,sep="")
     
    result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=z,theta=theta)
	return (result) 
	} 
	
"sim.poly.ideal.npl" <- 
function (nvar = 5 ,n = 500, low=-2,high=2,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1,cat=5,theta=NULL) 
	{ cat <- cat -1  #binomial is based upon one fewer than categories
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} else {if(length(d)==1) d <- rep(d,nvar)}
	if(is.null(a)) {a <- rep(1,nvar)}
	if(is.null(theta)) {theta <- rnorm(n,mu,sd)}
	item <- 2*matrix((c+(z-c)*exp(a*t(-theta %+% t( d)))/(1+2*exp(a*t(-theta %+% t( d))) + exp(a*t(2*(-theta %+% t( d)))))),n,nvar,byrow=TRUE)
    p <- item
	item[] <- rbinom(n*nvar, cat, item)

      
    colnames(item) <- paste("V",1:nvar,sep="")
    result <- list(p=p,items=item,discrimination=a,difficulty=d,gamma=c,zeta=z,theta=theta)
	return (result) 
	} 
	


"sim.poly.ideal.npn" <- 
function (nvar = 5 ,n = 500, low=-2,high=2,a=NULL,c=0,z=1,d=NULL,mu=0,sd=1,cat=5) {
warning("Not ready for prime time")
	 cat <- cat -1  #binomial is based upon one fewer than categories
	if(is.null(d)) {d <- seq(low,high,(high-low)/(nvar-1))} else {if(length(d)==1) d <- rep(d,nvar)}
	if(is.null(a)) {a <- rep(1,nvar)}
	theta <- rnorm(n,mu,sd) # the latent variable
	
	item <- matrix(t(c+(z-c)*pnorm(abs(a*t(theta %+% t(- d))))),n,nvar)  #need to transpose and retranpose to get it right
	#now convert these probabilities to outcomes
	item[] <- rbinom(n*nvar, cat, item)
   
	colnames(item) <- paste("V",1:nvar,sep="")        
   result <- list(items=item,discrimination=a,difficulty=d,gamma=c,zeta=z,theta=theta)
	return (result) 
	}  	



#simulate a particular polychoric structure R with specified marginals (m)
#Modified December 17,2013 to actually work
"sim.poly.mat" <- function(R,m,n) {
 e <- eigen(R)
 v <- pmax(e$values,0)
 nvar <- ncol(R)
 ncat <- ncol(m) 
X <- matrix(rnorm(nvar*n),n)
X <- t(e$vectors %*% sqrt(diag(v)) %*% t(X))
marg <- m
Y <- matrix(0,ncol=n,nrow=nvar)
for(i in 1:(ncat)) {
Y[t(X) > marg[,i]] <- i }
return(t(Y))
}