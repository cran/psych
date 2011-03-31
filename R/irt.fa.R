"irt.fa" <- 
function(x,...) {
cl <- match.call()
if (is.matrix(x) | is.data.frame(x)) {
	n.obs <- dim(x)[1]
	tx <- table(as.matrix(x))
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x)
	    typ = "tet"} else {tet <- polychoric(x)
	    typ = "poly"}

	r <- tet$rho
	tau <- tet$tau}  else {if (!is.null(x$rho)) { r <- x$rho
   			tau <- x$tau
   			n.obs <- x$n.obs
   			typ <- class(x)[2]
   			if (typ == "irt.fa") typ <- "tet"
   			  }  else {stop("x must  be a data.frame or matrix or the result from tetra or polychoric")}
              }
t <- fa(r,n.obs=n.obs,...)
nf <- dim(t$loadings)[2]
 diffi <- list() 
     for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-t$loadings[,i]^2)
     }
     
discrim <- t$loadings/sqrt(1-t$loadings^2)
class(diffi) <- NULL
class(discrim) <- NULL
tl <- t$loadings
class(tl) <- NULL
irt <- list(difficulty=diffi,discrimination=discrim)
nlevels <- dim(diffi[[1]])[2]
#if(!is.null(nlevels)) {
#colnames(coeff) <- c(paste("Location",1:nlevels,sep=""),"Discrimination",paste("tau",1:nlevels,sep=""),"Loading") } else {
#colnames(coeff) <- c("Location","Discrimination","tau","Loading")}
result <- list(irt=irt,fa = t,rho=r,tau=tau,n.obs=n.obs,Call=cl)
if(typ =="tet") {
class(result) <- c("psych","irt.fa")} else {class(result) <- c("psych","irt.poly")}
return(result)
}



"plot.irt" <- 
function(x,xlab,ylab,main,D,type=c("ICC","IIC","test"),cut=.3,labels,...) {

item <- x
if((is.data.frame(x)) | (is.matrix(x))) {nf <- dim(x)[2] -1} else {
           
nf <- length(x$irt$difficulty)}
#if there was more than 1 factor, repeat the figure nf times
for(i in 1:nf) {if((is.data.frame(item)) | (is.matrix(item))) {discrimination <- item[,1]
     location <- item[,i+1]
    } else {
  discrimination=item$irt$discrimination[,i]
  location=item$irt$difficulty[[i]] }
x <- NULL 
nvar <- length(discrimination)
if(missing(labels) ) labels = 1:nvar
if(missing(type)) {type = "IIC"}   

if(missing(D)) {D <- 1.702
if(missing(xlab)) xlab <- "Latent Trait (normal scale)"
x <- seq(-3,3,.1)}

if(D==1) {if(missing(xlab)) xlab <- "Latent Trait (logistic scale)"}
if(missing(xlab)) xlab <- "Latent Trait"

if(is.null(x)) x <- seq(-4,4,.1)
if(type=="ICC") {
if(missing(main)) main <- "Item parameters from factor analysis"
if(missing(ylab)) ylab <- "Probability of Response"
ii <- 1 
while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
plot(x,logistic(x,a=discrimination[ii]*D,d=location[ii]),ylim=c(0,1),ylab=ylab,xlab=xlab,type="l",main=main)
text(location[ii],.53,labels[ii])
for(i in (ii+1):nvar) {
   if(abs(discrimination[i])  > cut) {
	lines(x,logistic(x,a=discrimination[i]*D,d=location[i]),lty=c(1:6)[(i %% 6) + 1 ])
	text(location[i],.53,labels[i])}
	}
	}  else {
	tInfo <- matrix(0,ncol=nvar,nrow=length(x))
	for(i in 1:nvar) {
	   if(abs(discrimination[i])  > cut) {
	tInfo[,i] <- logisticInfo(x,a=discrimination[i]*D,d=location[i])} else {tInfo[,i] <- 0}
	}
	if(type=="test") {
		if(missing(main)) main <- "Test information -- item parameters from factor analysis"
		testInfo <- rowSums(tInfo)
		if(missing(ylab)) ylab <- "Test Information"
		plot(x,testInfo,typ="l",ylim=c(0,max(testInfo)),ylab="Test Information",main=main)} else {
		if(missing(ylab)) ylab <- "Item Information"
	if(missing(main)) main <- "Item information from factor analysis"
	ii <- 1 
while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
	plot(x,logisticInfo(x,a=discrimination[ii]*D,d=location[ii]),ylim=c(0,max(tInfo)+.03),ylab=ylab,xlab=xlab,type="l",main=main)
text(location[ii],max(tInfo[,ii])+.03,labels[1])
for(i in (ii+1):nvar) {
    if(abs(discrimination[i])  > cut) {
	lines(x,logisticInfo(x,a=discrimination[i]*D,d=location[i]),lty=c(1:6)[(i %% 6) + 1 ])
	text(location[i],max(tInfo[,i])+.03,labels[i])
	}}} }
	devAskNewPage(ask = TRUE)}	
	devAskNewPage(ask = FALSE)
}


"logisticInfo" <-  
function(x,d=0, a=1,c=0,z=1) {c + (z-c)*exp(a*(d-x))*a^2/(1+exp(a*(d-x)))^2}



"irt.fa" <- 
function(x,...) {
cl <- match.call()
if (is.matrix(x) | is.data.frame(x)) {
	n.obs <- dim(x)[1]
	tx <- table(as.matrix(x))
	if(dim(tx)[1] ==2) {tet <- tetrachoric(x)
	typ = "tet"} else {tet <- polychoric(x)
	    typ = "poly"}

	r <- tet$rho
	tau <- tet$tau}  else {if (!is.null(x$rho)) { r <- x$rho
   			tau <- x$tau
   			n.obs <- x$n.obs
   			typ <- class(x)[2]
   			  }  else {stop("x must  be a data.frame or matrix or the result from tetra or polychoric")}
              }
t <- fa(r,n.obs=n.obs,...)
nf <- dim(t$loadings)[2]
 diffi <- list() 
     for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-t$loadings[,i]^2)
     }
     
discrim <- t$loadings/sqrt(1-t$loadings^2)
class(diffi) <- NULL
class(discrim) <- NULL
tl <- t$loadings
class(tl) <- NULL
irt <- list(difficulty=diffi,discrimination=discrim)
nlevels <- dim(diffi[[1]])[2]
#if(!is.null(nlevels)) {
#colnames(coeff) <- c(paste("Location",1:nlevels,sep=""),"Discrimination",paste("tau",1:nlevels,sep=""),"Loading") } else {
#colnames(coeff) <- c("Location","Discrimination","tau","Loading")}
result <- list(irt=irt,fa = t,rho=r,tau=tau,n.obs=n.obs,Call=cl)
if(typ=="tet") {
class(result) <- c("psych","irt.fa")} else {class(result) <- c("psych","irt.poly")}
return(result)
}



"plot.irt" <- 
function(x,xlab,ylab,main,D,type=c("ICC","IIC","test"),cut=.3,labels=NULL,...) {

item <- x
if((is.data.frame(x)) | (is.matrix(x))) {nf <- dim(x)[2] -1} else {
           
nf <- length(x$irt$difficulty)}
#if there was more than 1 factor, repeat the figure nf times
for(i in 1:nf) {if((is.data.frame(item)) | (is.matrix(item))) {discrimination <- item[,1]
     location <- item[,i+1]
    } else {
  discrimination=item$irt$discrimination[,i]
  location=item$irt$difficulty[[i]] }
x <- NULL 
nvar <- length(discrimination)
if(is.null(labels) ) labels = 1:nvar
if(missing(type)) {type = "IIC"}   

if(missing(D)) {D <- 1.702
if(missing(xlab)) xlab <- "Latent Trait (normal scale)"
x <- seq(-3,3,.1)}

if(D==1) {if(missing(xlab)) xlab <- "Latent Trait (logistic scale)"}
if(missing(xlab)) xlab <- "Latent Trait"

if(is.null(x)) x <- seq(-4,4,.1)
if(type=="ICC") {
if(missing(main)) main <- "Item parameters from factor analysis"
if(missing(ylab)) ylab <- "Probability of Response"
ii <- 1 
while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
plot(x,logistic(x,a=discrimination[ii]*D,d=location[ii]),ylim=c(0,1),ylab=ylab,xlab=xlab,type="l",main=main)
text(location[ii],.53,labels[ii])
for(i in (ii+1):nvar) {
   if(abs(discrimination[i])  > cut) {
	lines(x,logistic(x,a=discrimination[i]*D,d=location[i]),lty=c(1:6)[(i %% 6) + 1 ])
	text(location[i],.53,labels[i])}
	}
	}  else {
	tInfo <- matrix(0,ncol=nvar,nrow=length(x))
	for(i in 1:nvar) {
	   if(abs(discrimination[i])  > cut) {
	tInfo[,i] <- logisticInfo(x,a=discrimination[i]*D,d=location[i])} else {tInfo[,i] <- 0}
	}
	if(type=="test") {
		if(missing(main)) main <- "Test information -- item parameters from factor analysis"
		testInfo <- rowSums(tInfo)
		if(missing(ylab)) ylab <- "Test Information"
		plot(x,testInfo,typ="l",ylim=c(0,max(testInfo)),ylab="Test Information",main=main)} else {
		if(missing(ylab)) ylab <- "Item Information"
	if(missing(main)) main <- "Item information from factor analysis"
	ii <- 1 
while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
	plot(x,logisticInfo(x,a=discrimination[ii]*D,d=location[ii]),ylim=c(0,max(tInfo)+.03),ylab=ylab,xlab=xlab,type="l",main=main)
text(location[ii],max(tInfo[,ii])+.03,labels[1])
for(i in (ii+1):nvar) {
    if(abs(discrimination[i])  > cut) {
	lines(x,logisticInfo(x,a=discrimination[i]*D,d=location[i]),lty=c(1:6)[(i %% 6) + 1 ])
	text(location[i],max(tInfo[,i])+.02,labels[i])
	}}} }
	devAskNewPage(ask = TRUE)}	
	devAskNewPage(ask = FALSE)
}


"logisticInfo" <-  
function(x,d=0, a=1,c=0,z=1) {c + (z-c)*exp(a*(d-x))*a^2/(1+exp(a*(d-x)))^2}



"plot.poly" <- 
function(x,D,xlab,ylab,ylim,main,type=c("ICC","IIC","test"),cut=.3,labels,...) {
item <- x
if((is.data.frame(x)) | (is.matrix(x))) {nf <- dim(x)[2] -1} else {
           
nf <- length(x$irt$difficulty)}
temp <- list()
#if there was more than 1 factor, repeat the figure nf times
  
x <- NULL
if(missing(ylim)) ylim <- c(0,1)
nvar <- length(item$irt$discrimination[,1])
ncat <- dim(item$irt$difficulty[[1]])[2]
if(missing(type)) {type = "IIC"}   

if(missing(D)) {D <- 1.702
if(missing(xlab)) xlab <- "Latent Trait (normal scale)"
x <- seq(-3,3,.1)}


if(D==1) {if(missing(xlab)) xlab <- "Latent Trait (logistic scale)"}
if(missing(xlab)) xlab <- "Latent Trait"

if(is.null(x)) x <- seq(-4,4,.1)

for(f in 1:nf) {discrimination=item$irt$discrimination[,f]
  location=item$irt$difficulty[[f]]
  difficulty <- location[,1:ncat]
  
if(type=="ICC") {
if(missing(main)) main <- "Item parameters from factor analysis"
if(missing(ylab)) ylab <- "Probability of Response"

for(i in 1:nvar) {
 if (abs(discrimination[i]) > cut) {
	if(discrimination[i] > 0 ) {
		plot(x,logistic(x,a=-D*discrimination[i],d=location[i,1]),ylim=ylim,ylab=ylab,xlab=xlab,type="l",main=main,...) 
		text(0,.70,colnames(item$rho)[i])} else { 
		plot(x,logistic(x,a=-D*discrimination[i],d=-location[i,1]),ylim=ylim,ylab=ylab,xlab=xlab,type="l",main=main,...)
		text(max(0),.7,paste("-",colnames(item$rho)[i],sep=""))
			}
  for (j in 2:(ncat-1))  { 
		if(discrimination[i] > 0 ) {
			lines(x,(-logistic(x,a=D*discrimination[i],d=location[i,j])+logistic(x,a=D*discrimination[i],d=location[i,j-1])),lty=c(1:6)[(j %% 6) + 1 ])
			} else {lines(x,(-logistic(x,a=-D*discrimination[i],d= -location[i,j])+logistic(x,a=-D*discrimination[i],d=-location[i,j-1])),lty=c(1:6)[(j %% 6) + 1 ])}
			}
	if(discrimination[i] > 0 ) {
	lines(x,(logistic(x,a=D*discrimination[i],d=location[i,ncat])))
	} else {lines(x,(logistic(x,a=-D*discrimination[i],d=-location[i,ncat]))) }
	}}

	}  else {
	
	#item and test information
	x <- as.matrix(x,ncol=1)
	
	tInfo <- apply(x,1,logisticInfo,a=discrimination,d=sign(discrimination) * difficulty)  
	tInfo <- array(unlist(tInfo),dim=c(nvar,ncat,length(x)))  #this is now an array with items levels and x 

	testInfo <- matrix(NA,ncol=nvar,nrow=length(x))
	for (xi in 1:length(x)) {
	for (i in 1:nvar) {

	if (abs(discrimination[i]) > cut) {
	testInfo[[xi,i]] <- sum(tInfo[i,,xi]) } else {testInfo[[xi,i]] <- 0}}
	}
	
	
	if(type=="test") { 

	
	if(missing(main)) main <- "Test information for factor "
	main1 <- paste(main,'  ',f)
	if(missing(ylab)) ylab <- "Test Information"
	rsInfo <- rowSums(testInfo)
	 plot(x,rsInfo,typ="l",ylim=c(0,max(rsInfo)),ylab=ylab,xlab=xlab,main=main1)
	 } else { 
	if(missing(ylab)) ylab <- "Item Information"
	if(missing(main)) main <- "Item information from factor analysis"
	
	ii <- 1 
while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
	plot(x,testInfo[,ii],ylim=c(0,max(testInfo,na.rm=TRUE)+.03),ylab=ylab,xlab=xlab,type="l",main=main)
	#xmax <- which
if(discrimination[ii] > 0 ) {text(x[which.max(testInfo[,ii])],max(testInfo[,ii])+.03,colnames(item$rho)[ii])} else {text(x[which.max(testInfo[,ii])],max(testInfo[,ii])+.03,paste("-",colnames(item$rho)[ii],sep=""))}
for(i in (ii+1):nvar) { if (abs(discrimination[i]) > cut) {
	lines(x,testInfo[,i],lty=c(1:6)[(i %% 6) + 1 ])
	if(discrimination[i] > 0 ) {
	text(x[which.max(testInfo[,i])],max(testInfo[,i])+.03,colnames(item$rho)[i]) } else {text(x[which.max(testInfo[,i])],max(testInfo[,i])+.03,paste("-",colnames(item$rho)[i],sep=""))}
	}}
	} }
	if (type !="ICC")  {temp[[f]] <- testInfo}
	}
	
	
	#This next piece of code is being developed and currently does nothing
 if (type !="ICC")  {result <- matrix(NA,nvar,nf)
 
    	for (f in 1:nf) {
     		 for (j in 1:nvar) {
    			result[i,f] <- NA  
    				}
    				}
     class(result) <- c("psych","polyinfo")
     return(result)}
     
}


"irt.select" <- function(x,y) {
  if(is.null(dim(x$tau))) {typ="tet"} else {typ="poly"}
  rho <- x$rho[y,y]
  tau <- x$tau[y]
  n.obs <- x$n.obs
  result <- list(rho=rho,tau=tau,n.obs=n.obs)
  class(result) <- c("psych",typ)
return(result)
}