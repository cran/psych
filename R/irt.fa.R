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
diff <- tau/sqrt(1-t$loadings[,1]^2)
discrim <- t$loadings/sqrt(1-t$loadings[,1]^2)
class(diff) <- NULL
class(discrim) <- NULL
tl <- t$loadings
class(tl) <- NULL
coeff <- data.frame(diff,discrim,tau,tl)
nlevels <- dim(diff)[2]
if(!is.null(nlevels)) {
colnames(coeff) <- c(paste("Location",1:nlevels,sep=""),"Discrimination",paste("tau",1:nlevels,sep=""),"Loading") } else {
colnames(coeff) <- c("Location","Discrimination","tau","Loading")}
result <- list(coefficients=coeff,stats = t,rho=r,tau=tau,n.obs=n.obs,Call=cl)
if(typ=="tet") {
class(result) <- c("psych","irt.fa")} else {class(result) <- c("psych","irt.poly")}
return(result)
}



"plot.irt" <- 
function(x,xlab,ylab,main,D,type=c("ICC","IIC","test"),...) {
item <- x
x <- NULL
nvar <- dim(item$coefficients)[1]
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
plot(x,logistic(x,a=item$coefficients[1,2]*D,d=item$coefficients[1,1]),ylim=c(0,1),ylab=ylab,xlab=xlab,type="l",main=main)
text(item$coefficients[1,1],.53,1)
for(i in 2:nvar) {
	lines(x,logistic(x,a=item$coefficients[i,2]*D,d=item$coefficients[i,1]),lty=c(1:6)[(i %% 6) + 1 ])
	text(item$coefficients[i,1],.53,i)
	}
	}  else {
	tInfo <- matrix(0,ncol=nvar,nrow=length(x))
	for(i in 1:nvar) tInfo[,i] <- logisticInfo(x,a=item$coefficients[i,2]*D,d=item$coefficients[i,1])
	
	if(type=="test") {
		if(missing(main)) main <- "Test information -- item parameters from factor analysis"
		testInfo <- rowSums(tInfo)
		if(missing(ylab)) ylab <- "Test Information"
		plot(x,testInfo,typ="l",ylim=c(0,max(testInfo)),ylab="Test Information",main=main)} else {
		if(missing(ylab)) ylab <- "Item Information"
	if(missing(main)) main <- "Item information from factor analysis"
	plot(x,logisticInfo(x,a=item$coefficients[1,2]*D,d=item$coefficients[1,1]),ylim=c(0,max(tInfo)+.03),ylab=ylab,xlab=xlab,type="l",main=main)
text(item$coefficients[1,1],max(tInfo[,1])+.03,1)
for(i in 2:nvar) {
	lines(x,logisticInfo(x,a=item$coefficients[i,2]*D,d=item$coefficients[i,1]),lty=c(1:6)[(i %% 6) + 1 ])
	text(item$coefficients[i,1],max(tInfo[,i])+.03,i)
	}
	} 
	}
	
}


"logisticInfo" <-  
function(x,d=0, a=1,c=0,z=1) {c + (z-c)*exp(a*(d-x))*a^2/(1+exp(a*(d-x)))^2}



"plot.poly" <- 
function(x,D,xlab,ylab,ylim,main,type=c("ICC","IIC","test"),...) {
item <- x
x <- NULL
if(missing(ylim)) ylim <- c(0,1)
nvar <- dim(item$coefficients)[1]
ncat <- (dim(item$coefficients)[2] -2)/2
if(missing(type)) {type = "IIC"}   

if(missing(D)) {D <- 1.702
if(missing(xlab)) xlab <- "Latent Trait (normal scale)"
x <- seq(-3,3,.1)}
difficulty <- item$coefficients[,1:ncat]
discrimination <- item$coefficients[ncat+1]
if(D==1) {if(missing(xlab)) xlab <- "Latent Trait (logistic scale)"}
if(missing(xlab)) xlab <- "Latent Trait"

if(is.null(x)) x <- seq(-4,4,.1)
if(type=="ICC") {
if(missing(main)) main <- "Item parameters from factor analysis"
if(missing(ylab)) ylab <- "Probability of Response"

for(i in 1:nvar) {
	if(discrimination[i,1] > 0 ) {
		plot(x,logistic(x,a=-D*discrimination[i,1],d=difficulty[i,1]),ylim=ylim,ylab=ylab,xlab=xlab,type="l",main=main,...) 
		text(0,.70,colnames(item$rho)[i])} else {
		plot(x,logistic(x,a=D*discrimination[i,1],d=difficulty[i,1]),ylim=ylim,ylab=ylab,xlab=xlab,type="l",main=main,...)
		text(max(0),.7,paste("-",colnames(item$rho)[i],sep=""))
			}
  for (j in 2:(ncat-1))  { 
		if(discrimination[i,1] > 0 ) {
			lines(x,(-logistic(x,a=D*discrimination[i,1],d=difficulty[i,j])+logistic(x,a=D*discrimination[i,1],d=difficulty[i,j-1])),lty=c(1:6)[(j %% 6) + 1 ])
			} else {lines(x,(-logistic(x,a=-D*discrimination[i,1],d=difficulty[i,j])+logistic(x,a=-D*discrimination[i,1],d=difficulty[i,j-1])),lty=c(1:6)[(j %% 6) + 1 ])}
			}
	if(discrimination[i,1] > 0 ) {
	lines(x,(logistic(x,a=D*discrimination[i,1],d=difficulty[i,ncat])))
	} else {lines(x,(logistic(x,a=-D*discrimination[i,1],d=difficulty[i,ncat]))) }
	}

	}  else {
	
	#item and test information
	x <- as.matrix(x,ncol=1)
	discrimination <- discrimination[,1]  #convert to a vector so we can do the next step
	tInfo <- apply(x,1,logisticInfo,a=discrimination,d=difficulty)  #fix this to include items correctly
	tInfo <- array(unlist(tInfo),dim=c(nvar,ncat,length(x)))  #this is now an array with items levels and x 

	testInfol <- list()
	for (xi in 1:length(x)) testInfol[[xi]] <- rowSums(tInfo[,,xi])
	testInfo <- matrix(unlist(testInfol),ncol=nvar,byrow=TRUE)
	
	if(type=="test") { 
	if(missing(main)) main <- "Test information -- item parameters from factor analysis"
	
	
	if(missing(ylab)) ylab <- "Test Information"
	rsInfo <- rowSums(testInfo)
	 plot(x,rsInfo,typ="l",ylim=c(0,max(rsInfo)),ylab=ylab,xlab=xlab,main=main)
	 } else { 
	if(missing(ylab)) ylab <- "Item Information"
	if(missing(main)) main <- "Item information from factor analysis"
	plot(x,testInfo[,1],ylim=c(0,max(testInfo,na.rm=TRUE)+.03),ylab=ylab,xlab=xlab,type="l",main=main)
	#xmax <- which
if(discrimination[1] > 0 ) {text(x[which.max(testInfo[,1])],max(testInfo[,1])+.03,colnames(item$rho)[1])} else {text(x[which.max(testInfo[,1])],max(testInfo[,1])+.03,paste("-",colnames(item$rho)[1],sep=""))}
for(i in 2:nvar) {
	lines(x,testInfo[,i],lty=c(1:6)[(i %% 6) + 1 ])
	if(discrimination[i] > 0 ) {
	text(x[which.max(testInfo[,i])],max(testInfo[,i])+.03,colnames(item$rho)[i]) } else {text(x[which.max(testInfo[,i])],max(testInfo[,i])+.03,paste("-",colnames(item$rho)[i],sep=""))}
	}
	} 
	}
	
}

