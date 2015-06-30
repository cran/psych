"plot.irt" <- 
function(x,xlab,ylab,main,D,type=c("ICC","IIC","test"),cut=.3,labels=NULL,keys=NULL,xlim,ylim,y2lab,lncol="black",...) {
if(class(x)[2] == "irt.poly") {
if(missing(type)) type = "IIC" 

plot.poly(x=x,D=D,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main,type=type,cut=cut,labels=labels,keys=keys,y2lab=y2lab,lncol=lncol,...)} else {

item <- x
temp <- list()
sumtemp <- list()
byKeys <- FALSE
if((is.data.frame(x)) | (is.matrix(x))) {nf <- dim(x)[2] -1} else {
           
nf <- length(x$irt$difficulty)}
#if there was more than 1 factor, repeat the figure nf times
#or, if there is one factor but multiple keys
if(!is.null(keys)) {
   nkeys = ncol(keys)
   if (nf < nkeys) {byKeys <- TRUE
                nf <- nkeys}
   }
for(f in 1:nf) {if((is.data.frame(item)) | (is.matrix(item))) {if(byKeys) {discrimination <- item[,1]} else {discrimination <- item[,f]}
     if(byKeys) {location <- item[,2]} else {location <- item[,f+1]}
    } else {
  if(byKeys) {discrimination=item$irt$discrimination[,1]} else {discrimination=item$irt$discrimination[,f]}
  if(!is.null(keys)) discrimination <- discrimination *abs( keys[,f])
  if(byKeys) {location=item$irt$difficulty[[1]]} else {location=item$irt$difficulty[[f]] }}
x <- NULL 
nvar <- length(discrimination)
if(is.null(labels)) {if(!is.null(rownames(item$irt$discrimination)))  {labels = rownames(item$irt$discrimination)} else {labels <- 1:nvar}}
if(missing(type)) {type = "IIC"}   

if(missing(D)) {D <- 1.702
if(missing(xlab)) xlab <- "Latent Trait (normal scale)"
x <- seq(-3,3,.1)
summaryx <- seq(-3,3,1)}

if(D==1) {if(missing(xlab)) xlab <- "Latent Trait (logistic scale)"}
if(missing(xlab)) xlab <- "Latent Trait"

if(is.null(x))  {x <- seq(-4,4,.1)
   summaryx <- seq(-3,3,1)}
lenx <- length(x)
sumInfo <- matrix(NA,ncol=nvar,nrow=length(summaryx))
summaryx <- as.matrix(summaryx,ncol=1)
if(type=="ICC") {
 		summtInfo <- NULL 
		if(missing(main)) main <- "Item parameters from factor analysis"
		if(missing(ylab)) ylab <- "Probability of Response"
		if(length(lncol)<2) lncol <- rep(lncol,nvar)
		ii <- 1 
		while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
		plot(x,logistic(x,a=discrimination[ii]*D,d=location[ii]),ylim=c(0,1),ylab=ylab,xlab=xlab,type="l",main=main,col=lncol[1],...)
		text(location[ii],.53,labels[ii])
		for(i in (ii+1):nvar) {
  		 if(abs(discrimination[i])  > cut) {
			lines(x,logistic(x,a=discrimination[i]*D,d=location[i]),lty=c(1:6)[(i %% 6) + 1 ],col=lncol[i],...)
			text(location[i],.53,labels[i])}
		}
	}  else {  #not ICC
	tInfo <- matrix(0,ncol=nvar,nrow=length(x))
	for(i in 1:nvar) {
	   if(abs(discrimination[i])  > cut) {
		tInfo[,i] <- logisticInfo(x,a=discrimination[i]*D,d=location[i])
		sumInfo[,i] <- logisticInfo(summaryx,a=discrimination[i]*D,d=location[i]) } else {tInfo[,i] <- 0
		sumInfo[,i] <- 0}
		}
	AUC <- colSums(tInfo)
	max.info <- apply(tInfo,2,which.max)
	if(type=="test") {
		if(missing(main)) main <- "Test information -- item parameters from factor analysis"
		if(missing(y2lab)) y2lab <- "Reliability"
		testInfo <- rowSums(tInfo)
		if(missing(ylab)) ylab <- "Test Information"
		if(missing(xlab)) xlab <- "Latent Trait (normal scale)"
		if(length(lncol)< 2) lncol <- rep(lncol,nvar)
		op <- par(mar=c(5,4,4,4))  #set the margins a bit wider
		plot(x,testInfo,typ="l",ylim=c(0,max(testInfo)),ylab=ylab,xlab=xlab,main=main,col=lncol[1],...)
		 ax4 <- seq(0,max(testInfo),max(testInfo)/4)
	 rel4 <- round(1-1/ax4,2)
	 rel4[1] <- NA
	 axis(4,at=ax4,rel4) 
	 	
	 mtext(y2lab,side=4,line=2)
	 op <- par(op)   #set them back to what we had before
		} else {
		if(missing(ylab)) ylab <- "Item Information"
	if(missing(main)) main <- "Item information from factor analysis"
	if(length(lncol) <2) lncol <- rep(lncol,nvar)
	ii <- 1 
while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
    if(missing(ylim)) {ylimit=c(0,max(tInfo)+.03)} else {ylimit <- ylim}
	plot(x,logisticInfo(x,a=discrimination[ii]*D,d=location[ii]),ylim=ylimit,ylab=ylab,xlab=xlab,type="l",main=main,col=lncol[1],...)
text(location[ii],max(tInfo[,ii])+.03,labels[ii])

for(i in (ii+1):nvar) {
    if(abs(discrimination[i])  > cut) {
	lines(x,logisticInfo(x,a=discrimination[i]*D,d=location[i]),lty=c(1:6)[(i %% 6) + 1 ],col=lncol[i])
	text(location[i],max(tInfo[,i])+.02,labels[i])
	}}} 
	if  (type !="ICC") {temp[[f]] <- list(AUC=AUC,max.info=max.info)
	                   sumInfo <- t(sumInfo)
	                   colnames(sumInfo) <- summaryx
	                   rownames(sumInfo) <- rownames(item$rho)
	                  sumtemp[[f]] <- sumInfo}
	}
	devAskNewPage(ask = TRUE)}   #end of f loop	
	devAskNewPage(ask = FALSE)
	if(type!="ICC") {
     AUC <- matrix(NA,ncol=nf,nrow=nvar)
     max.info <- matrix(NA,ncol=nf,nrow=nvar) 
       for(f in 1:nf) {
     AUC[,f] <- temp[[f]]$AUC

	max.info[,f] <- temp[[f]]$max.info}
	AUC <- AUC/lenx  #quasi normalize it 
	max.info <- (max.info - lenx/2)*6/(lenx-1)
	max.info[max.info < -2.9] <- NA
	if(byKeys) {colnames(AUC) <- colnames(max.info) <- colnames(keys)} else {colnames(AUC) <- colnames(max.info) <- colnames(item$irt$discrimination)}
    rownames(AUC) <- rownames(max.info) <- rownames(item$rho)
    
   
    
    
    
    result <- list(AUC=AUC,max.info=max.info,sumInfo =sumtemp)
	invisible(result)				
     class(result) <- c("psych","polyinfo")
   invisible(result)}
}	
}


"logisticInfo" <-  
function(x,d=0, a=1,c=0,z=1) {c + (z-c)*exp(a*(d-x))*a^2/(1+exp(a*(d-x)))^2}



"plot.poly" <- 
function(x,D,xlab,ylab,xlim,ylim,main,type=c("ICC","IIC","test"),cut=.3,labels=NULL,keys=NULL,y2lab,lncol="black",...) {

item <- x
byKeys <- FALSE
if((is.data.frame(x)) | (is.matrix(x))) {nf <- dim(x)[2] -1} else {
           
nf <- length(x$irt$difficulty)}
#if there was more than 1 factor, repeat the figure nf times
#or, if there is one factor but multiple keys
if(!is.null(keys)) {
   nkeys = ncol(keys)
   if (nf < nkeys) {byKeys <- TRUE
                nf <- nkeys}
                }

temp <- list()
sumtemp <- list()
 
x <- NULL

nvar <- length(item$irt$discrimination[,1])
ncat <- dim(item$irt$difficulty[[1]])[2]
if(length(lncol) < 2) lncol <- rep(lncol,nvar)
if(missing(type)) {type = "IIC"}   

if(missing(D)) {D <- 1.702
if(missing(xlab)) xlab <- "Latent Trait (normal scale)"
if(missing(xlim)) {x <- seq(-3,3,.1) } else {x <- seq(xlim[1],xlim[2],.1)} #used for item summary table
summaryx <- x 
}		 
if(D==1) {if(missing(xlab)) xlab <- "Latent Trait (logistic scale)"}
if(missing(xlab)) xlab <- "Latent Trait"

if(is.null(x)) {x <- seq(-4,4,.1)
              summaryx <- seq(-3,3,1)} #used for item summary table
if(is.null(labels)) {if(!is.null(rownames(item$irt$discrimination)))  {labels = rownames(item$irt$discrimination)} else {labels <- 1:nvar}}

lenx <- length(x)
sumInfo <- matrix(NA,ncol=nvar,nrow=length(summaryx))

#if there was more than 1 factor, repeat the figure nf times
for(f in 1:nf) {if(byKeys) {discrimination <- item$irt$discrimination[,1]} else {discrimination <-item$irt$discrimination[,f]}
 if(!is.null(keys)) discrimination <- discrimination * abs(keys[,f])
 if(any(is.nan(discrimination))) {bad <- which(is.nan(discrimination))
        discrimination[is.nan(discrimination)] <- max(discrimination,na.rm=TRUE) 
     warning("An discrimination  with a NaN value was replaced with the maximum discrimination for factor = ",f, " and item ",labels[bad], "\nexamine the factor analysis object (fa)  to identify the Heywood case") }
     if(byKeys) {location=item$irt$difficulty[[1]]} else { location=item$irt$difficulty[[f]]}
  difficulty <- location[,1:ncat]
  

if(type=="ICC") { #this draws the item characteristic curves
  #summtInfo <- NULL 

if(missing(main)) {main <- "Item parameters from factor analysis"
                  main1 <- paste(main,'  ',f)} else {if ( length(main) > 1) main1 <- main[f]}
if(missing(ylab)) ylab <- "Probability of Response"
if(missing(ylim)) ylim <- c(0,1)
	
for(i in 1:nvar) {
 if (abs(discrimination[i]) > cut) {
	if(discrimination[i] > 0 ) {
		plot(x,logistic(x,a=-D*discrimination[i],d=location[i,1]),ylim=ylim,ylab=ylab,xlab=xlab,type="l",main=main1,col=lncol[1],...) 
		text(0,.70,labels[i])} else { 
		plot(x,logistic(x,a=D*discrimination[i],d=location[i,1]),ylim=ylim,ylab=ylab,xlab=xlab,type="l",main=main1,col=lncol[1],...)
		text(max(0),.7,paste("-",labels[i],sep=""))
			}
  for (j in 2:(ncat))  {
		if(discrimination[i] > 0 ) {                              
			lines(x,(-logistic(x,a=D*discrimination[i],d=location[i,j])+logistic(x,a=D*discrimination[i],d=location[i,j-1])),lty=c(1:6)[(j %% 6) + 1 ],col=lncol[i])
			} else {lines(x,(-logistic(x,a=-D*discrimination[i],d= location[i,j])+logistic(x,a=-D*discrimination[i],d=location[i,j-1])),lty=c(1:6)[(j %% 6) + 1 ],col=lncol[i])}
			}
	if(discrimination[i] > 0 ) {                                     
	lines(x,(logistic(x,a=D*discrimination[i],d=location[i,ncat])),col=lncol[i])
	} else {lines(x,(logistic(x,a=-D*discrimination[i],d=location[i,ncat])),col=lncol[i]) }
	}}
	}  #now do the summary stuff for all cases
	
	
   #summaryx <- as.matrix(summaryx,ncol=1)
   # summtInfo <- apply(summaryx,1,logisticInfo,a=discrimination,d=sign(discrimination) * difficulty)  #notice that we just need to add the logistics, not the differences 
   # summtInfo <- array(unlist(summtInfo),dim=c(nvar,ncat,length(summaryx)))  #this is now an array with items levels and summaryx 
	  
	
	#item and test information
	x <- as.matrix(x,ncol=1)
	summaryx <- as.matrix(summaryx,ncol=1)
	tInfo <- apply(x,1,logisticInfo,a=discrimination,d=sign(discrimination) * difficulty)  
	tInfo <- array(unlist(tInfo),dim=c(nvar,ncat,length(x)))  #this is now an array with items levels and x 
	
	
	
    summtInfo <- apply(summaryx,1,logisticInfo,a=discrimination,d=sign(discrimination) * difficulty) 
    summtInfo <- array(unlist(summtInfo),dim=c(nvar,ncat,length(summaryx)))  #this is now an array with items levels and summaryx 
    
    tInfo[is.nan(tInfo)] <- 0    #this gets around the problem of no values 
    summtInfo[is.nan(summtInfo)] <- 0    #this gets around the problem of no values 
	testInfo <- matrix(NA,ncol=nvar,nrow=length(x))
	sumInfo <- matrix(NA,ncol=nvar,nrow=length(summaryx))
	
	for (xi in 1:length(x)) {
	  for (i in 1:nvar) {
	     if (abs(discrimination[i]) > cut) {
	      testInfo[[xi,i]] <- sum(tInfo[i,,xi]) } else {testInfo[[xi,i]] <- 0 }}
	      }
	    for (xi in 1:length(summaryx)) {
	  for (i in 1:nvar) {
	     if (abs(discrimination[i]) > cut) {  
	    sumInfo[[xi,i]] <- sum(summtInfo[i,,xi]) } else {sumInfo[[xi,i]] <- 0}}
	    }
	

	if(type=="test") { 

	
	if(missing(main)) {main <- "Test information for factor "
	                  main1 <- paste(main,' ',f)} else {if (length(main) > 1) {main1 <- main[f]} else {main1 <- paste(main,'',f)}}
	if(missing(ylab)) ylab <- "Test Information"
	if(missing(y2lab)) y2lab <- "Reliability"
	
	rsInfo <- rowSums(testInfo)
	if(missing(ylim)) ylim = c(0,max(rsInfo))
	op <- par(mar=c(5,4,4,4))  #set the margins a bit wider
	 plot(x,rsInfo,typ="l",ylim=ylim,ylab=ylab,xlab=xlab,main=main1,col=lncol[1],...)
	 ax4 <- seq(0,ylim[2],ylim[2]/4)
	 rel4 <- round(1-1/ax4,2)
	 rel4[1] <- NA
	
	 axis(4,at=ax4,rel4)
	 	 mtext(y2lab,side=4,line=2)
	 op <- par(op)
	 } else { if(type != "ICC") {
	if(missing(ylab)) ylab <- "Item Information"
	if(missing(main)) {main1 <- paste("Item information from factor analysis for factor" ,f)} else {if (length(main) > 1) {main1 <- main[f]} else {main1 <- main}}
	if(missing(ylim)) ylim <- c(0,1)
	ii <- 1 
while((abs(discrimination[ii]) < cut) && (ii < nvar)) {ii <- ii + 1} 
	plot(x,testInfo[,ii],ylim=c(0,max(testInfo,na.rm=TRUE)+.03),ylab=ylab,xlab=xlab,type="l",main=main1,col=lncol[1],...)

if(discrimination[ii] > 0 ) {text(x[which.max(testInfo[,ii])],max(testInfo[,ii])+.03,labels[ii])} else {text(x[which.max(testInfo[,ii])],max(testInfo[,ii])+.03,paste("-",labels[ii],sep=""))}
for(i in (ii+1):nvar) { if (abs(discrimination[i]) > cut) {
	lines(x,testInfo[,i],lty=c(1:6)[(i %% 6) + 1 ],col=lncol[i])
	if(discrimination[i] > 0 ) {
	text(x[which.max(testInfo[,i])],max(testInfo[,i])+.03,labels[i]) } else {text(x[which.max(testInfo[,i])],max(testInfo[,i])+.03,paste("-",labels[i],sep=""))}
	}}
	} }
	
	#if (type !="ICC")  {
	temp[[f]] <- testInfo
	 
	 sumInfo <- t(sumInfo) 
	 rownames(sumInfo) <- labels
	 colnames(sumInfo) <- summaryx
	 
	sumtemp[[f]] <- sumInfo    #keep the summary information for each factor

	}  #end of 1:nf loop
	
     AUC <- matrix(NaN,ncol=nf,nrow=nvar)
     max.info <- matrix(NaN,ncol=nf,nrow=nvar) 
       for(f in 1:nf) {
     AUC[,f] <- colSums(temp[[f]])
     
   max.info[,f] <- apply(temp[[f]],2,which.max) }
	AUC <- AUC/lenx  #quasi normalize it 
	#max.info[is.nan(AUC) ] <- NA

	                                              
	max.info <- (max.info - lenx/2)*6/(lenx-1)
	max.info[max.info < -2.9] <- NA

	if(byKeys) {colnames(AUC) <- colnames(max.info) <- colnames(keys)} else {colnames(AUC) <- colnames(max.info) <- colnames(item$irt$discrimination)}
    rownames(AUC) <- rownames(max.info) <- rownames(item$rho)
    
    result <- list(AUC=AUC,max.info=max.info,sumInfo=sumtemp)
	invisible(result)				
     class(result) <- c("psych","polyinfo")
   invisible(result)
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


"irt.stats.like" <-function(items,stats=NULL,keys=NULL,cut=.3) {
     results <- list()
     tau <- irt.tau(items)
 if(!is.null(stats)) { nf <- dim(stats$loadings)[2]
          diffi <- list() 
          for (i in 1:nf) {diffi[[i]]  <- tau/sqrt(1-stats$loadings[,i]^2) }
         discrim <- stats$loadings/sqrt(1-stats$loadings^2)
        }
if(!is.null(keys)) {
        if(is.null(dim(keys))) { nf <- 1 } else {nf <- dim(keys)[2]}
        
          diffi <- list() 
          for (i in 1:nf) {diffi[[i]]  <- tau }
         discrim <- keys
        }
        
   class(diffi) <- NULL
   class(discrim) <- NULL
   irt <- list(difficulty=diffi,discrimination=discrim)
   results$irt <- irt
  
class(results) <- c("psych","irt.poly")
return(results)
}
   
   