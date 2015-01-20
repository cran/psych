"error.crosses.by" <-
function (x,y,z,labels=NULL,main=NULL,xlim=NULL,ylim= NULL,xlab=NULL,ylab=NULL,pos=NULL,offset=1,arrow.len=.2,alpha=.05,sd=FALSE,...)  # x  and y are data frame or descriptive stats
    {if(is.null(x$mean)) {x <- describe.by(x,z,mat=TRUE)
           }
     xmin <- min(x$mean)
     xmax <- max(x$mean)
     if(sd) {max.sex <- max(x$sd,na.rm=TRUE)
                      if(is.null(xlim))  {xlim=c(xmin - max.sex,xmax + max.sex) }}  else {max.sex <- max(x$se,na.rm=TRUE)}       
     if(is.null(y$mean)) {y <- describe(y)}
     ymin <- min(y$mean)
     ymax <- max(y$mean)
     if(sd) {max.sey <- max(y$sd,na.rm=TRUE)
              if(is.null(ylim))  {ylim=c(ymin - max.sey,ymax +max.sey)}} else {   max.sey <- max(y$se,na.rm=TRUE)  } 
     
     if(is.null(xlim))  xlim=c(xmin - 2*max.sex,xmax +2*max.sex)
     if(is.null(ylim))  ylim=c(ymin - 2*max.sey,ymax +2*max.sey)
     
     if(is.null(main)) {if(!sd) { main = paste((1-alpha)*100,"% confidence limits",sep="") } else {main= paste("Means and standard deviations")} }
     if(is.null(xlab)) xlab <- "Group 1"
     if(is.null(ylab)) ylab <- "Group 2"
     plot(x$mean,y$mean,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,...)
     
    cix <- qt(1-alpha/2,x$n-1) #modified Sept 11, 2013
    ciy <- qt(1-alpha/2,y$n-1)
     z <- dim(x)[1]
    if(sd) {x$se <- x$sd
            y$se <- y$sd
            cix <- ciy <- rep(1,z)
           }
    
     if (is.null(pos)) {locate <- rep(1,z)} else {locate <- pos}
     if (is.null(labels))  {labels <- rownames(x)}
    if (is.null(labels))  {lab <- paste("V",1:z,sep="")}  else {lab <-labels}
    
        for (i in 1:z)  
    	{xcen <- x$mean[i]
    	 ycen <- y$mean[i]
    	 xse  <- x$se[i]
    	 yse <-  y$se[i]
    	 arrows(xcen-cix[i]* xse,ycen,xcen+ cix[i]* xse,ycen,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 arrows(xcen,ycen-ciy[i]* yse,xcen,ycen+ ciy[i]*yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	text(xcen,ycen,labels=lab[i],pos=locate[i],cex=1,offset=offset)     #puts in labels for all points
    	}	
   }
   
   
   
   "ellipse" <-    function (x,y,r1,r2,...) { 
#code adapted from John Fox
    segments=51
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
   
    xs <- r1
    #ys <- e.size * yrange
    ellipse <- unit.circle 
    ellipse[,1] <- ellipse[,1]*r1 + x
    ellipse[,2] <- ellipse[,2]*r2+ y  #ys?
    lines(ellipse, ...)
    return(xs)
}
   
   "errorCircles" <-
function (x,y,data,ydata=NULL,group=NULL,paired=FALSE, labels=NULL,main=NULL,xlim=NULL,ylim= NULL,xlab=NULL,ylab=NULL,add=FALSE,pos=NULL,offset=1,arrow.len=.2,alpha=.05,sd=FALSE,bars=TRUE,circles=TRUE,...) { # x  and y are data frame or descriptive stats
     
     xvar <- x
     yvar <- y
    # if((length(xvar) ==1) && (length(yvar)==1) && !is.null(group)) {data <- statsBy(data[,c(group,xvar,yvar)],group=group)} else {
     if(!is.null(group)) {data <- statsBy(data,group=group)}
  #  }
    x <- list()
     if(paired) {
          	x$mean <- t(data$mean[,xvar])
     		x$sd <- t(data$sd[,xvar])
     		x$n <- t(data$n[,xvar]) 
     	} else {  #the normal case
    	 x$mean <- data$mean[,xvar]
     	x$sd <- data$sd[,xvar]
     	x$n <- data$n[,xvar]}  
     
     xmin <- min(x$mean,na.rm=TRUE)
     xmax <- max(x$mean,na.rm=TRUE)
     x$se <- x$sd/sqrt(x$n)
     
     if(sd) {max.sex <- max(x$sd,na.rm=TRUE)
                      if(is.null(xlim))  {xlim=c(xmin - max.sex,xmax + max.sex) }}  else {max.sex <- max(x$se,na.rm=TRUE)}       
     
     y <- list()
     if(!is.null(ydata)) {
          	y$mean <- ydata$mean[,yvar]
         	y$sd <- ydata$sd[,yvar]
        	y$n <- ydata$n[,yvar]
        	} else {
       	 	y$mean <- data$mean[,yvar]
    		 y$sd <- data$sd[,yvar]
    	 	y$n <- data$n[,yvar]}
     
     ymin <- min(y$mean,na.rm=TRUE)
     ymax <- max(y$mean,na.rm=TRUE)
     y$se <- y$sd/sqrt(y$n)
     if(sd) {max.sey <- max(y$sd,na.rm=TRUE)
              if(is.null(ylim))  {ylim=c(ymin - max.sey,ymax +max.sey)}} else {   max.sey <- max(y$se,na.rm=TRUE)  } 
     
     if(is.null(xlim))  xlim=c(xmin - 2*max.sex,xmax +2*max.sex)
     if(is.null(ylim))  ylim=c(ymin - 2*max.sey,ymax +2*max.sey)
     
     if(is.null(main)) {if(!sd) { main = paste((1-alpha)*100,"% confidence limits",sep="") } else {main= paste("Means and standard deviations")} }
     if(paired) {if(is.null(xlab)) xlab <- "Group 1"
            if(is.null(ylab)) ylab <- "Group 2"
        }  else {
        if(is.null(xlab)) xlab <- colnames(data$mean)[xvar]
        if(is.null(ylab)) ylab <- colnames(data$mean)[yvar]
     }
     if(add)  {  
      if(paired) {points(x$mean,typ="p",...) } else {points(x$mean,y$mean,typ="p",...)}
    } else {
     if(paired) {plot(x$mean,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,typ="p",...) } else {plot(x$mean,y$mean,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,typ="p",...)}
    }
     N <-x$n
     Nmax <- max(N)
    cix <- qt(1-alpha/2,x$n-1)
    ciy <- qt(1-alpha/2,y$n-1)
     if(paired) {z <- nrow(x$mean) } else {z <- length(x$mean)}
    if(sd) {x$se <- x$sd
            y$se <- y$sd
            cix <- ciy <- rep(1,z)
           }
    
     if (is.null(pos)) {locate <- rep(1,z)} else {locate <- pos}
     if (is.null(labels))  {labels <- rownames(x$mean)}
    if (is.null(labels))  {lab <- paste("V",1:z,sep="")}  else {lab <-labels}
    
        for (i in 1:z)  
    	{ if(paired) { xcen <- x$mean[i,1]
    	 ycen <- x$mean[i,2]
    	 xse  <- x$se[i,1]
    	 yse <-  x$se[i,2]
    	  } else {
    	xcen <- x$mean[i]
    	 ycen <- y$mean[i]
    	 xse  <- x$se[i]
    	 yse <-  y$se[i]}
    if(bars) {if(max(x$se,na.rm=TRUE) > 0) 	 arrows(xcen-cix[i]* xse,ycen,xcen+ cix[i]* xse,ycen,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	if(max(y$se,na.rm=TRUE) >0 )  arrows(xcen,ycen-ciy[i]* yse,xcen,ycen+ ciy[i]*yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 }
    	
    	text(xcen,ycen,labels=lab[i],pos=locate[i],cex=1,offset=offset)     #puts in labels for all points
    if(circles) { xrange <- xlim[2] - xlim[1]
                yrange <- ylim[2] - ylim[1]
                xscale <-max(x$se)
                yscale <-max(y$se)
    ellipse(xcen,ycen,sqrt(xscale*x$n[i]/Nmax),sqrt( yscale*x$n[i]/Nmax))
    	}
    	}	
    if(!is.null(group)) return(invisible(data))
   }
   
   
  

"cor.wt" <- function(data,vars=NULL, w=NULL,sds=NULL, cor=TRUE) {
   if(is.list(data) && !is.data.frame(data)) {w <- data$n   #use the output from statsBy
                   sds <- data$sd
                   x <- data$mean} else {x <- data}
    if(!is.null(vars)) {x <- x[,vars] 
                        w <- w[,vars]
                        sds <- sds[,vars] }              
   if(is.null(w)) w <- matrix(rep(rep(1/nrow(x),nrow(x)),ncol(x)),nrow=nrow(x),ncol=ncol(x))
   wt <- t(t(w)/colSums(w))
   cnames <- colnames(x)
    for (i in 1:ncol(x)) {if(is.factor(x[,i]) || is.logical(x[,i])) {
             x[,i] <- as.numeric(x[,i])
             colnames(x)[i] <- paste(cnames[i],"*",sep="")
             }}
   means <- colSums(x * wt,na.rm=TRUE) 
   xc  <-  scale(x,center=means,scale=FALSE)  #these are the weighted centered data
   if(is.null(sds)) {xs <- xc /sqrt(w) } else {xs  <- xc * sds/sqrt(w)}
   xwt <- sqrt(wt) * xc
       cov <- crossprod(xwt)             #/(1-colSums(wt^2,na.rm=TRUE))
    if(cor) {r <- cov2cor(cov)} else {r <- cov}
    xw <- wt * xc
return(list(r=r,xwt = xwt,wt=wt,mean=means,xc=xc,xs=xs))}
   
   
  
        
