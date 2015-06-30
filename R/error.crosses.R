"error.crosses" <-
function (x,y,labels=NULL,main=NULL,xlim=NULL,ylim= NULL,xlab=NULL,ylab=NULL,pos=NULL,offset=1,arrow.len=.2,alpha=.05,sd=FALSE,add=FALSE,...)  # x  and y are data frame or descriptive stats
    {if(is.vector(x)) {x <- describe(x)}
     xmin <- min(x$mean)
     xmax <- max(x$mean)
     if(sd) {max.sex <- max(x$sd,na.rm=TRUE)
                      if(is.null(xlim))  {xlim=c(xmin - max.sex,xmax + max.sex) }}  else {max.sex <- max(x$se,na.rm=TRUE)}       
     if(is.vector(y)) {y <- describe(y)}
     ymin <- min(y$mean)
     ymax <- max(y$mean)
     if(sd) {max.sey <- max(y$sd,na.rm=TRUE)
              if(is.null(ylim))  {ylim=c(ymin - max.sey,ymax +max.sey)}} else {   max.sey <- max(y$se,na.rm=TRUE)  } 
     
     if(is.null(xlim))  xlim=c(xmin - 2*max.sex,xmax +2*max.sex)
     if(is.null(ylim))  ylim=c(ymin - 2*max.sey,ymax +2*max.sey)
     
     if(is.null(main)) {if(!sd) { main = paste((1-alpha)*100,"% confidence limits",sep="") } else {main= paste("Means and standard deviations")} }
     if(is.null(xlab)) xlab <- "Group 1"
     if(is.null(ylab)) ylab <- "Group 2"
     if(!add) plot(x$mean,y$mean,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,...)
     
    cix <- qt(1-alpha/2,x$n-1)
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
#Sept 11, 2013 changed n to n-1 in call to qt  (following a suggestion by Trevor Dodds)

