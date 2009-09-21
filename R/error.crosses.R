"error.crosses" <-
function (x,y,labels=NULL,pos=NULL,arrow.len=.2,...)  # x  and y are data frame or descriptive stats
    {if(is.null(x$mean)) {x <- describe(x)
     xmin <- min(x$min)
     xmax <- max(x$max)}
     if(is.null(y$mean)) {y <- describe(y)
      ymin <- min(y$min)
     ymax <- max(y$max)
     plot(x$mean,y$mean,xlim=c(xmin,xmax),ylim=c(ymin,ymax),...)}
     
    
    z <- dim(x)[1]
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     if (length(labels)==0) lab <- rep("",z) else lab <-labels
     
        for (i in 1:z)  
    	{xcen <- x$mean[i]
    	 ycen <- y$mean[i]
    	 xse  <- x$se[i]
    	 yse <-  y$se[i]
    	 arrows(xcen-xse,ycen,xcen+xse,ycen,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 arrows(xcen,ycen-yse,xcen,ycen+yse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	text(xcen,ycen,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
    	}	
   }


