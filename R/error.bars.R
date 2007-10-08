"error.bars" <-
function (x,ylab ="Dependent Variable",xlab="", ylim= NULL, ci=1.96, labels=NULL,pos=NULL,arrow.len=.05,add=FALSE,...)  # x   data frame with 
    {
    x.stats <- describe(x)
    min.x <- min(x.stats$mean)
    max.x <- max(x.stats$mean)
    max.se <- max(x.stats$se)
    if(is.null(ylim)) {ylim=c(min.x - 2*max.se,max.x+2*max.se)}
    if(!add) {plot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,...) } else {points(x.stats$mean,...) }
     z <- dim(x)[2]
    if(!is.null(labels)) {lab <- labels} else {lab <- paste("V",1:z,sep="")}
   
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     if (length(labels)==0) lab <- rep("",z) else lab <-labels
        for (i in 1:z)  
    	{xcen <- x.stats$mean[i]
    	# ycen <- y$mean[i]
    	 xse  <- x.stats$se[i]
    	# yse <-  y$se[i]
    	
    	 arrows(i,xcen-ci*xse,i,xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 
    	#text(xcen,i,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
    	}	
   }