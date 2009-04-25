
   
   "error.bars" <-
function (x,stats=NULL,ylab ="Dependent Variable",xlab="Independent Variable",main=NULL,ylim= NULL, alpha=.05, labels=NULL,pos=NULL,arrow.len=.05,add=FALSE,bars=FALSE,...)  # x   data frame with 
    {
    if(is.null(stats)) {
    	x.stats <- describe(x)
    	z <- dim(x)[2]
    	names <- colnames(x)
    	}  else { x.stats <- stats
    	          z <- dim(x.stats)[1]
    	          names <- rownames(stats)
    	}
    	min.x <- min(x.stats$mean)
    	max.x <- max(x.stats$mean)
    	max.se <- max(x.stats$se)
   		ci <- qt(1-alpha/2,x.stats$n)
    if(is.null(main)) main = paste((1-alpha)*100,"% confidence limits",sep="") 
    if(is.null(ylim)) {if(bars) {ylim=c(0,max.x+2*max.se) } else {ylim=c(min.x - 2*max.se,max.x+2*max.se)} }
    if(bars) {mp =barplot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,main=main,...)
     axis(1,mp[1:z],names)
     axis(2)
     box()
    } else {
    if(!add) {plot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,axes=FALSE,main=main,...)
     axis(1,1:z,names)
     axis(2)
     box()
     } else {points(x.stats$mean,...) }
     }  #end of if(bars)
    if(!is.null(labels)) {lab <- labels} else {lab <- paste("V",1:z,sep="")}
   
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     if (length(labels)==0) lab <- rep("",z) else lab <-labels
        for (i in 1:z)  
    	{xcen <- x.stats$mean[i]
    	# ycen <- y$mean[i]
    	 xse  <- x.stats$se[i]
    	# yse <-  y$se[i]
    	if(bars) {arrows(mp[i],xcen-ci*xse,mp[i],xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} else {
    	 arrows(i,xcen-ci*xse,i,xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 }
    	#text(xcen,i,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
    	}	
   }