"error.bars" <-
function (x,stats=NULL,ylab ="Dependent Variable",xlab="Independent Variable",main=NULL,ylim= NULL, alpha=.05, labels=NULL,pos=NULL,arrow.len=.05,add=FALSE,bars=FALSE,within=FALSE,...)  # x   data frame with 
    {
    if(is.null(stats)) {
    	x.stats <- describe(x)
    	if (within) { x.smc <- smc(x,covar=TRUE)
    	              x.stats$se <- sqrt((x.stats$sd^2 - x.smc)/x.stats$n)
    	               }
    	if(is.null(dim(x))) { z <- 1} else {z <- dim(x)[2]}  #if we just have one variable
    	names <- colnames(x)
    	}  else { x.stats <- stats
    	          z <- dim(x.stats)[1]
    	          names <- rownames(stats)
    	}
    	min.x <- min(x.stats$mean,na.rm=TRUE)
    	max.x <- max(x.stats$mean,na.rm=TRUE)
    	max.se <- max(x.stats$se,na.rm=TRUE)
   		ci <- qt(1-alpha/2,x.stats$n)
    if(is.null(main)) main = paste((1-alpha)*100,"% confidence limits",sep="") 
    if(is.null(ylim)) {if(is.na(max.x) | is.na(max.se) | is.na(min.x) | is.infinite(max.x)| is.infinite(min.x) | is.infinite(max.se)) {
                        ylim=c(0,1)} else {
                          if(bars) {
                                   ylim=c(0,max.x+2*max.se)
                                   } else {
                                           ylim=c(min.x - 2*max.se,max.x+2*max.se)
                                           }}
                        }
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
      s <- c(1:z)
    	        if(bars) {arrows(mp[s],x.stats$mean[s]-ci[s]* x.stats$se[s],mp[s],x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} else {
    	 arrows(s[s],x.stats$mean[s]-ci[s]* x.stats$se[s],s[s],x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL) }
    	
   }
   #corrected July 25, 2009 to fix bug reported by Junqian Gordon Xu and then modified to be cleaner code
   #modified Sept 5, 2009 to handle data with all missing values (why we would want to that is a mystery, but was requested by Junqian Gordon Xu.)
   #April 5, 2010: the within parameter was added to allow for error bars in repeated measure designs 