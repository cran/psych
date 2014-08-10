"error.bars" <-
function (x,stats=NULL,ylab ="Dependent Variable",xlab="Independent Variable",main=NULL,eyes=TRUE,ylim= NULL,xlim=NULL, alpha=.05, sd=FALSE, labels=NULL,pos=NULL,arrow.len=.05,arrow.col="black",add=FALSE,bars=FALSE,within=FALSE,col="blue",...)  # x   data frame with 
    {
    SCALE=.5   #scale the width of the cats eyes
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
    	 {if(!sd) {
   		          if(is.null(stats)) {ci <- qt(1-alpha/2,x.stats$n-1) } else {ci <- rep(1,z) }}  else {ci <- sqrt(x.stats$n) 
   		           max.se <- max(ci * x.stats$se,na.rm=TRUE)} }
    if(is.null(main)) {if(!sd) { main = paste((1-alpha)*100,"% confidence limits",sep="") } else {main= paste("Means and standard deviations")} }
    if(is.null(ylim)) {if(is.na(max.x) | is.na(max.se) | is.na(min.x) | is.infinite(max.x)| is.infinite(min.x) | is.infinite(max.se)) {
                        ylim=c(0,1)} else {
                          if(bars) {
                                   ylim=c(min(0,min.x- 2*max.se),max.x+2*max.se)
                                   } else {
                                           ylim=c(min.x - 2*max.se,max.x+2*max.se)
                                           }}
                        }
    if(bars) {mp =barplot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,main=main,...)
     axis(1,mp[1:z],names)
     axis(2)
     box()
    } else {
    if(!add){
     if(missing(xlim)) xlim<- c(.5,z+.5)
       if(is.null(x.stats$values)) {
     plot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,xlim=xlim,axes=FALSE,main=main,...) 
       axis(1,1:z,names,...)
       axis(2)
     box()} else { plot(x.stats$values,x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,main=main,...) }
     
     
     } else {points(x.stats$mean,...) }
     }  #end of if(bars)
     
    
    if(!is.null(labels)) {lab <- labels} else {lab <- paste("V",1:z,sep="")}
     
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     if (length(labels)==0) lab <- rep("",z) else lab <-labels
      s <- c(1:z)
    	        if(bars) {arrows(mp[s],x.stats$mean[s]-ci[s]* x.stats$se[s],mp[s],x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} else { 
    	       
    	        if(is.null(x.stats$values)) {
    	 arrows(s[s],x.stats$mean[s]-ci[s]* x.stats$se[s],s[s],x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3,col=arrow.col) } else {
    	 arrows(x.stats$values,x.stats$mean[s]-ci[s]* x.stats$se[s],x.stats$values,x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3, col=arrow.col)}
 	 
   	 
  if(eyes) { if(length(col) == 1) col <- rep(col,z)  
    	     ln <- seq(-3,3,.1)
    	     rev <- (length(ln):1)	   
    for (s in 1:z){ if(!is.null(x.stats$n[s] )) {catseyes(x=s,y=x.stats$mean[s],se=x.stats$se[s],n=x.stats$n[s],alpha=alpha,density=-10,col=col[s])}}
    	  }
    	    }
   }
   #corrected July 25, 2009 to fix bug reported by Junqian Gordon Xu and then modified to be cleaner code
 
   #modified Sept 5, 2009 to handle data with all missing values (why we would want to that is a mystery, but was requested by Junqian Gordon Xu.)
   #April 5, 2010: the within parameter was added to allow for error bars in repeated measure designs 
   #modified June, 2010 to allow for easier use of stats input
   #modified June 15, 2010 to allow color bars to match color of lines and to have standard deviations as an option
   #modified Sept 11, 2013 to pass n -1 to the qt function (reported by Trevor Dodds)
   #modified March 10, 2014 add the capability to draw "cats eyes" 
   
   
   
  "catseyes" <- function(x,y,se,n,alpha=alpha,density=density,col=col) {
     SCALE=.7
    	     ln <- seq(-3,3,.1)
    	     rev <- (length(ln):1)
    if(n >1) {
    norm <-  dt(ln,n-1)
    norm <- c(norm,-norm[rev])
    ln <- seq(-3,3,.1)
    clim <- qt(alpha/2,n-1)
    cln <- seq(clim,-clim,.01)
    cnorm <- dnorm(cln)
    cnorm <- c(0,cnorm,0,-cnorm,0)  #this closes the probability interval	  
    polygon(norm*SCALE+x,c(ln,ln[rev])*se+y)
    polygon(cnorm*SCALE+x,c(clim,cln,-clim,-cln,clim)*se+y,density=density,col=col)}}