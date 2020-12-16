"error.bars" <-
function (x,stats=NULL,data=NULL,group=NULL,ylab ="Dependent Variable",xlab="Independent Variable", main=NULL, eyes=TRUE, ylim=NULL, xlim=NULL, alpha=.05, sd=FALSE, labels=NULL,  
pos=NULL,arrow.len=.05,arrow.col="black",add=FALSE,bars=FALSE,within=FALSE,col=c("black","blue","red"),density=-10,...)  # x   data frame with 
    {    
     if(!missing(x) && (inherits(x, "formula"))) {#if(!is.null(data))   #if there is a grouping variable, call error.bars.by 
   # cat("\nFormula input detected, calling error.bars.by")
error.bars.by(x,data=data,x.cat=TRUE,ylab =ylab,xlab=xlab,main=main,ylim= ylim,
 eyes=eyes,alpha=alpha,sd=sd,labels=labels, v.labels=NULL, pos=pos, 
 arrow.len=arrow.len,add=add,bars=bars,within=within,colors=col, 
 legend=0,density=density,...)
     } else {
      if(!missing(group)) {  error.bars.by(x,group=group,x.cat=TRUE,ylab =ylab,xlab=xlab,main=main,ylim= ylim, eyes=eyes,alpha=.05,sd=sd,labels=labels, 
     v.labels=NULL, pos=pos, arrow.len=arrow.len, add=add,bars=bars, within=within, 
     colors=col, legend=0,density=density,...)}  else {
     
  #no grouping variable, this is the normal case
     
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
                                   ylim=c(min(0,min.x- 3*max.se),max.x+3*max.se)
                                   } else {
                                           ylim=c(min.x - 3*max.se,max.x+3*max.se)
                                           }}
                        }
    if(bars) {mp =barplot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,main=main,col=col,...)
     axis(1,mp[1:z],names)
     axis(2)
     box()
    } else {
    if(!add){
     if(missing(xlim)) {if (is.null(x.stats$values)) {xlim<- c(.5,z+.5) } else {xlim <- c(min(x.stats$values)-.5,max(x.stats$values)+.5)}}
       if(is.null(x.stats$values)) {
     plot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,xlim=xlim,axes=FALSE,main=main,...) 
       axis(1,1:z,names,...)
       axis(2)
     box()} else { plot(x.stats$values,x.stats$mean,ylim=ylim,xlim=xlim,xlab=xlab,ylab=ylab,main=main,...) }
     
     
     } else {points(x.stats$mean,...) }
     }  #end of if(bars)
     
    
    if(!is.null(labels)) {lab <- labels} else {lab <- paste("V",1:z,sep="")}
     
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     if (length(labels)==0) lab <- rep("",z) else lab <-labels
      s <- c(1:z)   #this allows us to address each one separately
    	        if(bars) {arrows(mp[s],x.stats$mean[s]-ci[s]* x.stats$se[s],mp[s],x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} else { 
    	       
   if(is.null(x.stats$values)) {
    	 arrows(s[s],x.stats$mean[s]-ci[s]* x.stats$se[s],s[s],x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3,col=arrow.col) } else {
    	 arrows(x.stats$values,x.stats$mean[s]-ci[s]* x.stats$se[s],x.stats$values,x.stats$mean[s]+ci[s]* x.stats$se[s],length=arrow.len, angle = 90, code=3, col=arrow.col)}
 	 
   	 
  if(eyes) { if(length(col) == 1) col <- rep(col,z)  
    	     ln <- seq(-3,3,.1)
    	     rev <- (length(ln):1)	   
    for (s in 1:z){ if(!is.null(x.stats$values[s] )) {
            catseyes(x=x.stats$values[s],y=x.stats$mean[s],se=x.stats$se[s],n=x.stats$n[s],alpha=alpha,density=density,col=col[s])} else {
            catseyes(x=s,y=x.stats$mean[s],se=x.stats$se[s],n=x.stats$n[s],alpha=alpha,density=density,col=col[s])}}
    	  }
    	    }
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
   #corrected April 22, 2016 to correctly handle the "stats" option with cats eyes
   
   
   
  "catseyes" <- function(x,y,se,n,alpha=alpha,density=density,col=col) {
     SCALE=.7
    	     ln <- seq(-3,3,.1)
    	     rev <- (length(ln):1)
    if(is.null(n) || is.na(n))  {norm <- dnorm(ln) 
     qt <- qnorm(alpha/2)
     clim <- qnorm(alpha/2)} else {if(n >1) {
    norm <-  dt(ln,n-1) 
    clim <- qt(alpha/2,n-1)}}
    norm <- c(norm,-norm[rev])
    ln <- seq(-3,3,.1)
    
    cln <- seq(clim,-clim,.01)
    cnorm <- dnorm(cln)
    cnorm <- c(0,cnorm,0,-cnorm,0)  #this closes the probability interval	  
    polygon(norm*SCALE+x,c(ln,ln[rev])*se+y)
    polygon(cnorm*SCALE+x,c(clim,cln,-clim,-cln,clim)*se+y,density=density,col=col)}
    
    
  #added May 30, 2016   
 error.bars.tab <- function(t,way="columns",raw=FALSE,col=c('blue','red'),...) {
rnames <- rownames(t)
cnames <- colnames(t)
t <- as.matrix(t)
switch(way,
columns = {p <- t(t(t)/colSums(t)) #convert to probabilities
       if(!raw) {standard.error  <- t(sqrt(t(p * (1-p))/colSums(t)))} else {standard.error  <- t(sqrt(t(p * (1-p))*colSums(t)))}},
rows  = {
      t <- as.matrix(t)
     p <- t/rowSums(t) #convert to probabilities
     if(!raw) {standard.error  <-sqrt(p * (1-p)/rowSums(t)) } else {standard.error  <-sqrt(p * (1-p)* rowSums(t))}},
both = {p <- t(t(t)/sum(t)) #convert to probabilities
       if(!raw) {standard.error  <- t(sqrt(t(p * (1-p))/sum(t)))} else{standard.error  <- t(sqrt(t(p * (1-p))*sum(t)))}  })
    

colnames(p) <-colnames(t)
rownames(p) <- rownames(t)
nv <- ncol(t)
ng <- nrow(t)
if(raw) {p <- t}
stats <- data.frame(mean=as.vector(p),
                        se=as.vector(standard.error))
rownames(stats) <- paste(rnames,rep(cnames,each=ng))
space <- rep(.1,nv*ng)
for(i in 1:(nv-1)) {space[ng*i + 1] <- 1}
error.bars(stats=stats,bars=TRUE,space=space,
 density=c(20,-10,20,-10),col=col,...)
 invisible(list(p=p,stats=stats))
}