"error.bars.by" <-
function (x,group,by.var=FALSE,x.cat=TRUE,ylab =NULL,xlab=NULL,main=NULL,ylim= NULL, alpha=.05, labels=NULL,pos=NULL,arrow.len=.05,add=FALSE,bars=FALSE,...)  # x   data frame with 
    {
    lty = "solid"
    all.stats <- describe(x)
        min.x <- min(all.stats$min)
   		max.x <- max(all.stats$max)
    	max.se <- max(all.stats$se)
    if(is.null(ylim)) {ylim=c(min.x - 2*max.se,max.x+2*max.se)}
    if(is.null(main)) main <- paste(1-alpha,"% confidence limits",sep="")
    
    if (bars) { #draw a bar plot and add error bars -- this is ugly but some people like it
    group.stats <- describe.by(x,group,mat=TRUE)   
           n.var <- dim(all.stats)[1]
           n.groups <- length(group.stats[[1]])/n.var
           group.means <- matrix(group.stats$mean,ncol=n.groups,byrow=TRUE)
           group.se <- matrix(group.stats$se,ncol=n.groups,byrow=TRUE)
           group.n <- matrix(group.stats$n,ncol=n.groups,byrow=TRUE)
           rownames(group.means) <- rownames(all.stats)
           if(is.null(labels)) {colnames(group.means) <- paste("Group",1:n.groups)} else {colnames(group.means) <- labels }
           
           if (is.null(ylab)) ylab <- "Dependent Variable"
          
           ylim=c(0,max.x+2*max.se)
           if(by.var)  {
           if (is.null(xlab)) xlab <- "Variables"
           mp <- barplot(t(group.means),beside=TRUE,ylab=ylab,xlab=xlab,ylim=ylim,main=main)
            for(i in 1:n.var) {
             for (j in 1:n.groups) {
         	 xcen <- group.means[i,j]
    	 	 xse  <- group.se[i,j]
    	 	ci <- qt(1-alpha,group.n[i,j])
    	 	
       if(xse>0)    arrows(mp[j,i],xcen-ci*xse,mp[j,i],xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
          }} } else {
           if (is.null(xlab)) xlab <- "Grouping Variable"
           mp <- barplot(group.means,beside=TRUE,ylab=ylab,xlab=xlab,ylim=ylim,main=main)
            for(i in 1:n.var) {
             for (j in 1:n.groups) {
         	 xcen <- group.means[i,j]
    	 	 xse  <- group.se[i,j]
    	 	ci <- qt(1-alpha,group.n[i,j])
        if(xse>0)     arrows(mp[i,j],xcen-ci*xse,mp[i,j],xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
          }} }
         
         
          axis(2)
          box()
           
           
           
           } else {   #the normal case is to not use bars
           
           group.stats <- describe.by(x,group)
    n.group <- length(group.stats)
     z <- dim(x)[2]
     if(is.null(z)) z <- 1 
     
    if (is.null(ylab)) ylab <- "Dependent Variable"
     
     if(!by.var) {
      
      if (is.null(xlab)) xlab <- "Independent Variable"
    for (g in 1:n.group) {
   	 	x.stats <- group.stats[[g]]
   	 	

    	if(!add) {plot(x.stats$mean,ylim=ylim,xlab=xlab,ylab=ylab,main=main,typ="b",lty=((g-1) %% 8 +1),axes=FALSE,...)
    	axis(1,1:z,colnames(x))
    	axis(2)
    	box()
    	} else {points(x.stats$mean,typ="b",lty=((g-1) %% 8 +1),...) }
    
    	
    	if(!is.null(labels)) {lab <- labels} else {lab <- paste("V",1:z,sep="")}
   
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     	if (length(labels)==0) lab <- rep("",z) else lab <-labels
        for (i in 1:z)  
    	{xcen <- x.stats$mean[i]
    	# ycen <- y$mean[i]
    	 xse  <- x.stats$se[i]
    	# yse <-  y$se[i]
    	 ci <- qt(1-alpha,x.stats$n[i])
    	  if(xse>0)  arrows(i,xcen-ci*xse,i,xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 
    	#text(xcen,i,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
    	}
     add <- TRUE
     lty <- "dashed"
     }   #end of g loop 
    }    else  { # end of not by var loop
    
    #alternatively, do it by variables rather than by groups
     if (is.null(xlab)) xlab <- "Grouping Variable"
    
    n.vars <- dim(x)[2]
    if(is.null(n.vars)) n.vars <- 1  #if we just have one variable to plot
     var.means <- matrix(NA,nrow=n.vars,ncol=n.group)
     var.se <- matrix(NA,nrow=n.vars,ncol=n.group)
     

    for (g in 1:n.group) {
   	 	var.means[,g] <- group.stats[[g]]$mean
   	 	var.se[,g] <- group.stats[[g]]$se 
   	 	}
   	 	
   	 if(x.cat) {x.values <- 1:n.group}  else {
   	   x.values <- as.numeric(names(group.stats))  }  
   	 
    for (i in 1:n.vars) {	
   
    	if(!add) {
    	 plot(x.values,var.means[1,],ylim=ylim,xlab=xlab,ylab=ylab,main=main,typ="b",axes=FALSE,...)
    		if(x.cat) {axis(1,1:n.group,unlist(dimnames(group.stats))) } else {axis(1)}
    		axis(2)
    		box() 
    		add <- TRUE
    		} else {points(x.values,var.means[i,],typ="b",lty=((i-1) %% 8 +1),...) }
    	
    	if(!is.null(labels)) {lab <- labels} else {lab <- paste("G",1:z,sep="")}
   
    	if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     	if (length(labels)==0) lab <- rep("",z) else lab <-labels
       
        for (g in 1:n.group)  {
      
    		xcen <- var.means[i,g]
    	 	xse  <- var.se[i,g]
    	    ci <- qt(1-alpha,group.stats[[g]]$n)
    	   if(x.cat)  {arrows(g,xcen-ci*xse,g,xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)}  else {
    	    
    	            arrows(x.values[g],xcen-ci*xse,x.values[g],xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)} 
    	  #text(xcen,i,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
    		}
  
     lty <- "dashed"
     }   #end of i loop 
    }  #end of by var is true loop
    } # end of if not bars condition
   }