"error.bars.by" <-
function (x,group,data=NULL,by.var=FALSE,x.cat=TRUE,ylab =NULL,xlab=NULL,main=NULL,
	ylim= NULL, xlim=NULL, eyes=TRUE,alpha=.05,sd=FALSE,labels=NULL, 
	v.labels=NULL,v2.labels=NULL,add.labels=NULL,
 	pos=NULL,arrow.len=.05,min.size=1, add=FALSE,bars=FALSE, within=FALSE,
 	colors=c("black","blue","red"),lty,lines=TRUE, legend=0, pch=16,
 	density=-10,stats=NULL,...)  # x   data frame with 
    {

    if(!lines) {typ <- "p"} else {typ <- "b"}
    n.color <- length(colors)
n.pch =length(pch)
      #first, see if they are in formula mode   added August 18, 2018
  formula <- FALSE
   if(inherits(x, "formula")) {  ps <- fparse(x) 
     if(missing(data) | (length(data)==0)) {x <- get(ps$y)
	group <- ps$x } else {x <- data[ps$y]
	group <- data[ps$x]}
   
   formula <- TRUE
  #  if(is.null(data)) stop("You must specify the data if you are using formula input") 
#      x <- data[ps$y]
#    group <- data[ps$x]
    }
    
   grp2 <- n.grp2 <- 1   #just to set an initial value
   if(is.null(ylab)) ylab <- colnames(x)
   if(is.null(xlab)) xlab <- colnames(group)
   if(missing(by.var)) by.var=TRUE
   if(missing(lines)) lines <- FALSE
   if(NROW(group) <2 ) group <- x[,group]    
   if(NCOL(group)==1) {n.grp1 <- length(table(group))
     n.grp2 <- 1} else {n.grp1 <- (dim(table(group))[1])
     n.grp2 <- dim(table(group))[2] }
  
    nvar <- NCOL(x)
   # if(is.null(nvar)) nvar <- 1  #added May 21, 2016 to handle the case of a single variable
    if(by.var & (nvar > n.color)) {colors <- rainbow(nvar)}
    if(!missing(density)) {col12 <- col2rgb(colors,TRUE)/255
    colors <- rgb(col12[1,],col12[2,],col12[3,],.5)
     n.color <- length(colors)
     if(length(density==1)) density=rep(density,n.color)}
    #density = -10
    
    if(missing(lty)) lty <- 1:8
 
    legend.location <- c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right",  "center","none")
   if(is.null(stats)) {
   
    all.stats <- describe(x) } else {all.stats <- stats
    	          z <- dim(all.stats)[1]
    	          names <- rownames(stats)
    	          if(is.null(all.stats$se) ) stats$se <-  all.stats$se <- stats$sd/sqrt(stats$n-1)
    	          if(is.null(all.stats$max)) all.stats$max <- max(stats$mean,na.rm=TRUE)
    	           if(is.null(all.stats$min)) all.stats$min <- min(stats$mean,na.rm=TRUE)
    	}
        min.x <- min(all.stats$min,na.rm=TRUE)
   		max.x <- max(all.stats$max,na.rm=TRUE)
    	max.se <- max(all.stats$se,na.rm=TRUE)
    	if(sd) max.se <- max(all.stats$sd,na.rm=TRUE)
    if(is.null(ylim)) {if(is.na(max.x) | is.na(max.se) | is.na(min.x) | is.infinite(max.x)| is.infinite(min.x) | is.infinite(max.se)) {
                        ylim=c(0,1)} else {
     if(sd) { ylim <- c(min.x - max.se,max.x+max.se) } else {
    ylim=c(min.x - 2*max.se,max.x+2*max.se)}} }
    
    
    if(is.null(main)) {if(sd) {main <- paste("Means + Standard Deviations") } else {main <- paste((1-alpha)*100,"% confidence limits",sep="")} }

##########
# bar plot (ugly but some people like them) 
##########   
 if (bars) { #draw a bar plot and add error bars -- this is ugly but some people like it
   if(is.null(stats)){ group.stats <- describeBy(x,group,mat=TRUE,skew=FALSE)
      n.var <- dim(all.stats)[1]
       n.group <- length(group.stats[[1]])/n.var  #this is total number of groups but it may  be 2 x 2 or n x m
    } else {n.group <- NROW(stats)   #we have the stats as input, we just use them here
        group.stats <- list()
        group.stats$mean <- stats$mean
        group.stats$n <- stats$n
        group.stats$sd <- stats$sd
        group.stats$se <- stats$se
     } 
   # group.stats <- describeBy(x,group,mat=TRUE,skew=FALSE)   
          # n.var <- dim(all.stats)[1]
          # n.group <- length(group.stats[[1]])/n.var
          
          if(is.null(stats)){group.means <- matrix(group.stats$mean,ncol=n.group,byrow=TRUE)} else {group.means <- as.matrix(group.stats$mean,drop=FALSE)}
           
            if(missing(pch)) pch <- seq(15,(15+n.grp2))
          
           if(sd) {group.se <- matrix(group.stats$sd,ncol=n.group,byrow=TRUE)} else { group.se <- matrix(group.stats$se,ncol=n.group,byrow=TRUE)}
           group.n <- matrix(group.stats$n,ncol=n.group,byrow=TRUE)
           if(within) {group.smc <- matrix(unlist(by(x,group,smc)),nrow=n.group,byrow=TRUE)
                       group.sd <- matrix(group.stats$sd,ncol=n.group,byrow=TRUE)
                       if(sd) {group.se <- sqrt(group.se^2 * (1-group.smc))} else { group.se <- sqrt(group.sd^2 *(1-group.smc)/group.n) }}
           
           rownames(group.means) <- rownames(all.stats)
           if(is.null(labels)) {colnames(group.means) <- paste("Group",1:n.group)} else {colnames(group.means) <- labels }
           
           if (is.null(ylab)) ylab <- "Dependent Variable"
          
           if(missing(ylim)) ylim=c(0,max.x+2*max.se)
           if(by.var)  {
           if (is.null(xlab)) xlab <- "Variables"
           mp <- barplot(t(group.means),beside=TRUE,ylab=ylab,xlab=xlab,ylim=ylim,main=main,col=colors,...)
            for(i in 1:n.var) {
             for (j in 1:n.group) {
         	 xcen <- group.means[i,j]
    	 	 xse  <- group.se[i,j]
    	 	if(sd) {ci <- 1} else {ci <- qt(1-alpha/2,group.n[i,j])}
    	 	
       if(is.finite(xse) && xse>0)    arrows(mp[j,i],xcen-ci*xse,mp[j,i],xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
          }} } else {
           if (is.null(xlab)) xlab <- "Grouping Variable"
           mp <- barplot(group.means,beside=TRUE,ylab=ylab,xlab=xlab,ylim=ylim,main=main,col=colors,...)
            for(i in 1:n.var) {
             for (j in 1:n.group) {
         	 xcen <- group.means[i,j]
    	 	 
    	 	xse  <- group.se[i,j]
    	 	xn <- group.n[i,j]
    	 	if(sd) {ci <- 1} else {ci <- qt(1-alpha/2,group.n[i,j])}
        if(is.finite(xse) && xse>0 && xn > 0  && is.finite(ci))     arrows(mp[i,j],xcen-ci*xse,mp[i,j],xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = par("fg"), lty = NULL, lwd = par("lwd"), xpd = NULL)
          }} }
         
         
          axis(2,...)
          box()
            if(legend >0  ){
        
       if(!is.null(v.labels)) {lab <- v.labels} else {lab <- paste("V",1:n.var,sep="")} 
       legend(legend.location[legend], legend=lab, col = colors[(1: n.color)],pch=pch[1: n.var],
       text.col = "green4", lty = lty[1:n.var],
       merge = TRUE, bg = 'gray90')}
               
    } else {  


###############   
#the normal case is to not use bars
# could use some of the code from the bars section to be tighter code
###############
   if(is.null(stats)){ group.stats <- describeBy(x,group)
    n.group <- length(group.stats)   #this is total number of groups but it may  be 2 x 2 or n x m
    } else {n.group <- NROW(stats)   #we have the stats as input, we just use them here
        group.stats <- list()
      for(g in (1:n.group)) {group.stats[[g]] <- stats[g,]}
   
     }    #why not use matrix output for this
   
    n.var <- ncol(x)
     if(is.null(n.var)) n.var <- 1 
    #first set up some defaults to allow the specification of colors, lty, and pch dynamically and with defaults
    if(missing(pch)) pch <- seq(15,(15+n.group))
    if(missing(lty)) lty <- 1:8
     if(within) {group.smc <- by(x,group,smc) }
     z <- dim(x)[2]
     if(is.null(z)) z <- 1 
     
    if (is.null(ylab)) ylab <- "Dependent Variable"
     
     if(!by.var) {
     
        if(is.null(v.labels)) v.labels <- names(group.stats)   #use the values of the grouping variable
      	if (is.null(xlab)) xlab <- "Independent Variable"
    	for (g in 1:n.group) {
   	 	x.stats <- group.stats[[g]]   
   	 	if (within) { x.smc <- group.smc[[g]]  
    	              if(sd) {x.stats.$se <- sqrt(x.stats$sd^2* (1- x.smc))} else { x.stats$se <- sqrt((x.stats$sd^2* (1- x.smc))/x.stats$n)}
    	               }
        if (missing(xlim)) xlim <- c(.5,n.var+.5)
    	if(!add) {plot(x.stats$mean,ylim=ylim,xlim=xlim, xlab=xlab,ylab=ylab,main=main,typ=typ,lty=(lty[((g-1) %% 8 +1)]),axes=FALSE,col = colors[(g-1) %% n.color +1], pch=pch[g],...)
    	#axis(1,1:z,colnames(x),...)
    	axis(1,1:z,v.labels,...)
    	axis(2,...)
    	box()
    	} else {points(x.stats$mean,typ = typ,lty=lty[((g-1) %% 8 +1)],col = colors[(g-1) %% n.color +1],  pch=pch[g]) 
    	       }
    	if(!is.null(labels)) {lab <- labels} else {lab <- paste("V",1:z,sep="")}
   
   
   
     if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     	if (length(labels)==0) lab <- rep("",z) else lab <-labels
        for (i in 1:z)  {xcen <- x.stats$mean[i]
    	
    	 if(sd) {xse <- x.stats$sd[i] } else {xse  <- x.stats$se[i]}
    	
    	if(sd) {ci <- 1} else { if(x.stats$n[i] >1) {ci <- qt(1-alpha/2,x.stats$n[i]-1)} else {ci <- 0}}  #corrected Sept 11, 2013
    	  if(is.finite(xse) & (xse > 0))  {
    	  arrows(i,xcen-ci*xse,i,xcen+ci* xse,length=arrow.len, angle = 90, code=3,col = colors[(g-1) %% n.color +1], lty=(lty[((g-1) %% 8 +1)]), lwd = par("lwd"), xpd = NULL)
    	 if (eyes) {catseyes(i,xcen,xse,x.stats$n[i],alpha=alpha,density=density[(g-1) %% n.color +1],col=colors[(g-1) %% n.color +1] )}
    	#text(xcen,i,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
    	}
    	}
     add <- TRUE
   #  lty <- "dashed"
     }   #end of g loop 
     
       if(legend >0  ){
     
       
       if(!is.null(labels)) {lab <- labels} else {lab <- paste("G",1:n.group,sep="")} 
       legend(legend.location[legend], legend=lab, col = colors[(1: n.color)],pch=pch[1: n.group],
       text.col = "green4", lty = lty[1:8],
       merge = TRUE, bg = 'gray90')
   }
   #    #now, if v2.labels, get the var.means, etc 
#        n.vars <- dim(x)[2]
#     if(is.null(n.vars)) n.vars <- 1  #if we just have one variable to plot
#      var.means <- matrix(NA,nrow=n.vars,ncol=n.group)
#      var.n <- var.se <- ci <- matrix(0,nrow=n.vars,ncol=n.group)
#      for (g in 1:n.group) {
#    	 	var.means[,g] <- as.numeric (group.stats[[g]]$mean)    #problem with dimensionality -- if some grouping variables are empty
#    	 	if(sd) {var.se[,g] <- as.numeric(group.stats[[g]]$sd)}  else {var.se[,g] <- as.numeric(group.stats[[g]]$se) }
#    	 	var.n [,g] <-as.numeric( group.stats[[g]]$n)
#    	 	}

   } else  { # end of not by var loop
    
    #alternatively, do it by variables rather than by groups, or if we have two grouping variables, treat them as two variables
     if (is.null(xlab)) xlab <- "Grouping Variable"

    n.vars <- dim(x)[2]
    if(is.null(n.vars)) n.vars <- 1  #if we just have one variable to plot
     var.means <- matrix(NA,nrow=n.vars,ncol=n.group)
     var.n <- var.se <- ci <- matrix(0,nrow=n.vars,ncol=n.group)
     

 #if there are two or more grouping variables,this strings them out and treats them as two separate variables. 



    for (g in 1:n.group) {
   	 	var.means[,g] <- as.numeric (group.stats[[g]]$mean)    #problem with dimensionality -- if some grouping variables are empty
   	 	if(sd) {var.se[,g] <- as.numeric(group.stats[[g]]$sd)}  else {var.se[,g] <- as.numeric(group.stats[[g]]$se) }
   	 	var.n [,g] <-as.numeric( group.stats[[g]]$n)
   	 	}
     var.n[is.na(var.n)] <- 0
   	 var.n[var.n < min.size] <- 0
   	    if(x.cat) {x.values <- 1:n.grp1}  else {
   	      x.values <- names(table(group[,1]))
          x.values <- as.numeric(x.values)    }  
   	   if(missing(xlim)) xlim <- c(.5,max(x.values) + .5)
   	   if(is.null(v.labels)) v.labels <- names(unlist(dimnames(group.stats)[1]))

   	   
    for (i in 1:n.vars) {	


    	if(!add) {
    	
    	#although index by i, this just plots the first variable using plot, the other one is plotted using points
    	 plot(x.values,var.means[i,1:n.grp1],ylim=ylim,xlim = xlim, xlab=xlab,ylab=ylab,main=main,typ = typ,axes=FALSE,lty=lty[i],pch=pch[i],col = colors[(i-1) %% n.color +1],...)
    		  #we need to index by what grp1 we are in, not by i
    		  grp2 <- 1   #just in case we only have one IV
    	 if(n.grp1 < n.group) { #do we have multiple IVs
    	           # x.values1 <- rep(x.values,(n.group/n.grp1 -1))
          for(grp2 in 1:(n.grp2-1)){
    		    points(x.values,var.means[i,(n.grp1 *grp2 +1):(n.grp1*(grp2+1))],typ=typ,lty=lty[(grp2) %% 8 +1],col = colors[grp2 %% n.color + 1],pch=pch[(grp2 %% n.pch + 1)],...)
    		}
    	# points(x.values1,var.means[i,(n.grp1 +1):n.group],typ = typ,lty=lty[((i-1) %% 8 +1)],col = colors[(i) %% n.color + 1], pch=pch[i],...) #the first grouping variable
    	 }
    		add <- TRUE
    		} else {
    		  #now, draw the points with different line types and colors for each value of the second grouping variable
    		 #we want to adjust the line type for the fact that we have two different grouping variables 
    		
    		  points(x.values,var.means[i,1:(n.grp1)],typ = typ,lty=lty[((i-1) %% 8 +1)],col = colors[(i-1) %% n.color + 1], pch=pch[i],...) 
   if(n.grp2>1)  {		  
     for(grp2 in 1:(n.grp2-1)){
    		    points(x.values,var.means[i,(n.grp1 *grp2 +1):(n.grp1*(grp2+1))],typ=typ,lty=lty[(grp2-1) %% 8 +1],col = colors[(grp2-1) %% n.color + 1], pch=pch[((grp2-1) %% n.pch + 1)],...)
    		# if(n.grp1 < n.group) { points(x.values1,var.means[i,(n.grp1 +1):(n.group)],typ = typ,lty=lty[((1:(n.grp1)-1) %% 8 +1)],col = colors[(1:(n.grp1)) %% n.color + 1], pch=pch[i],...) }
    		   }     # points(x.values,var.means[i,],typ = typ,lty=lty,...)
    		}
    		}
    
    	if(!is.null(labels)) {lab <- labels} else {lab <- paste("G",1:z,sep="")}
   
    	if (length(pos)==0) {locate <- rep(1,z)} else {locate <- pos}
     	if (length(labels)==0) lab <- rep("",z) else lab <-labels

      #for (g in 1:n.group)  {
    		xcen <- var.means[i,]
    	 	xse  <- var.se[i,]
    	
    	   if(sd) {ci <- rep(1,n.group)} else { ci[(var.n > 1)] <- qt(1 - alpha/2, var.n[(var.n > 1)] - 1)
 }
    	  # }
    	  
   

  
    	   for (g in 1:n.grp1) {
    	   x.stats <- group.stats[[g]]  
    	   if(x.cat)  {   #the normal case
    	      # arrows(g,xcen[g]-ci[g]*xse[g],g,xcen[g]+ci[g]* xse[g],length=arrow.len, angle = 90, code=3, col = colors[(i-1) %% n.color +1], lty = NULL, lwd = par("lwd"), xpd = NULL)
    	 if(!is.na(var.n[g]) & (var.n [g] >1)) {if (eyes) { 
    	        catseyes(g,xcen[g],xse[g],group.stats[[g]]$n[i],alpha=alpha,density=density[(i-1) %% n.color +1],col=colors[(i-1) %% n.color +1] )   #why is this here?
    	                              }  else {
    	     
    	           if(ci[g]> 0) arrows(x.values[g],xcen[g]-ci[g]*xse[g],x.values[g],xcen[g]+ci[g]* xse[g],length=arrow.len, angle = 90, code=3,col = colors[(i-1) %% n.color +1], lty=(lty[((grp2-1) %% 8 +1)]), lwd = par("lwd"), xpd = NULL)
    	            if (eyes) {catseyes(x.values[g],xcen[g],xse[g],x.stats$n[g],alpha=alpha,density=density[(i-1) %% n.color +1],col=colors[(i-1) %% n.color +1] )}} 
    	            }
    	            }  else {#consider the !x.cat condition
    	                      }
    	        
    	  #text(xcen,i,labels=lab[i],pos=pos[i],cex=1,offset=arrow.len+1)     #puts in labels for all points
   		}
          x.values <- as.numeric(x.values)
   	if(n.grp1 < n.group) {
   		  for (g in 1: n.grp1) {
   		 
   		    for(grp2 in 1:(n.grp2-1)){if(!is.na(xse[g + n.grp1]) && (xse[g + n.grp1] > 0) && (ci[g + n.grp1]>0 )) {
   		    
   		   
   		    
   		    
    	              if(x.cat)  {
    	              arrows(g,xcen[g+n.grp1*grp2]-ci[g+n.grp1*grp2]*xse[g+n.grp1*grp2],g,xcen[g+n.grp1*grp2]+ci[g+n.grp1*grp2]* xse[g+n.grp1*grp2],length=arrow.len, angle = 90, code=3, col = colors[(grp2) %% n.color +1], lty=(lty[((grp2-1) %% 8 +1)]), lwd = par("lwd"), xpd = NULL)
    	                   } else {#consider the none cat condition 
    			                  x.values <- as.numeric(x.values)
    			                  arrows(x.values[g],xcen[g+n.grp1*grp2]-ci[g+n.grp1*grp2]*xse[g+n.grp1*grp2],x.values[g],xcen[g+n.grp1*grp2]+ci[g+n.grp1*grp2]* xse[g+n.grp1*grp2],length=arrow.len, angle = 90, code=3, col = colors[(grp2) %% n.color +1], lty=(lty[((grp2-1) %% 8 +1)]), lwd = par("lwd"), xpd = NULL)
    	                    }  
    	         	 if (eyes) { 
    	        catseyes(x.values[g],xcen[g+n.grp1*grp2],xse[g+n.grp1 * grp2],group.stats[[g+n.grp1*grp2]]$n[i],alpha=alpha,density=density[(i) %% n.color +1],col=colors[(i) %% n.color +1] )
    	                   }  else {
    	     
    	          # if(xse[g + n.grp1]> 0) 
    	          if(ci[g]>0) {arrows(x.values[g],xcen[g+n.grp1]-ci[g+n.grp1]*xse[g+n.grp1],x.values[g+n.grp1],xcen[g + n.grp1]+ci[g+n.grp1]* xse[g+n.grp1],length=arrow.len, angle = 90, code=3,col = colors[(grp2-1) %% n.color +1], lty=(lty[((grp2-1) %% 8 +1)]), lwd = par("lwd"),
    	        xpd = NULL)
    	            if (eyes) {catseyes(x.values[g],xcen[g+n.grp1],xse[g+n.grp1],x.stats$n[g+n.grp1],alpha=alpha,density=density[(i-1) %% n.color +1],col=colors[(i-1) %% n.color +1] )}} 
    	}
    	}
    	}  #for grp2 loop
   }
     }
     if(x.cat) {axis(1,1:n.grp1,v.labels,...) } else {axis(1)}
    		axis(2,...)
    		box() 
    		}   #end of i loop 
    }  #end of by var is true loop
    
   
    #### now add labels and lengeds if desired	
    
      #now, if v2.labels, get the var.means, etc 
       n.vars <- dim(x)[2]
    if(is.null(n.vars)) n.vars <- 1  #if we just have one variable to plot
     var.means <- matrix(NA,nrow=n.vars,ncol=n.group)
     var.n <- var.se <- ci <- matrix(0,nrow=n.vars,ncol=n.group)
     for (g in 1:n.group) {
   	 	var.means[,g] <- as.numeric (group.stats[[g]]$mean)    #problem with dimensionality -- if some grouping variables are empty
   	 	if(sd) {var.se[,g] <- as.numeric(group.stats[[g]]$sd)}  else {var.se[,g] <- as.numeric(group.stats[[g]]$se) }
   	 	var.n [,g] <-as.numeric( group.stats[[g]]$n)
   	 	}
    if(!by.var) var.means <- t(var.means)	
      if(!is.null(add.labels)) {
        
               if(is.null(v2.labels)) v2.labels <- names(unlist(dimnames(group.stats)[2]))
        if(is.null(pos)) {

        if(x.cat) {text.values <- data.frame(x=rep(.5 + (add.labels=="right")*NCOL(var.means),NROW(var.means)),y=var.means[,(add.labels=="right")*(NCOL(var.means)-1)+1],name=v2.labels)} else {
                   text.values <- data.frame(x=rep((.5+ (add.labels == "right") * x.values[n.grp1]),n.grp2),y=var.means[seq(1+(add.labels == "right")* (n.grp1),n.grp2 * n.grp1,(n.grp1))],name=v2.labels)}
        #if(x.cat){text.values <- data.frame(x=rep((.5+ (add.labels == "right")* n.grp1),n.grp2),y=var.means[seq(1+(add.labels == "right")* (n.grp1-1),n.grp2 * n.grp1,)],name=v2.labels)} else {
                  # text.values <- data.frame(x=rep((.5+ (add.labels == "right") * x.values[n.grp1]),n.grp2),y=var.means[seq(1+(add.labels == "right")* (n.grp1),n.grp2 * n.grp1,(n.grp1))],name=v2.labels)}
                   
                   text(text.values$x,text.values$y,text.values$name)
        } else {
            text.values <- data.frame(x=rep(.5 + (add.labels=="right")*NCOL(var.means)),y=var.means[,NCOL(var.means)],name=v2.labels ,position=pos)
           # text.values <- data.frame(x=rep((.5+ (add.labels == "right")* x.values[n.grp1]),n.grp2),y=var.means[,1+(add.labels == "right")* (n.grp1-1) ],name=v2.labels,position=pos)
        text(text.values$x,text.values$y,text.values$name,pos=text.values$position)
               }

          }
      if(legend > 0  ){
  
       if(!is.null(labels)) {lab <- labels} else {lab <- paste("V",1:n.grp2,sep="")} 
   if ((nvar > 1) & (n.grp2 ==1)) { legend(legend.location[legend], legend=lab, col = colors[(0: nvar)%% n.color + 1],pch=pch[1: n.grp2],
       text.col = "green4", lty = lty[1:nvar],
       merge = TRUE, bg = 'gray90')} else {
        legend(legend.location[legend],legend= lab, col = colors[(0: n.grp2) %% n.color + 1],pch=pch[1: n.grp2],
       text.col = "green4", lty = lty[1:n.grp2],
       merge = TRUE, bg = 'gray90') }   
                        }
     
    }  
   
    
  invisible(group.stats) }
   
   #corrected Feb 2, 2011 to plot alpha/2 rather than alpha 
   #modifed Feb 2, 2011 to not plot lines if they are not desired.
   #modified May 21, 2016 to handle a case of a single vector having no columns
   #modified April 9, 2019 to include v.labels for plots
   #modified July 11, 2020 to allow formula input
   #modified Novemember, 2020 to  handle NULL sample sizes
   #further modified November 2020 to properly label the lines and fix the pch to match the legend
   #And yet some more, December 31st to better handle non-=categorial variables
   #8/27/23  added the ability to plot v.labels for both by.var=TRUE and by.var=FALSE