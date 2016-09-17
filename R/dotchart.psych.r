#completely taken from dotchart with some minor addition of confidence intervals


#July 17, 2016
#input is the mean + standard errors, and (optionally, alpha)

"error.dots" <- 
function (x,var=NULL, se=NULL, head = 12, tail = 12, sort=TRUE,decreasing=FALSE,alpha=.05,eyes=FALSE, min.n = NULL,max.labels =40, labels = NULL, groups = NULL, gdata = NULL, cex = par("cex"), 
    pt.cex = cex, pch = 21, gpch = 21, bg = par("bg"), color = par("fg"), 
    gcolor = par("fg"), lcolor = "gray", xlim = range(x[is.finite(x)]), 
    main = NULL, xlab = NULL, ylab = NULL, ...) 
{
    opar <- par("mai", "mar", "cex", "yaxs")
    on.exit(par(opar))
    par(cex = cex, yaxs = "i")
    #first, see if the data come from a psych object with sd and n or se 
    if(length(class(x)) > 1 ) {if (class(x)[1] == "psych") {obj <- class(x)[2]
    switch(obj,
        statsBy = {if(is.null(min.n)) {se <- x$sd[,var]/sqrt(x$n[,var]) 
                   x <- x$mean[,var] } else {se <- x$sd[,var]
                                             n.obs <- x$n[,var]
                                              x <- x$mean[,var]
                                              x <- subset(x,n.obs > min.n)
                                              se <- subset(se,n.obs > min.n)
                                              n.obs <- subset(n.obs, n.obs > min.n)
                                              se <- se/sqrt(n.obs)
                                              }},
        describe = {se <- x$se
        			labels <- rownames(x)
                    x <- x$mean
                    names(x) <- labels
                    } 
        )
      }  #end switch
      
     n.var <- length(x)
    if(!is.null(se)) {ci <- qnorm((1-alpha/2))*se} else {ci <- NULL}
    if (sort) { ord <- order(x,decreasing=decreasing) } else {ord <- n.var:1}   
    		 x <- x[ord]
   		 se <- se[ord]
   		 
   	 
   	
   	temp <- temp.se <- rep(NA,min(head+tail,n.var))
   	if((head+tail) < n.var) {
   	if (head > 0 ){ temp[1:head] <- x[1:head]
   	                temp.se[1:head] <- se[1:head]
   	                names(temp) <- names(x)[1:head]
   	                }
   	if(tail > 0 ) {temp[(head + 1):(head + tail)] <- x[(length(x)-tail+1):length(x)]
   	               temp.se[(head + 1):(head + tail)] <- se[(length(x)-tail+1):length(x)] 
   	                names(temp)[(head + 1):(head + tail)] <- names(x)[(length(x)-tail+1):length(x)]
   	               }
   	   
   	 x <- temp
   	 se <- temp.se
   	  }
   	 
   	 
   	  
   	 
   	 labels <- names(x)
   	 if(!is.null(se)) {ci <- qnorm((1-alpha/2))*se} else {ci <- NULL}
   	  if(!is.null(ci)) xlim <- c(min(x - ci),max(x + ci)) 
           
   	 labels <- substr(labels,1,max.labels)
   	 
   	   if(eyes) {  #get ready to draw catseyes
    	     ln <- seq(-3,3,.1)
    	     rev <- (length(ln):1)
    	     }
   	
   	        if (!is.numeric(x)) 
        stop("'x' must be a numeric vector or matrix")
    n <- length(x)
    if (is.matrix(x)) {
        if (is.null(labels)) 
            labels <- rownames(x)
        if (is.null(labels)) 
            labels <- as.character(1L:nrow(x))
        labels <- rep_len(labels, n)
        if (is.null(groups)) 
            groups <- col(x, as.factor = TRUE)
        glabels <- levels(groups)
    }
    else {
        if (is.null(labels)) 
            labels <- names(x)
        glabels <- if (!is.null(groups)) 
            levels(groups)
        if (!is.vector(x)) {
            warning("'x' is neither a vector nor a matrix: using as.numeric(x)")
            x <- as.numeric(x)
        }
    }
    plot.new()
    linch <- if (!is.null(labels)) 
        max(strwidth(labels, "inch"), na.rm = TRUE)
    else 0
    if (is.null(glabels)) {
        ginch <- 0
        goffset <- 0
    }
    else {
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- 0.4
    }
    if (!(is.null(labels) && is.null(glabels))) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + goffset, ginch) + 
            0.1
        par(mai = nmai)
    }
    if (is.null(groups)) {
        o <- 1L:n
        y <- o
        ylim <- c(0, n + 1)
    }
    else {
        o <- sort.list(as.numeric(groups), decreasing = TRUE)
        x <- x[o]
        groups <- groups[o]
        color <- rep_len(color, length(groups))[o]
        lcolor <- rep_len(lcolor, length(groups))[o]
        offset <- cumsum(c(0, diff(as.numeric(groups)) != 0))
        y <- 1L:n + 2 * offset
        ylim <- range(0, y + 2)
    }
    plot.window(xlim = xlim, ylim = ylim, log = "")
    lheight <- par("csi")
    if (!is.null(labels)) {
        linch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        loffset <- (linch + 0.1)/lheight
        labs <- labels[o]
        mtext(labs, side = 2, line = loffset, at = y, adj = 0, 
            col = color, las = 2, cex = cex, ...)
    }
    abline(h = y, lty = "dotted", col = lcolor)
    points(x, y, pch = pch, col = color, bg = bg, cex = pt.cex/cex)
    if(!is.null(ci)) {if(!eyes) {
    segments(x - ci, y, x+ci, y,
         col = par("fg"), lty = par("lty"), lwd = par("lwd"))
         } else {catseyes(x,y,se=se,n=NULL,alpha=alpha,density=-10)}}
    if (!is.null(groups)) {
        gpos <- rev(cumsum(rev(tapply(groups, groups, length)) + 
            2) - 1)
        ginch <- max(strwidth(glabels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(glabels, side = 2, line = goffset, at = gpos, adj = 0, 
            col = gcolor, las = 2, cex = cex, ...)
        if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
            points(gdata, gpos, pch = gpch, col = gcolor, bg = bg, 
                cex = pt.cex/cex, ...)
        }
    }
    axis(1)
    box()
    title(main = main, xlab = xlab, ylab = ylab, ...)
    invisible()
}
}


  
 #not used here  --   
  "catseyes" <- function(x,y,ci,density=density,col=col) {
     SCALE=.7
    	     ln <- seq(-3,3,.1)
    	     rev <- (length(ln):1)
     norm <- dnorm(ln) 
     clim <- qnorm(alpha/2)
   
    norm <-  dt(ln,n-1) 
    clim <- qt(alpha/2,n-1)
    norm <- c(norm,-norm[rev])
    ln <- seq(-3,3,.1)
   
    
    cln <- seq(clim,-clim,.01)
    cnorm <- dnorm(cln)
    cnorm <- c(0,cnorm,0,-cnorm,0)  #this closes the probability interval	  
    polygon(norm*SCALE+x,c(ln,ln[rev])*se+y)
    polygon(cnorm*SCALE+x,c(clim,cln,-clim,-cln,clim)*se+y,density=density,col=col)}
    
    
			