"dia.rect" <- function(x, y = NULL, labels =NULL,  cex = 1, xlim=c(0,1),ylim=c(0,1),...) {
     text(x=x, y = y, labels = labels,  cex = cex,   ...)
      xrange = (xlim[2] - xlim[1])
    yrange = (ylim[2] - ylim[1])
    xs <- .10 * xrange
    ys <- .10 * xrange
     len <- max(nchar(labels)*cex*.12/2,cex*.25)*xs
     vert <- min(cex*.3 * ys,.3)
     rect(x-len,y-vert,x+len,y+vert)
     left <- c(x-len,y)
     right <- c(x+len,y)
     top <- c(x,y+vert)
     bottom <- c(x,y-vert)
    dia.rect <- list(left=left,right=right,top=top,bottom=bottom)
     }
     
 "dia.ellipse" <-  function(x, y = NULL, labels = NULL, cex = 1,e.size=.05,xlim=c(0,1),ylim=c(0,1), ...) {
     text(x=x, y = y, labels = labels,cex = cex, ...)
     len <- max(nchar(labels)*cex*.14/2,cex*.2)
     vert <- cex*.3
     xs <- dia.ellipse1(x,y,xlim=xlim,ylim=ylim,e.size=e.size,...)
     left <- c(x-xs,y)
     right <- c(x+xs,y)
     top <- c(x,y+xs)
     bottom <- c(x,y-xs)
     center <- c(x,y)
    dia.ellipse <- list(left=left,right=right,top=top,bottom=bottom,center=center,radius=xs)
     }
     
     
     
"dia.ellipse1" <-    function (x,y,segments=51,e.size=.05,xlim=c(0,1),ylim=c(0,1),...) { 
angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    xrange = (xlim[2] - xlim[1])
    yrange = (ylim[2] - ylim[1])
    xs <- e.size * xrange
    ys <- e.size * xrange
    ellipse <- unit.circle 
    ellipse[,1] <- ellipse[,1]*xs + x
    ellipse[,2] <- ellipse[,2]*ys + y
    lines(ellipse, ...)
    return(xs)
}


"dia.arrow.1" <- function(from,to,labels=NULL,scale=1,radius=0,cex=1,...) {
       theta <- atan((from[2] - to[2])/(from[1] - to[1]))
        x <- (from[1] + to[1])/2
        y <- (from[2] + to[2])/2
        if(is.null(labels)) {h.size <- 0 } else{ h.size <- nchar(labels)*cex*.4}
       if(is.null(labels)) {v.size <- 0 } else{ v.size <- cex * .2}
        
        if(from[1] <  to[1] ) h.size <- -h.size
        x0 <- from[1] - cos(theta) * radius
        y0 <- from[2] - sin(theta) * radius
        xr <- x + h.size * cos(theta)* v.size
        xl <- x - h.size * cos(theta)* v.size
        yr <- y + v.size * sin(theta) * h.size 
        yl <- y - v.size * sin(theta)  *h.size
       if(!is.null(labels))  text(x,y,labels,cex=cex,...)
        arrows(x0,y0,xr,yr, length = 0, angle = 30, code = 2, ...)
        arrows(xl,yl,to[1],to[2], length = 0.1*scale, angle = 30, code = 2,...)
       }
       
       
"dia.arrow" <- function(from,to,labels=NULL,scale=1,radius1=0,radius2=0,cex=1,...) {
       theta <- atan((from[2] - to[2])/(from[1] - to[1]))
        x <- (from[1] + to[1])/2
        y <- (from[2] + to[2])/2
        if(is.null(labels)) {h.size <- 0 } else{ h.size <- nchar(labels)*cex*.4}
        if(is.null(labels)) {v.size <- 0 } else{ v.size <- cex * .2}
     
        if(from[1] <  to[1] ) h.size <- -h.size
        x0 <- from[1] - cos(theta) * radius1
        y0 <- from[2] - sin(theta) * radius1
        xr <- x + h.size * cos(theta)* v.size
        xl <- x - h.size * cos(theta)* v.size
        yr <- y + v.size * sin(theta) * h.size 
        yl <- y - v.size * sin(theta)  *h.size
        xe <- to[1] + cos(theta) * radius2
        ye <- to[2] + sin(theta) * radius2
       if(!is.null(labels))  text(x,y,labels,cex=cex,...)
        arrows(x0,y0,xr,yr, length = 0, angle = 30, code = 2, ...)
        arrows(xl,yl,xe,ye, length = 0.1*scale, angle = 30, code = 2,...)
       }
       

 
 
 
 
       
"dia.curve" <- function(from,to,labels=NULL,scale=1,...) {
    n <- 40
    scale <- .8 * scale
    if(is.null(labels)) {shift<- 0} else {shift <- 2}
    x <- c(from[1],(from[1]+to[1])/2+scale,to[1])
    y <- c(from[2],(from[2]+to[2])/2,to[2])
    sp <- spline(y,x,n)
 	 lines(sp$y[1:(n/2-shift)],sp$x[1:(n/2-shift)])
     lines(sp$y[(n/2+shift):n],sp$x[(n/2+shift):n])
  	text(sp$y[n/2],sp$x[n/2],labels,...)  
}

 "dia.self" <- function(location,labels=NULL,scale=.8,side=2,...) {
if(side %%2 < 1) {scale <- scale*(location$right[1] - location$left[1]) } else {scale <- scale*(location$top[2] - location$bottom[2]) }
 if (side==1) {x <- c(location$bottom[1]-scale/2,location$bottom[1],location$bottom[1]+scale/2)
               y <- c(location$bottom[2],location$bottom[2]-scale/2,location$bottom[2])
               sp <- spline(x,y,n=20)
               lines(sp$x,sp$y)}
 if (side==2) {x <- c(location$left[1],location$left[1]-scale/2,location$left[1])
               y <- c(location$left[2]-scale/2,location$left[2],location$left[2]+scale/2)
               sp <- spline(y,x,n=20)
               lines(sp$y,sp$x)}
 if (side==3) {x <- c(location$top[1]-scale/2,location$top[1],location$top[1]+scale/2)
               y <- c(location$top[2],location$top[2]+scale/2,location$top[2])
               sp <- spline(x,y,n=20)
                 lines(sp$x,sp$y)}
 if (side==4) {x <- c(location$right[1],location$right[1]+scale/2,location$right[1])
               y <- c(location$right[2]-scale/2,location$right[2],location$right[2]+scale/2)
               sp <- spline(y,x,n=20)
               lines(sp$y,sp$x)}
               
    }