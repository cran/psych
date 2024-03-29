#The diagram functions 

#The diagram function is the generic function (actually empty) that allows a clearer help file
#the entire set of functions for diagramming are all prefixed as dia. 
"diagram" <- function(fit,...) {
# figure out which psych function was called and then call the appropriate diagram
value <- NULL
omega <- bassAck <- extend <- extension <-  lavaan <- fa.reg <-  NULL    #weird requirement to avoid being global
 names <- cs(lm, lavaan,fa,principal,iclust,omega,bassAck, extension, extend, fa.reg, esem,mediate, lmCor, setCor, fa.multi,lm)
if(length(class(fit)) == 1) {
   value <- inherits(fit,names,which=TRUE) #if (inherits(fit,"lm")) value <- "lm" } 
   if(any(value > 0) ) { value <- names[which(value > 0)]} else {value <- "other"}
  }

if(length(class(fit)) > 1)  {
    names <- cs(lm, lavaan,fa,principal,iclust,omega,bassAck, extension, extend, fa.reg, esem,mediate, lmCor, setCor, fa.multi,lm)
    value <- inherits(fit,names,which=TRUE)  # value <- class(x)[2]
    if(any(value > 0) ) { value <- names[which(value > 0)]} else {value <- "other"}
    
     } 

switch(value,  #replaced a series of if statements to switch  12/23/18
fa  = {fa.diagram(fit,...) },
principal = {fa.diagram(fit,...) },
iclust = {iclust.diagram(fit,...)},
omega = {omega.diagram(fit,...)},
lavaan =  {lavaan.diagram(fit,...)}, 
bassAck = {bassAckward.diagram(fit,...)} ,
extend = {extension.diagram(fit ,...)},
extension ={extension.diagram(fit,...)},
fa.reg ={extension.diagram(fit,...)},
esem = {esem.diagram(fit,...)},
setCor = {lmDiagram(fit,...)},
lmCor = {lmDiagram(fit,...)},
lm = {lmDiagram(fit, ...)},
mediate ={mediate.diagram(fit,...)},
fa.multi ={fa.multi.diagram(fit,...)}
)
}

# 6/6/20 added the draw option to allow for speeding up higher level calls
# added the pmax to handle multiple calls 
#modified April 19, 2012 to handle long names more gracefully
"dia.rect" <- function(x, y = NULL, labels =NULL,  cex = 1, xlim=c(0,1),ylim=c(0,1),draw=TRUE,...) {
    if(draw)  text(x=x, y = y, labels = labels,  cex = cex,  ...)
      xrange = (xlim[2] - xlim[1])
    yrange = (ylim[2] - ylim[1])
    xs <- .10 * xrange
    ys <- .10 * yrange
     #len <- max(nchar(labes)*cex*.2/2,cex*.25)*xs
     len <- pmax(strwidth(labels,units="user",cex=cex,...),strwidth("abc",units="user",cex=cex,...)) /1.8   #was 1.8
     vert <- pmax(strheight(labels,units="user",cex=cex,...),strheight("ABC",units="user",cex=cex,...))/1.
    # vert <- min(cex*.3 * ys,.3)
    if(draw) rect(x-len,y-vert,x+len,y+vert)
     left <- c(x-len,y)
     right <- c(x+len,y)
     top <- c(x,y+vert)
     bottom <- c(x,y-vert)
     radius <- sqrt(len^2+vert^2)
    dia.rect <- list(left=left,right=right,top=top,bottom=bottom,center=c(x,y),radius=radius)
     }
     
"dia.rect1" <- function(x,y,xleft,ybot,xright,ytop, labels =NULL,  cex = 1, xlim=c(0,1),ylim=c(0,1),draw=TRUE,...) {
       text(x,y,labels,cex=cex,...)
       rect(xleft,ybot,xright,ytop)
}
     
 "dia.ellipse" <-  function(x, y = NULL, labels = NULL, cex = 1,e.size=.05,xlim=c(0,1),ylim=c(0,1),draw=TRUE, ...) {
    if(draw) text(x=x, y = y, labels = labels,cex = cex, ...)
     len <- max(strwidth(labels),strwidth("abc"))/1.6
     ht <-  max(strheight(labels),strheight("abc"))/1.0
     #vert <- cex*.3
     xs <- dia.ellipse1(x,y,xlim=xlim,ylim=ylim,e.size=e.size*len,draw=draw,...)
    # ys <- xs[2]
    # xs <- xs[1]
    ys <- xs* ht/len
    
     left <- c(x-xs,y)   #xs is a little too big
     right <- c(x+xs,y)
     top <- c(x,y+ys)
     bottom <- c(x,y-ys)
     center <- c(x,y)
    dia.ellipse <- list(left=left,right=right,top=top,bottom=bottom,center=center,radius=xs)
     }
          
     
"dia.ellipse1" <-    function (x,y,e.size=.05,xlim=c(0,1),ylim=c(0,1),draw=TRUE,...) { 
#code adapted from John Fox
    segments=51
    angles <- (0:segments) * 2 * pi/segments
    angles <- c(angles,NA)   #this allows us to draw multiple ellipses
    unit.circle <- cbind(cos(angles), sin(angles))
    xrange = (xlim[2] - xlim[1])
    yrange = (ylim[2] - ylim[1])
    xs <- e.size * xrange
    ys <- e.size * yrange
    ellipse <- unit.circle 
    ellipse[,1] <- ellipse[,1]*xs + x
    ellipse[,2] <- ellipse[,2]*ys + y  #ys?
   if(draw) lines(ellipse, ...)
   # return(list(xs=xs,ys=ys))
   return(xs)
}

"dia.ellipse.fill" <-    function (x,y,e.size=.05,xlim=c(0,1),ylim=c(0,1),draw=TRUE,density=NULL,angle=NULL,...) { 
#code adapted from John Fox
    segments=51
    angles <- (0:segments) * 2 * pi/segments
    angles <- c(angles,NA)   #this allows us to draw multiple ellipses
    unit.circle <- cbind(cos(angles), sin(angles))
    xrange = (xlim[2] - xlim[1])
    yrange = (ylim[2] - ylim[1])
    xs <- e.size * xrange
    ys <- e.size * yrange
    ellipse <- unit.circle 
    ellipse[,1] <- ellipse[,1]*xs + x
    ellipse[,2] <- ellipse[,2]*ys + y  #ys?
   if(draw) lines(ellipse, ...)
   polygon(ellipse,density=density,angle=angle,...)
   # return(list(xs=xs,ys=ys))
   return(xs)
}

#added to allow for fast drawing of multiple ellipses in fa.diagram and bassAckward diagram
dia.multi.ellipse <-  function(x, y = NULL, labels = NULL, cex = 1,e.size=.4,xlim=c(0,1),ylim=c(0,1),draw=TRUE, ...) {
 segments=51
    text(x=x, y = y, labels = labels,cex = cex, ...)
    angles <- (0:segments) * 2 * pi/segments
    angles <- c(angles,NA)   #this allows us to draw multiple ellipses
    unit.circle <- cbind(cos(angles), sin(angles))
    xrange = (xlim[2] - xlim[1])
    yrange = (ylim[2] - ylim[1])
  
  
    ys <- max(strheight(labels,units="user"))  #*.55
  
    # xs <- e.size * xrange
   #ys <- e.size * yrange
    ellipse <- unit.circle 

  
   xs  <- strwidth(labels,units="user")*.6 #this makes the ellipse closer to left and right

   #xs <- xs * e.size
   xs <- xs * cex
   ys <- ys * cex
     ellipsex <- rep(x,each=(segments + 2)) + unit.circle[,1] * rep(xs,each=segments +2)
     ellipsey <- rep(y,each=(segments + 2))  + unit.circle[,2] *ys
     
      lines(ellipsex,ellipsey)
}
       
"dia.triangle" <- function(x, y = NULL, labels =NULL,  cex = 1, xlim=c(0,1),ylim=c(0,1),...) {
     text(x=x, y = y, labels = labels,  cex = cex,   ...)
    STRETCH=.25
      xrange = (xlim[2] - xlim[1])
    yrange = (ylim[2] - ylim[1])
    xs <- .10 * xrange
    ys <- .10 * xrange
     #len <- max(nchar(labels)*cex*.2/2,cex*.25)*xs
     len <- max(strwidth(labels)*.7,strwidth("abc"))
     #vert <- min(cex*.3 * ys,.3)
     vert <- .7*len
     left <- c(x-len/2,y+vert/2-STRETCH/2)
     right <- c(x+len/2,y+vert/2-STRETCH/2)
     top <- c(x,y+vert)
    # bottom <- c(x,y-vert/2)
    bottom <- c(x,y-vert*STRETCH)
     xs <- c(x-len,x+len,x)
     ys <- c(y-vert*STRETCH,y-vert*STRETCH,y+vert)
     #ys <- c(y,y,y+vert)
     polygon(xs,ys)
     radius <- sqrt(len^2+vert^2)
    dia.rect <- list(left=left,right=right,top=top,bottom=bottom,center=c(x,y),radius=radius)
     }
     

     

       
       
"dia.curve" <- function(from,to,labels=NULL,scale=1,...) {
   #arrow at beginning and end of curve
   #first, figure out the boundaries of the curve
    radius1 <- radius2 <- 0
 	if(is.list(from)) {if(!is.null(from$radius)) {radius1 <- from$radius
 	        radius2 <- 0
    		from <- from$center}}
       if(is.list(to)) {if(!is.null(to$radius)) {radius2<- to$radius 
       			to <- to$center}}
       			
       	theta <- atan((from[2] - to[2])/(from[1] - to[1]))
        if((from[1] <  to[1])) {
                             radius1 <- -radius1
                           radius2 <- -radius2}
       	from[1] <-  from[1] - cos(theta) * radius1
       	from[2] <-  from[2] - sin(theta) * radius1
       	to[1] <- to[1] + cos(theta) * radius2
        to[2] <- to[2] + sin(theta) * radius2
        
      
    n <- 40
    scale <- .8 * scale
     if(is.null(labels)) {shift<- 0} else {shift <- 4}
    if(abs(from[1] - to[1]) < abs(from[2] - to[2]))  { #is it primarily up and down or left to right?
   		 x <- c(from[1],(from[1]+to[1])/2+scale,to[1])
   		 y <- c(from[2],(from[2]+to[2])/2,to[2])
   		sp <- spline(y,x,n)
 	 	lines(sp$y[1:(n/2-shift)],sp$x[1:(n/2-shift)])
     	lines(sp$y[(n/2+shift):n],sp$x[(n/2+shift):n])
     	arrows(sp$y[3],sp$x[3],sp$y[1],sp$x[1],length = .5*abs(sp$x[2]-sp$x[1]))
     	arrows(sp$y[n-3],sp$x[n-3],sp$y[n],sp$x[n],length = .5*abs(sp$x[2]-sp$x[1]))
  		text(sp$y[n/2],sp$x[n/2],labels,...)  
  		} else {
  		 x <- c(from[1],(from[1]+to[1])/2,to[1])
   		 y <- c(from[2],(from[2]+to[2])/2+ scale,to[2])
   		sp <- spline(x,y,n)
 	 	lines(sp$x[1:(n/2-shift)],sp$y[1:(n/2-shift)])
     	lines(sp$x[(n/2+shift):n],sp$y[(n/2+shift):n])
     	arrows(sp$x[3],sp$y[3],sp$x[1],sp$y[1],length = 1*abs(sp$y[2]-sp$y[1]))
     	arrows(sp$x[n-3],sp$y[n-3],sp$x[n],sp$y[n],length = 1*abs(sp$y[2]-sp$y[1]))
  		text(sp$x[n/2],sp$y[n/2],labels,...)  
  		}
}

"dia.curved.arrow" <- function(from,to,labels=NULL,scale=1,both=FALSE,dir=NULL,draw=TRUE,...) {
     #note that the splines seem to go from left to right no matter what!
     #first determine whether to add or subtract the radius
     result <- NULL
      radius1 <- radius2 <- 0
 	if(is.list(from)) {if(!is.null(from$radius)) {radius1 <- from$radius
 	        radius2 <- 0
    		from <- from$center}}
       if(is.list(to)) {if(!is.null(to$radius)) {radius2<- to$radius 
       			to <- to$center}}
       			
       	theta <- atan((from[2] - to[2])/(from[1] - to[1]))
       	if( from[1] < to[1] ) {
                             radius1 <- -radius1
                             radius2 <- -radius2}
       	from[1] <-  from[1] - cos(theta) * radius1
       	from[2] <-  from[2] - sin(theta) * radius1
       	to[1] <- to[1] + cos(theta) * radius2
        to[2] <- to[2] + sin(theta) * radius2
    
     scale <- .8 * scale
      n <- 40
     if(is.null(labels)) {shift<- 0} else {shift <- 4}
     #3 cases: dir="u", dir = "l", dir = NULL
      if(is.null(dir)) {
            if((abs(from[1] - to[1]) < abs(from[2] - to[2]))) {dir <- "u"} else {dir <- "l"}}
      switch(dir,
   # if(is.null(dir) & (abs(from[1] - to[1]) < abs(from[2] - to[2]))) 
      u =  { #is it primarily up and down or left to right?
   		 x <- c(from[1],(from[1]+to[1])/2+scale,to[1])
   		 y <- c(from[2],(from[2]+to[2])/2,to[2])
   		sp <- spline(y,x,n)
 	 if(draw) {	lines(sp$y[1:(n/2-shift)],sp$x[1:(n/2-shift)])
     	lines(sp$y[(n/2+shift):n],sp$x[(n/2+shift):n])
     	arrows(sp$y[3],sp$x[3],sp$y[1],sp$x[1],length = .5*abs(sp$x[2]-sp$x[1]))
     	if(both) arrows(sp$y[n-3],sp$x[n-3],sp$y[n],sp$x[n],length = .5*abs(sp$x[2]-sp$x[1]))
  		text(sp$y[n/2],sp$x[n/2],labels,...)  
  		} else {	result <- list(firstx=c(sp$y[1:(n/2-shift)],NA),
  	               firsty=c(sp$x[1:(n/2-shift)],NA),
  	               secondx = c(sp$y[(n/2+shift):n],NA),
  	               secondy=c(sp$x[(n/2+shift):n],NA),
  	               firstarrow = c(sp$y[3],sp$x[3],sp$y[1],sp$x[1],length = .5*abs(sp$x[2]-sp$x[1]) ),
  	               text = c(sp$y[n/2],sp$x[n/2],labels)
  	               )}
  		}, 
  	   l =  {
  		 x <- c(from[1],(from[1]+to[1])/2,to[1])
   		 y <- c(from[2],(from[2]+to[2])/2+ scale,to[2])
   		sp <- spline(x,y,n)
 	 	if(draw) {lines(sp$x[1:(n/2-shift)],sp$y[1:(n/2-shift)])
     	lines(sp$x[(n/2+shift):n],sp$y[(n/2+shift):n])
     	if(both) {
     	  arrows(sp$x[3],sp$y[3],sp$x[1],sp$y[1],length = 1*abs(sp$y[2]-sp$y[1])) 
     	  arrows(sp$x[n-3],sp$y[n-3],sp$x[n],sp$y[n],length = 1*abs(sp$y[2]-sp$y[1]))
     	   } else {if((from[1] > to[1] ) ) {arrows(sp$x[3],sp$y[3],sp$x[1],sp$y[1],length = 1*abs(sp$y[2]-sp$y[1])) } else {
     	arrows(sp$x[n-3],sp$y[n-3],sp$x[n],sp$y[n],length = 1*abs(sp$y[2]-sp$y[1]))}
  		text(sp$x[n/2],sp$y[n/2],labels,...)  }
  		} else {result <- list(firstx=c(sp$x[1:(n/2-shift)],NA),
  	               firsty=c(sp$y[1:(n/2-shift)],NA),
  	               secondx = c(sp$x[(n/2+shift):n],NA),
  	               secondy=c(sp$y[(n/2+shift):n],NA),
  	               firstarrow = c(sp$x[3],sp$y[3],sp$x[1],sp$y[1],length = .5*abs(sp$y[2]-sp$x[1]) ),
  	               text = c(sp$x[n/2],sp$y[n/2],labels)
  	               )
  		}
  		})
  		
  	
return(result) 		
}


 "dia.self" <- function(location,labels=NULL,scale=.8,side=2,draw=TRUE,...) {
   n <- 20
   
if(side %%2 > 0) {scale <- scale*(location$right[1] - location$left[1]) } else {scale <- scale*(location$top[2] - location$bottom[2]) }
 # scale <- scale *.8
 result <- NULL
if(side ==1)   
  {x <- c(location$bottom[1]-scale/2,location$bottom[1],location$bottom[1]+scale/2)
               y <- c(location$bottom[2],location$bottom[2]-scale,location$bottom[2])
               sp <- spline(x,y,n=20)
              if(draw) {lines(sp$x,sp$y)
               arrows(sp$x[3],sp$y[3],sp$x[1],sp$y[1],length = 2*abs(sp$x[3]-sp$x[1]))
     	       arrows(sp$x[n-3],sp$y[n-3],sp$x[n],sp$y[n],length = 2*abs(sp$x[3]-sp$x[1]))
     	        text(sp$x[n/2],sp$y[n/2]-scale,labels,...)}
     	        result <- list(sp$x[n/2],sp$y[n/2]-scale,labels, c(sp$x,NA),c(sp$y,NA),sp$x[3],sp$y[3],sp$x[1],sp$y[1],sp$x[n-3],sp$y[n-3],sp$x[n],sp$y[n], 2*abs(sp$x[3]-sp$x[1]))  }
if(side == 2)    {x <- c(location$left[1],location$left[1]-scale,location$left[1])
               y <- c(location$left[2]-scale/2,location$left[2],location$left[2]+scale/2)
               sp <- spline(y,x,n=20)
              if(draw) {lines(sp$y,sp$x)
               arrows(sp$y[3],sp$x[3],sp$y[1],sp$x[1],length = 2*abs(sp$x[2]-sp$x[1]))
     	      arrows(sp$y[n-3],sp$x[n-3],sp$y[n],sp$x[n],length = 2*abs(sp$x[2]-sp$x[1]))
     	      text(sp$y[n/2]-scale,sp$x[n/2],labels,...) }
     	      result <- list(sp$x[n/2],sp$y[n/2]-scale,labels, c(sp$x,NA),c(sp$y,NA),sp$y[3],sp$x[3],sp$y[1],sp$x[1],sp$y[n-3],sp$x[n-3],sp$y[n],sp$x[n],2*abs(sp$x[2]-sp$x[1]))  
     	       }
     	       
if(side == 3)   {x <- c(location$top[1]-scale/2,location$top[1],location$top[1]+scale/2)
               y <- c(location$top[2],location$top[2]+scale,location$top[2])
               sp <- spline(x,y,n=20)
                if(draw) {lines(sp$x,sp$y)
                 arrows(sp$x[3],sp$y[3],sp$x[1],sp$y[1],length = 2*abs(sp$x[3]-sp$x[1]))
     	         arrows(sp$x[n-3],sp$y[n-3],sp$x[n],sp$y[n],length = 2*abs(sp$x[3]-sp$x[1]))
     	           text(sp$x[n/2],sp$y[n/2]+scale,labels,...)}
     	            result <- list(sp$x[n/2],sp$y[n/2]-scale,labels, c(sp$x,NA),c(sp$y,NA),sp$x[3],sp$y[3],sp$x[1],sp$y[1],sp$x[n-3],sp$y[n-3],sp$x[n],sp$y[n],2*abs(sp$x[3]-sp$x[1]))  
     	           }
     	           
if(side == 4)    {x <- c(location$right[1],location$right[1]+scale,location$right[1])
               y <- c(location$right[2]-scale/2,location$right[2],location$right[2]+scale/2)
               sp <- spline(y,x,n=20)
              
             if(draw)  {lines(sp$y,sp$x)
               arrows(sp$y[3],sp$x[3],sp$y[1],sp$x[1],length =2* abs(sp$x[3]-sp$x[1]))
     	arrows(sp$y[n-3],sp$x[n-3],sp$y[n],sp$x[n],length = 2*abs(sp$x[3]-sp$x[1]))
     	  text(sp$y[n/2]+scale,sp$x[n/2],labels,...)}
     	  result <- list(sp$x[n/2],sp$y[n/2]-scale,labels, c(sp$x,NA),c(sp$y,NA),sp$y[3],sp$x[3],sp$y[1],sp$x[1],sp$y[n-3],sp$x[n-3],sp$y[n],sp$x[n],2*abs(sp$x[3]-sp$x[1]))  }
     return(result)	             
    }
    
  "dia.shape" <- function(x, y = NULL, labels = NULL, cex = 1,e.size=.05,xlim=c(0,1),ylim=c(0,1),shape=1, ...) {
     if (shape==1) {dia.rect(x, y = NULL, labels =NULL,    cex = 1, xlim=c(0,1),ylim=c(0,1),...) }
     if (shape==2) {dia.ellipse(x, y = NULL, labels = NULL, cex = 1,e.size=.05,xlim=c(0,1),ylim=c(0,1), ...)}
     if (shape==3) {dia.triangle(x, y = NULL, labels =NULL,  cex = 1, xlim=c(0,1),ylim=c(0,1),...)}
     }
     
 "multi.curved.arrow" <- function(curved.list,l.cex=NULL,...) {

leng  <- length(curved.list)
if(leng > 0) {
#nvar <- leng/6
seq <- seq.int(1,leng,6)
firstX <- unlist(curved.list[seq])
firstY <- unlist(curved.list[seq +1])
secondX <- unlist(curved.list[seq +2])
secondY <- unlist(curved.list[seq+3])
textvalues <- matrix(unlist(curved.list[seq+5]),ncol=3,byrow=TRUE)

lines(firstX,firstY)
lines(secondX,secondY)
if(is.null(l.cex) ) l.cex=1
text(x= textvalues[,1],y = textvalues[,2],labels=textvalues[,3],cex=l.cex,...)
 }}
 
 "multi.self" <- function(self.list,...) {
 leng <- length(self.list)
 seq <- seq.int(1,leng,14)
 x <- unlist(self.list[seq+3])
 y <- unlist(self.list[seq+4])
 lines(x,y)
 arrow.df <- data.frame(x0 = unlist(self.list[seq+5]),y0 = unlist(self.list[seq+6]),x1=unlist(self.list[seq+7]),y1 = unlist(self.list[seq+8]),
        x00 = unlist(self.list[seq+9]),y00 = unlist(self.list[seq+10]),x11=unlist(self.list[seq+11]),y11 = unlist(self.list[seq+12]))
 arrow.length <- unlist(self.list[seq+13])
 arrows(arrow.df$x0,arrow.df$y0,arrow.df$x1,arrow.df$y1,arrow.length[1])
 arrows(arrow.df$x00,arrow.df$y00,arrow.df$x11,arrow.df$y11,arrow.length[1])
 } 
  
     
"multi.arrow" <- function(arrows.list,...){
 tv <- matrix(unlist(arrows.list),byrow=TRUE,ncol=21)
   #cname<- colnames(tv)
   tv  <- data.frame(tv,stringsAsFactors=FALSE)
   tv[,c(1:2,4:19)] <- nchar2numeric(tv[,c(1:2,4:19)])
   if(nchar(as.vector(tv[1,21])) ==1) tv[,21] <- as.numeric(tv[,21 ]) #sometimes this is "solid", sometimes "1"  #the as.vector gets around a problem in R 3.6
   #colnames(tv) <- cname
  textlocation <- tv[,1:5]
  #colnames(textlocation) <- c("x","y","labels","pos","cex")
  textlocation1 <- tv[tv[,4]>0,1:5]
  textlocation0 <- tv[tv[,4]==0,1:5]
 if(NROW(textlocation1)>0 ) text(textlocation1[1:2],labels=textlocation1[,3],pos=textlocation1[,4],cex=textlocation1[,5])
  text(textlocation0[1:2],labels=textlocation0[,3],cex=textlocation0[,5])
      # text(tv[,1],tv[,2],tv[,3],pos=tv[,4],cex=tv[,5])
        len1 <- as.numeric(tv[1,10])
        len2 <- as.numeric(tv[1,17]) 
        arrows(x0=tv[,6],y0=tv[,7],x1=tv[,8],y1=tv[,9],length=len1,angle=30,code=1,col=tv[,20],lty=tv[,21],...)
        arrows(x0=tv[,13],y0=tv[,14],x1=tv[,15],y1=tv[,16],length=len2,angle=30,code=2,col=tv[,20],lty=tv[,21],...)

}
   
 
 
 multi.rect <-  function (rect.list, ...) 
{
    nvar <- length(rect.list)/7
    tv <- matrix(unlist(rect.list), nrow = nvar, byrow = TRUE)
    all.rects.x <- as.numeric(tv[, 6])
    all.rects.y <- as.numeric(tv[, 3])
    all.rects.names <- tv[, 1]
    dia.rect(all.rects.x, all.rects.y, all.rects.names)
}
  
# "multi.rect" <- function(rect.list,...) {
#   nvar <- length(rect.list)
#  tv <- matrix(unlist(rect.list),nrow=nvar,byrow=TRUE)
#       all.rects.x <- as.numeric(tv[,5]) #the center of the figure
#       all.rects.y <- as.numeric(tv[,2])
#       all.rects.names <-  tv[,1]    # i dont have names
#       dia.rect(all.rects.x, all.rects.y,all.rects.names) 
# }  
