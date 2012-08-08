"dia.cone" <- function(x=0, y=-2, theta=45, arrow=TRUE,curves=TRUE,add=FALSE,labels=NULL,xlim = c(-1, 1), ylim=c(-1,1),... ) {
 segments = 51
 extend = 1.1

xrange=2
#yrange=1
height= xrange
xs <-  tan(theta*pi/180)  * height 
ys =.3 * xs

 angles <- (0:segments) * 2 * pi/segments
 unit.circle <- cbind(cos(angles), sin(angles))
   

  
    ellipse <- unit.circle
    ellipse[, 1] <- ellipse[, 1] * xs + x
    ellipse[, 2] <- ellipse[, 2] * ys + y + height
    if(!add) {plot.new()
              plot.window(xlim=xlim*2,ylim=ylim,...)}
    lines(ellipse, ...)
  
   if(arrow) {
    arrows(x,y,(x-xs),y+ height,lty="dashed")
    arrows(x,y,(x + xs),y+ height,lty="dashed")
    arrows(x,y,x,y + extend^2 * height)} else {
 #don't draw arrows, just the cone 
    coords <- matrix(c(x,x-xs,y,y+height),2,2)
    lines(coords)
    coords <- matrix(c(x,x+xs,y,y+height),2,2)
    lines(coords)}
  if(curves) {dia.curve(c(x,y+height/3),c(x-xs/3,y+height/3),scale=.2,labels=labels[1])
              dia.curve(c(x,y+height/3),c(x+xs/3,y+height/3),scale=.2,labels=labels[2])
              dia.curve(c(x-xs/2,y+height/2),c(x+xs/2,y+height/2),scale=.3,labels=labels[3])
              }
}



