"draw.tetra" <- 
function(r,t1,t2,shade=TRUE) {
def.par <- par(no.readonly = TRUE) # save default, for resetting...
if(missing(r)) r <- .5
if(missing(t1)) t1 <- 1
if(missing(t2)) t2 <- 1
segments = 101

nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
#layout.show(nf)
par(mar=c(3,3,1,1))
#x amd y coordinates for an ellipse (code taken from John Fox)
 angles <- (0:segments) * 2 * pi/segments
        unit.circle <- cbind(cos(angles), sin(angles))
        if (abs(r) > 0) {
            theta <- sign(r)/sqrt(2)} else {theta = 1/sqrt(2)}
        shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta, 
            theta, -theta, theta), ncol = 2, byrow = TRUE)
        ellipse <- unit.circle %*% shape
      
       
x <- t1
y <- t2
xrange <-  c(-3,3)
yrange <- c(-3,3)
xloc = x+(3-x)/2
yloc <- y+ (3-y)/2
plot(x, y, xlim=xrange, ylim=yrange, xlab="X", ylab="Y",type="n")
lines(ellipse, type = "l")
ellipse3 <- ellipse*3
lines(ellipse3,type="l")

abline(v=x)
abline(h=y)
if(shade) {
   poly <- matrix(c(rep(t1,segments+1),rep(t2,segments+1)),ncol=2)  
   poly[,1] <- pmax(ellipse3[,1],t1)
   poly[,2] <- pmax(ellipse3[,2],t2)
polygon(poly ,density=10,angle=90)
         }
text(0,0,paste("rho = ",r))

text(xloc,yloc+.1,expression(X > tau))
text(xloc,yloc-.1,expression(Y > Tau))
text(-xloc,yloc+.1,expression( X < tau))
text(-xloc,yloc-.1,expression(Y > Tau))
text(xloc, -yloc + .1,expression( X > tau))
text(xloc,-yloc-.1,expression(Y < Tau))
text(-xloc,-yloc+.1,expression( X < tau))
text(-xloc,-yloc-.1,expression(Y < Tau))

#now the x distribution
par(mar=c(0,3,1,1))
curve(dnorm(x),-3,3,axes=FALSE)
lines(c(x,x),c(0,dnorm(x)))
text(xloc,dnorm(xloc)+.05,expression( X > tau))
text(t1,dnorm(t1) + .05,expression(tau))
if(shade) {
	xvals <- seq(t1,3,.1)
	yvals <- dnorm(xvals)
	polygon(c(xvals,rev(xvals)),c(rep(0,length(xvals)),dnorm(rev(xvals))),density=10,angle=-45)
}

#and the y distribution
par(mar=c(3,0,1,1))
x1 <- seq(-3,3,6/segments)
y1 <- dnorm(x1)
plot(y1,x1,axes=FALSE,typ="l")
lines(c(0,dnorm(y)),c(y,y))
text(.1,yloc,expression( Y > Tau))
text(dnorm(t2)+.05,t2,expression(Tau))
if(shade) {
	yvals <- seq(t2,3,.02)
	xvals <- dnorm(yvals)
	polygon(c(xvals,rev(xvals)),c(yvals,rep(t2,length(xvals))),density=10,angle=45)
          }
par(def.par)  #reset to the original values
}

