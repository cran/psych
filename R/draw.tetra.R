"draw.tetra" <- 
function(r,t1,t2,shade=TRUE) {

 binBvn <- function(rho, rc, cc) {
        row.cuts <- c(-Inf, rc, Inf)
        col.cuts <- c(-Inf, cc, Inf)
        P <- matrix(0, 2, 2)
        R <- matrix(c(1, rho, rho, 1), 2, 2)
        for (i in 1:2) {
            for (j in 1:2) {
               # P[i, j] <- pmvnorm(lower = c(row.cuts[i], col.cuts[j]), upper = c(row.cuts[i + 1], col.cuts[j + 1]),  corr = R)
                  P[i, j] <- sadmvn(lower = c(row.cuts[i], col.cuts[j]), 
                upper = c(row.cuts[i + 1], col.cuts[j + 1]), mean=rep(0,2),
                varcov = R) 
            }
        }
        P
    }


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


HR <- 1 - pnorm(t1)
SR <- 1 - pnorm(t2)
P <- binBvn(r,HR,SR)
ph <- phi(P)
text(0,-.3,paste("phi = " ,round(ph,2)))

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
          
#figure out the equivalent phi 


par(def.par)  #reset to the original values
}



#show the density of a bivariate correlation
#adapted from the examples of persp
#meant to show how tetrachorics work
#8/4/14

draw.cor <- function(r=.5,expand=10,theta=30,phi=30,N=101,nbcol=30,box=TRUE,main="Bivariate density  rho = ",cuts=NULL,all=TRUE,ellipses=TRUE,ze=.15) {

sigma <- matrix(c(1,r,r,1),2,2)  #the covariance matrix
x <- seq(-3, 3, length= N)
y <- x
#f <- function(x, y,sigma=sigma) { r <- dmvnorm(cbind(x,y),sigma=sigma)}
f <- function(x, y,sigma=sigma) { r <- dmnorm(cbind(x,y),varcov=sigma)}
z <- outer(x,y,f,sigma=sigma)
nrz <- nrow(z)
ncz <- ncol(z)
jet.colors <- colorRampPalette( c("light blue", "red") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- jet.colors(nbcol)
color <- jet.colors(nbcol)
nbcol <- length(color)
#color <- terrain.colors(nbcol)
#color <- rainbow(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
if(is.null(cuts)) { #the default is to show a correlation
pmat <- persp(x, y, z, col = color[facetcol], phi = phi, theta = theta,scale=FALSE,expand=expand,box=box,main=paste(main,r))
} else {

if(!all) {z[,y > cuts[2] ]<- NA #just show the lower quarter
	z[ x > cuts[1], ]<- NA}
	z[,abs(y - cuts[2]) < 6/(N)] <- 0  #show all four quarters, with a notch in the distributions
	z[abs(x -cuts[1]) < 6/(N),] <- 0
 pmat <- 	persp(x, y, z, col = color[facetcol], phi = phi, theta = theta,scale=FALSE,expand=expand,box=box,main=paste(main,r))
	}
	
if(ellipses) {
 angles <- (0:N) * 2 * pi/N
        unit.circle <- cbind(cos(angles), sin(angles))
        if (abs(r) > 0) 
            theta1 <- sign(r)/sqrt(2)
        else theta1 = 1/sqrt(2)
        shape <- diag(c(sqrt(1 + r), sqrt(1 - r))) %*% matrix(c(theta1, 
            theta1, -theta1, theta1), ncol = 2, byrow = TRUE)
        ellipse <- unit.circle %*% shape*2
 lines (trans3d(ellipse[,1],ellipse[,2],z = ze, pmat = pmat),col = "red", lwd = 2)

}
#show it again, but without the ellipses
pmat <- 	persp(x, y, z, col = color[facetcol], phi = phi, theta = theta,scale=FALSE,expand=expand,box=box,main=paste(main,r))
}  
