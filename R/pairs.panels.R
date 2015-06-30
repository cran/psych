#Adapted from the help for pairs
#modified December 15, 2011 to add the rug option
#further modified March 30, 2012 to add the method of correlation option (suggested by Carsten Dormann).
#by moving all the little functions to be outside the main function, this allows method and rug to be passed to these lower order functions.
#this should allow for somewhat cleaner code for the other functions
#modified March 15, 2015 to add the ability to control the size of the correlation separately from the cex variable in the points
#also added the ability to set the number of breaks in the histograms
#Also completely reorganized the main function to much cleaner

"pairs.panels" <-
function (x, smooth = TRUE, scale = FALSE, density=TRUE,ellipses=TRUE,digits = 2, method="pearson",pch = 20,lm=FALSE,cor=TRUE,jiggle=FALSE,factor=2,hist.col="cyan",show.points=TRUE,rug=TRUE, breaks="Sturges", cex.cor = 1 ,...)   #combines a splom, histograms, and correlations
{

#first define all the "little functions"


"panel.hist.density" <-
function(x,...) {
 usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x,breaks=breaks, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y,col=hist.col)
 if(density) { tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
  if(class(tryd) != "try-error") {
     d$y <- d$y/max(d$y)
   lines(d)}}
  if(rug) rug(x)
}


"panel.cor" <-
function(x, y, digits=2, prefix="",...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r  <- cor(x, y,use="pairwise",method=method)
         txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         cex <- cex.cor*0.8/strwidth(txt)
         if(scale)  {cex1 <- cex  * abs(r)
         if(cex1 < .25) cex1 <- .25 #otherwise they just vanish
         text(0.5, 0.5, txt, cex = cex1) } else {
         text(0.5, 0.5, txt,cex=cex)}
     }
        
"panel.smoother" <- 
function (x, y,pch = par("pch"), 
    col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  
  xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise",method=method)
 if(jiggle) { x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)}
 if(show.points)  points(x, y, pch = pch, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)

  panel.ellipse1(xm,ym,xs,ys,r,col.smooth=col.smooth,...)
}

"panel.smoother.no.noellipse" <- 
function (x, y,pch = par("pch"), 
    col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  
  xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise",method=method)
 if(jiggle) { x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)}
  if(show.points) points(x, y, pch = pch, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
}


 "panel.lm" <- 
function (x, y,  pch = par("pch"), 
    col.lm = "red",  ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
    if(jiggle) { x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)}
    points(x, y, pch = pch,ylim = ylim, xlim= xlim,...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(lm(y[ok]~ x[ok]), 
            col = col.lm, ...)
}

 "panel.lm.ellipse" <- 
function (x, y, pch = par("pch"), 
    col.lm = "red",  ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
    if(jiggle) { x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)}
    points(x, y, pch = pch, ylim = ylim, xlim= xlim,...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(lm(y[ok]~ x[ok]), 
            col = col.lm, ...)
     xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise",method=method)
  panel.ellipse1(xm,ym,xs,ys,r,col.smooth=col.lm,...)
   
}

 "panel.ellipse1" <- 
function(x=0,y=0,xs=1,ys=1,r=0,col.smooth,add=TRUE,segments=51,...) {
 #based upon John Fox's ellipse functions
 angles <- (0:segments) * 2 * pi/segments
     unit.circle <- cbind(cos(angles), sin(angles))
   # shape <- diag(c(1+r,1-r)) %*% matrix(c(r,sqrt(1-r^2),-sqrt(1-r^2),r),ncol=2,byrow=TRUE)
   if(!is.na(r)) {
   if (abs(r)>0 )theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 
 
 shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
    ellipse <- unit.circle %*% shape 
    ellipse[,1] <- ellipse[,1]*xs + x
    ellipse[,2] <- ellipse[,2]*ys + y
    points(x,y,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
    lines(ellipse, ...)   }    
 }

"panel.ellipse" <-
function (x, y,   pch = par("pch"), 
     col.smooth = "red", ...) 
 { segments=51
 xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
   r = cor(x, y,use="pairwise",method=method)
   if(jiggle) { x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)}
 if(show.points) points(x, y, pch = pch, ...)
 angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
     if(!is.na(r)) {
  if (abs(r)>0 ) theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 

  shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
   ellipse <- unit.circle %*% shape 
   ellipse[,1] <- ellipse[,1]*xs + xm
   ellipse[,2] <- ellipse[,2]*ys + ym
   points(xm,ym,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
  if(ellipses) lines(ellipse, ...) 
     }    
}
 





#Beginning of the main function
#The original organization was very clunky, but has now been cleaned up
#the organization gets around the problem of passing parameters to pairs by calling different routines
#this is truly clunky, but I don't know how else to do it
#It is possible that the trick I use for rug would work for other options.  -- YES 
#rug is a variable that is local to the entire function but is used only in the hist and hist.density functions.
#having done the same trick for passing method, it is now time to change all of this to pass some of these other options
#Done 15/15/03

old.par <- par(no.readonly = TRUE) # save default, for resetting... 
on.exit(par(old.par))     #and when we quit the function, restore to original values
#op <- par()
#old.par <- op$mfrow
#on.exit(par(old.par))
#old.par<- par(mfrow)    #get the original parameters so we can 
#on.exit(par(old.par))  #set them back upon exiting

if(missing(cex.cor)) cex.cor <- 1   #this allows us to scale the points separately from the correlations 

for(i in 1:ncol(x)) {
        if(is.character(x[[i]] ))  { x[[i]] <- as.numeric(as.factor(x[[i]]) )
                            colnames(x)[i] <- paste(colnames(x)[i],"*",sep="")}
           }
      
#par(pch = pch)

 
 #this has been greatly cleaned up by includign the jiggle, show.points and ellipses options  inside the lower.panel functions, and scale is inside the panel.cor function 
    
 if(!lm) { #the basic default is here
      if (smooth) {
        if(ellipses)  {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.smoother, pch=pch, ...)} else {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.smoother.no.noellipse,pch=pch,  ...)}    
        }
        else { #don't smooth, don't scale, not lm
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor,lower.panel=panel.ellipse,pch=pch,  ...) } 
    
    
    } else { #lm is TRUE
        if(!cor)  { #this case does not show the correlations, but rather shows the regression lines above and below the diagonal
                 if(ellipses) {
     			     pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm.ellipse, lower.panel = panel.lm.ellipse,pch=pch,  ...)} else {
     			     pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm, lower.panel = panel.lm, pch=pch, ...)   }
     			
     			} else {  #the normal case is to show the regressions below and the rs above
     			
        	     if(ellipses) {
     			    pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.lm.ellipse,pch=pch,  ...)} else {
     			    pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.lm,pch=pch,  ...)   } 
     	
     	
   			   }
   	}
#par(old.par) #reset the parameters
}




# "panel.hist" <-
# function(x, ...)
# {
#     usr <- par("usr"); on.exit(par(usr))
#    par(usr = c(usr[1:2], 0, 1.5) )
#     h <- hist(x, breaks=breaks,plot = FALSE)
#     breaks <- h$breaks; nB <- length(breaks)
#   y <- h$counts; y <- y/max(y)
#     rect(breaks[-nB], 0, breaks[-1], y,col=hist.col)
#      if(density) { tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
#   if(class(tryd) != "try-error") {
#      d$y <- d$y/max(d$y)
#      lines(d)}}
#  if(rug) rug(x)
# }








# "panel.ellipse1" <- 
# function(x=0,y=0,xs=1,ys=1,r=0,col.smooth,add=TRUE,segments=51,...) {
# #based upon John Fox's ellipse functions
# angles <- (0:segments) * 2 * pi/segments
#     unit.circle <- cbind(cos(angles), sin(angles))
#   # shape <- diag(c(1+r,1-r)) %*% matrix(c(r,sqrt(1-r^2),-sqrt(1-r^2),r),ncol=2,byrow=TRUE)
#   if(!is.na(r)) {
#   if (abs(r)>0 )theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 
# 
#  shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
#    ellipse <- unit.circle %*% shape 
#    ellipse[,1] <- ellipse[,1]*xs + x
#    ellipse[,2] <- ellipse[,2]*ys + y
#    points(x,y,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
#    lines(ellipse, ...)   }    
# }
# 
# "panel.ellipse.no" <-
# function (x, y,   pch = par("pch"), 
#      col.smooth = "red", ...) 
#  { segments=51
#  xm <- mean(x,na.rm=TRUE)
#   ym <- mean(y,na.rm=TRUE)
#   xs <- sd(x,na.rm=TRUE)
#   ys <- sd(y,na.rm=TRUE)
#    r = cor(x, y,use="pairwise",method=method)
#  if(jiggle) { x <- jitter(x,factor=factor)
#   y <- jitter(y,factor=factor)}
#  if (show.points ) points(x, y, pch = pch, ...)
#  angles <- (0:segments) * 2 * pi/segments
#     unit.circle <- cbind(cos(angles), sin(angles))
#      if(!is.na(r)) {
#   if (abs(r)>0 ) theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 
# 
#   shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
#    ellipse <- unit.circle %*% shape 
#    ellipse[,1] <- ellipse[,1]*xs + xm
#    ellipse[,2] <- ellipse[,2]*ys + ym
#    points(xm,ym,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
#    lines(ellipse, ...)   }    
# }



# "panel.smoother.no" <- 
# function (x, y,pch = par("pch"), 
#     col.smooth = "red", span = 2/3, iter = 3, ...) 
# {
#   
#   xm <- mean(x,na.rm=TRUE)
#   ym <- mean(y,na.rm=TRUE)
#   xs <- sd(x,na.rm=TRUE)
#   ys <- sd(y,na.rm=TRUE)
#   r = cor(x, y,use="pairwise",method=method)
#   if(jiggle) { x <- jitter(x,factor=factor)
#   y <- jitter(y,factor=factor)}
#   if(show.points) points(x, y, pch = pch, ...)
#     ok <- is.finite(x) & is.finite(y)
#     if (any(ok)) 
#         lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
#             col = col.smooth, ...)
# 
#   panel.ellipse1(xm,ym,xs,ys,r,col.smooth=col.smooth,...)
# }


# "panel.smoother.noellipse" <- 
# function (x, y,pch = par("pch"), 
#     col.smooth = "red", span = 2/3, iter = 3, ...) 
# {
#   
#   xm <- mean(x,na.rm=TRUE)
#   ym <- mean(y,na.rm=TRUE)
#   xs <- sd(x,na.rm=TRUE)
#   ys <- sd(y,na.rm=TRUE)
#   r = cor(x, y,use="pairwise",method=method)
#  if(jiggle) { x <- jitter(x,factor=factor)
#   y <- jitter(y,factor=factor)}
#  if(show.points)  points(x, y, pch = pch, ...)
#     ok <- is.finite(x) & is.finite(y)
#     if (any(ok)) 
#         lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
#             col = col.smooth, ...)
# }

#  "panel.cor.scale" <-
# function(x, y, digits=2, prefix="",...)
#      {
#          usr <- par("usr"); on.exit(par(usr))
#          par(usr = c(0, 1, 0, 1))
#          r = cor(x, y,use="pairwise",method=method)
#          txt <- format(c(r, 0.123456789), digits=digits)[1]
#          txt <- paste(prefix, txt, sep="")
#          cex <- 0.8*cex.cor/strwidth(txt)
#          cex1 <- cex  * abs(r)
#          if(cex1 < .25) cex1 <- .25 #otherwise they just vanish
#          text(0.5, 0.5, txt, cex = cex1)
#      }
#  


# "panel.jiggle" <- 
# function (x, y,  pch = par("pch"), 
#    col.smooth = "red", span = 2/3, iter = 3, ...) 
# {
#   
#   xm <- mean(x,na.rm=TRUE)
#   ym <- mean(y,na.rm=TRUE)
#   xs <- sd(x,na.rm=TRUE)
#   ys <- sd(y,na.rm=TRUE)
#   r = cor(x, y,use="pairwise",method=method)
#   x <- jitter(x,factor=factor)
#   y <- jitter(y,factor=factor)
#   points(x, y, pch = pch, ...)
#     ok <- is.finite(x) & is.finite(y)
#     if (any(ok)) 
#         lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
#             col = col.smooth, ...)
# 
#   panel.ellipse1(xm,ym,xs,ys,r,col.smooth=col.smooth,...)
# }
    