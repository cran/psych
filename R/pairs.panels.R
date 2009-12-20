"pairs.panels" <-
function (x, smooth = TRUE, scale = FALSE, density=TRUE,ellipses=TRUE,digits = 2, pch = 20,lm=FALSE,jiggle=FALSE,factor=2,...)   #combines a splom, histograms, and correlations
{

"panel.jiggle" <- 
function (x, y,  pch = par("pch"), 
   col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  
  xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise")
  x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)
  points(x, y, pch = pch, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)

  panel.ellipse1(xm,ym,xs,ys,r,col=col.smooth,...)
}
    op <- par(no.readonly = TRUE)  # save the whole list of settable par's.
    par(pch = pch)
    
    if(!lm) {if (density) {
      if (smooth) {
        if (scale) { 
            if(ellipses) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor.scale, lower.panel = panel.smoother, ...) } else {
             pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor.scale, lower.panel = panel.smooth, ...)}
        }  #don't scale
        else { #smooth and density 
        if(ellipses)  { if (jiggle) {  pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.jiggle , ...)
                       } else {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.smoother , ...) }} else {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.smooth, ...) }
          }
      # }
    }  #don't smooth
    else {
        if (scale) {
             if(ellipses) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor.scale,lower.panel=panel.ellipse,  ...) } else {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor.scale,  ...)}
        }
        else { #don't smooth, don't scale
            if(ellipses) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor,lower.panel=panel.ellipse, ...) } else {
              pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, ...) } 
        }
    }
    } else { # no densities
    if (smooth) {
        if (scale) {
             if(ellipses) {
            	pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale, lower.panel = panel.smoother, ...) } else {
            	pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale, lower.panel = panel.smooth, ...)}
        }  #don't scale
        else {
             if(ellipses) {
                       pairs(x, diag.panel = panel.hist, upper.panel = panel.cor, lower.panel = panel.smoother, ...) } else {
                       pairs(x, diag.panel = panel.hist, upper.panel = panel.cor, lower.panel = panel.smooth, ...)}
        }
    }  #don't smooth
    else {
        if (scale) {
         if(ellipses) {
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale,lower.panel=panel.ellipse,  ...) } else { 
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale,  ...)} 
        }
        else {
           if(ellipses) {
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor,lower.panel=panel.ellipse, ...) } else {
            pairs(x, diag.panel = panel.hist, upper.panel = panel.cor, ...)} 
        }
    }
    }
    } else { #lm is TRUE
    
        if (density) {
        	if(ellipses) {if (jiggle) {  pairs(x, diag.panel = panel.hist.density, upper.panel = panel.jiggle, lower.panel = panel.lm.ellipse , ...)
                       } else {
     			pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm.ellipse, lower.panel = panel.lm.ellipse, ...)}} else {
     			pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm, lower.panel = panel.lm, ...)   }
     	
     		} else {if(ellipses) {pairs(x, diag.panel = panel.hist, upper.panel = panel.lm, lower.panel = panel.lm.ellipse, ...)} else {
     	                      pairs(x, diag.panel = panel.hist, upper.panel = panel.lm, lower.panel = panel.lm, ...)}
   			 }
   			   }
    op <- par(op)
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
    points(x, y, pch = pch, ylim = ylim, xlim= xlim,...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(lm(y[ok]~ x[ok]), 
            col = col.lm, ...)
     xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise")
  panel.ellipse1(xm,ym,xs,ys,r,col=col.lm,...)
   
}

"panel.hist.density" <-
function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan")
  tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
  if(class(tryd) != "try-error") {
  
     d$y <- d$y/max(d$y)
   lines(d)}
}

"panel.hist" <-
function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan")
}


"panel.cor" <-
function(x, y, digits=2, prefix="", cex.cor,...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r  <- cor(x, y,use="pairwise")
         txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) {cex <- 0.8/strwidth(txt)} else {cex <- cex.cor}
         text(0.5, 0.5, txt,cex=cex)
     }
     
 "panel.cor.scale" <-
function(x, y, digits=2, prefix="", cex.cor,...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise"))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         cex1 <- cex  * abs(r)
         if(cex1 < .25) cex1 <- .25 #otherwise they just vanish
         text(0.5, 0.5, txt, cex = cex1)
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
   r = (cor(x, y,use="pairwise"))
  points(x, y, pch = pch, ...)
 angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
     if(!is.na(r)) {
  if (abs(r)>0 ) theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 

  shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
   ellipse <- unit.circle %*% shape 
   ellipse[,1] <- ellipse[,1]*xs + x
   ellipse[,2] <- ellipse[,2]*ys + y
   points(x,y,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
   lines(ellipse, ...)   }    
}
 


"panel.smoother" <- 
function (x, y,pch = par("pch"), 
    col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  
  xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise")
 
  points(x, y, pch = pch, ...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)

  panel.ellipse1(xm,ym,xs,ys,r,col=col.smooth,...)
}


