"pairs.panels" <-
function (x, smooth = TRUE, scale = FALSE, density=TRUE,ellipses=TRUE,digits = 2, pch = 20,lm=FALSE,jiggle=FALSE,  ...)   #combines a splom, histograms, and correlations
{
    op <- par(no.readonly = TRUE)  # save the whole list of settable par's.
    par(pch = pch)
    
    if(!lm) {if (density) {
      if (smooth) {
        if (scale) { 
            if(ellipses) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor.scale, lower.panel = panel.smoothie, ...) } else {
             pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor.scale, lower.panel = panel.smooth, ...)}
        }  #don't scale
        else {
        if(ellipses)  { #if (jiggle) {  } else {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.smoothie , ...) } else {
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
        else {
            if(ellipses) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor,lower.panel=panel.ellipse, ...) } else {
              pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, ...) } 
        }
    }
    } else { # no densities
    if (smooth) {
        if (scale) {
             if(ellipses) {
            	pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale, lower.panel = panel.smoothie, ...) } else {
            	pairs(x, diag.panel = panel.hist, upper.panel = panel.cor.scale, lower.panel = panel.smooth, ...)}
        }  #don't scale
        else {
             if(ellipses) {
                       pairs(x, diag.panel = panel.hist, upper.panel = panel.cor, lower.panel = panel.smoothie, ...) } else {
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
        	if(ellipses) {
     			pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm.ellipse, lower.panel = panel.lm.ellipse, ...)} else {
     			pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm, lower.panel = panel.lm, ...)   }
     	
     		} else {if(ellipses) {pairs(x, diag.panel = panel.hist, upper.panel = panel.lm, lower.panel = panel.lm.ellipse, ...)} else {
     	                      pairs(x, diag.panel = panel.hist, upper.panel = panel.lm, lower.panel = panel.lm, ...)}
   			 }
   			   }
    op <- par(op)
}

 "panel.lm" <- 
function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.lm = "red",  ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
    points(x, y, pch = pch, col = col, bg = bg, cex = cex,ylim = ylim, xlim= xlim)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(lm(y[ok]~ x[ok]), 
            col = col.lm, ...)
}

 "panel.lm.ellipse" <- 
function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.lm = "red",  ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
    points(x, y, pch = pch, col = col, bg = bg, cex = cex,ylim = ylim, xlim= xlim)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        abline(lm(y[ok]~ x[ok]), 
            col = col.lm, ...)
     xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise")
  p.ellipse(xm,ym,xs,ys,r,col=col.lm,...)
   
}

"panel.hist.density" <-
function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
   d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2)
    d$y <- d$y/max(d$y)
   lines(d)
}

"panel.hist" <-
function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}


"panel.cor" <-
function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise"))
         txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         text(0.5, 0.5, txt, cex = cex )
     }

"p.ellipse" <- 
function(x=0,y=0,xs=1,ys=1,r=0,add=TRUE,segments=51, col = par("col"),...) {
#based upon John Fox's ellipse functions
angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
  # shape <- diag(c(1+r,1-r)) %*% matrix(c(r,sqrt(1-r^2),-sqrt(1-r^2),r),ncol=2,byrow=TRUE)
  if (abs(r)>0 )theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 

 shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
   ellipse <- unit.circle %*% shape 
   ellipse[,1] <- ellipse[,1]*xs + x
   ellipse[,2] <- ellipse[,2]*ys + y
   
    points(x,y,pch=19,cex=1.5,col=col )
   lines(ellipse,col=col, ...)       
}


"panel.ellipse" <-
function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", ...) 
 {xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
   r = (cor(x, y,use="pairwise"))
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  p.ellipse(xm,ym,xs,ys,r,col=col.smooth,...)
}
 

"panel.smoothie" <- 
function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  
  xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise")
 
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)

  p.ellipse(xm,ym,xs,ys,r,col=col.smooth,...)
}

"panel.jiggle" <- 
function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  
  xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise")
  x <- jitter(x,factor=2)
  y <- jitter(y,factor=2)
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)

  p.ellipse(xm,ym,xs,ys,r,col=col.smooth,...)
}
