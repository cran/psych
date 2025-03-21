#Adapted from the help for pairs
#modified December 15, 2011 to add the rug option
#further modified March 30, 2012 to add the method of correlation option (suggested by Carsten Dormann).
#Fixed a bug in show.points on  March 18, 2017  (reported by Matthew Labrum)
#August 11, 2017  Added confidence intervals and adjust the histogram so it lines up with wht data points (suggested by Julan Martin) 
#by moving all the little functions to be outside the main function, this allows method and rug to be passed to these lower order functions.
#this should allow for somewhat cleaner code for the other functions
#modified March 15, 2015 to add the ability to control the size of the correlation separately from the cex variable in the points
#also added the ability to set the number of breaks in the histograms
#added the hist.border parameter as suggested by Jordan Adamson  (2/11/24)
#added ci.col and line.col parameters following  suggestions by Jordan Adamson (1/10/24)
"pairs.panels" <-
function (x, smooth = TRUE, scale = FALSE, density=TRUE,ellipses=TRUE,digits = 2, method="pearson",pch = 20,
    lm=FALSE,cor=TRUE,jiggle=FALSE,factor=2,hist.col="cyan",show.points=TRUE,rug=TRUE, breaks="Sturges", 
    cex.cor = 1 ,wt=NULL,smoother=FALSE,stars=FALSE,ci=FALSE,alpha=.05,
    hist.border="black", line.col="blue",ci.col="light blue",...)   #combines a splom, histograms, and correlations
{

#First define all the "little functions" that are internal to pairs.panels.  This allows easier coding later


"panel.hist.density" <-    
function(x,...) {
 usr <- par("usr"); on.exit(par("usr"))
  # par(usr = c(usr[1]-abs(.05*usr[1]) ,usr[2]+ abs(.05*usr[2])  , 0, 1.5) )  
     par(usr = c(usr[1] ,usr[2]  , 0, 1.5) )  
   tax <- table(x)
   if(length(tax) < 11) {breaks <- as.numeric(names(tax))
    y <- tax/max(tax)
    interbreak <- min(diff(breaks))*(length(tax)-1)/41
    rect(breaks-interbreak,0,breaks + interbreak,y,col=hist.col,border=hist.border)
    } else {
   
    h <- hist(x,breaks=breaks, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y,col=hist.col,border=hist.border)
    }
 if(density) { tryd <- try( d <- density(x,na.rm=TRUE,bw="nrd",adjust=1.2),silent=TRUE)
  if(!inherits(tryd,"try-error")) {
     d$y <- d$y/max(d$y)
   lines(d)}}
  if(rug) rug(x)
}


"panel.cor" <-
function(x, y, prefix="",...)  {
     
         usr <- par("usr"); on.exit(par("usr"))
         par(usr = c(0, 1, 0, 1))
        if(is.null(wt)) { r  <- cor(x, y,use="pairwise",method=method)} else {
        r <- cor.wt(data.frame(x,y),w=wt[,c(1:2)])$r[1,2]}
         txt <- format(c(round(r,digits), 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(stars) {pval <- r.test(sum(!is.na(x*y)),r)$p
                  symp <- symnum(pval, corr = FALSE,cutpoints = c(0,  .001,.01,.05, 1),
                symbols = c("***","**","*"," "),legend=FALSE)
                txt <- paste0(txt,symp)}
         cex <- cex.cor*0.8/(max(strwidth("0.12***"),strwidth(txt)))
         if(scale)  {cex1 <- cex  * abs(r)
         if(cex1 < .25) cex1 <- .25 #otherwise they just vanish
         text(0.5, 0.5, txt, cex = cex1) } else {
         text(0.5, 0.5, txt,cex=cex)}
     }
        
"panel.smoother" <- 
function (x, y,pch = par("pch"), 
    col.smooth = line.col, span = 2/3, iter = 3,ci.cor="light grey", ...) 
{
  # usr <- par("usr"); on.exit(par(usr))
 #  par(usr = c(usr[1]-abs(.05*usr[1]) ,usr[2]+ abs(.05*usr[2])  , usr[3],usr[4]) )     #doensn't affect the axis correctly
  xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
  r = cor(x, y,use="pairwise",method=method)
 if(jiggle) { x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)}
 if(smoother) {smoothScatter(x,y,add=TRUE, nrpoints=0)} else {if(show.points)  points(x, y, pch = pch, ...)}
 
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {   
     if(smooth & ci) {   lml <- loess(y~x ,degree=1,family="symmetric") 
     tempx <- data.frame(x = seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length.out=47))
       pred <-  predict(lml,newdata=tempx,se=TRUE ) 

 if(ci) {  upperci <- pred$fit + confid*pred$se.fit
       lowerci <- pred$fit - confid*pred$se.fit 
      polygon(c(tempx$x,rev(tempx$x)),c(lowerci,rev(upperci)),col=adjustcolor(ci.col, alpha.f=0.8), border=NA)
        }
     lines(tempx$x,pred$fit,  col = col.smooth, ...)   #this is the loess fit
}  else {if(smooth)  lines(stats::lowess(x[ok],y[ok],f=span,iter=iter),col=col.smooth) }}
 if(ellipses)  draw.ellipse(xm,ym,xs,ys,r,col.smooth=col.smooth,...)  #this just draws the ellipse 
}

 "panel.lm" <- 
function (x, y,  pch = par("pch"), 
    col.lm = "red", ci.col="blue", ...) 
{   ymin <- min(y)
    ymax <- max(y)
    xmin <- min(x)
    xmax <- max(x)
    ylim <- c(min(ymin,xmin),max(ymax,xmax))
    xlim <- ylim
    if(jiggle) { x <- jitter(x,factor=factor)
                 y <- jitter(y,factor=factor)}
  if(smoother) {smoothScatter(x,y,add=TRUE, nrpoints=0)} else {if(show.points) {points(x, y, pch = pch,ylim = ylim, xlim= xlim, ...)}}# if(show.points) points(x, y, pch = pch,ylim = ylim, xlim= xlim,...)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) {
      lml <- lm(y ~ x)  
    
       
	if(ci) {
   		 tempx <- data.frame(x = seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE),length.out=47))
   		 pred <-  predict.lm(lml,newdata=tempx,se.fit=TRUE)  #from Julian Martins 
   		 upperci <- pred$fit + confid*pred$se.fit
   		 lowerci <- pred$fit - confid*pred$se.fit
     	 polygon(c(tempx$x,rev(tempx$x)),c(lowerci,rev(upperci)),col=adjustcolor(ci.col, alpha.f=0.8), border=NA)
           }
    if(ellipses) {
    	xm <- mean(x,na.rm=TRUE)
 		ym <- mean(y,na.rm=TRUE)
  		xs <- sd(x,na.rm=TRUE)
  		ys <- sd(y,na.rm=TRUE)
   		r = cor(x, y,use="pairwise",method=method)
  		draw.ellipse(xm,ym,xs,ys,r,col.smooth=col.lm,...)   #just draw the ellipse
                }
         abline(lml, col = col.lm, ...)
  }
}


 "draw.ellipse" <-  function(x=0,y=0,xs=1,ys=1,r=0,col.smooth,add=TRUE,segments=51,...) {
                     #based upon John Fox's ellipse functions
 angles <- (0:segments) * 2 * pi/segments
     unit.circle <- cbind(cos(angles), sin(angles))
   if(!is.na(r)) {
   if (abs(r)>0 )theta <- sign(r)/sqrt(2) else theta=1/sqrt(2) 
 
 shape <- diag(c(sqrt(1+r),sqrt(1-r))) %*% matrix(c(theta,theta,-theta,theta),ncol=2,byrow=TRUE)
    ellipse <- unit.circle %*% shape 
    ellipse[,1] <- ellipse[,1]*xs + x
    ellipse[,2] <- ellipse[,2]*ys + y
    if(show.points) points(x,y,pch=19,col=col.smooth,cex=1.5 )  #draw the mean
    lines(ellipse, ...)   }    
 }

"panel.ellipse" <-
function (x, y,   pch = par("pch"), 
     col.smooth = "red", ...) 
 { segments=51
  usr <- par("usr"); on.exit(par("usr"))
  par(usr = c(usr[1]-abs(.05*usr[1]) ,usr[2]+ abs(.05*usr[2])  , 0, 1.5) ) 
 xm <- mean(x,na.rm=TRUE)
  ym <- mean(y,na.rm=TRUE)
  xs <- sd(x,na.rm=TRUE)
  ys <- sd(y,na.rm=TRUE)
   r = cor(x, y,use="pairwise",method=method)
   if(jiggle) { x <- jitter(x,factor=factor)
  y <- jitter(y,factor=factor)}
  if(smoother) {smoothScatter(x,y,add=TRUE, nrpoints=0)} else {if(show.points) {points(x, y, pch = pch, ...)}}
 	
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
#######

 #Beginning of the main function
######
#The original organization was very clunky, but has finally been cleaned up with lots of extra comments removed  (8/13/17)

old.par <- par(no.readonly = TRUE) # save default, for resetting... 
on.exit(par(old.par))     #and when we quit the function, restore to original values


if(missing(cex.cor)) cex.cor <- 1   #this allows us to scale the points separately from the correlations 

for(i in 1:ncol(x)) {  #treat character data as numeric
        if(is.character(x[[i]] ))  { x[[i]] <- as.numeric(as.factor(x[[i]]) )
                            colnames(x)[i] <- paste(colnames(x)[i],"*",sep="")}
           }
 n.obs <- nrow(x)     
 confid <- qt(1-alpha/2,n.obs-2)   #used in finding confidence intervals for regressions and loess
  
 if(!lm) { #the basic default is here
           if(cor) {
            pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor
            , lower.panel = panel.smoother, pch=pch, ...)} else {
        	pairs(x, diag.panel = panel.hist.density, upper.panel = panel.smoother, lower.panel = panel.smoother, pch=pch, ...)} 
    
    } else { #lm is TRUE
        if(!cor)  { #this case does not show the correlations, but rather shows the regression lines above and below the diagonal
                    pairs(x, diag.panel = panel.hist.density, upper.panel = panel.lm, lower.panel = panel.lm, pch=pch, ...)   
     			  } else {  #the normal case is to show the regressions below and the rs above
     			    pairs(x, diag.panel = panel.hist.density, upper.panel = panel.cor, lower.panel = panel.lm,pch=pch,  ...)   
     	
     	}
   			   }

   	}   #end of pairs.panels 
###

    
 "histo" <- function(x,breaks="Sturges", ...) {  
     tax <- table(x)
   if(length(tax) < 11) {breaks <- as.numeric(names(tax))
    y <- tax/max(tax)
    interbreak <- min(diff(breaks))*(length(tax)-1)/21
    rect(breaks-interbreak,0,breaks + interbreak,y)
    } else {
   
    h <- hist(x,breaks=breaks)
}}