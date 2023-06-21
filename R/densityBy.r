"violin"  <- function(x,data=NULL,var=NULL,grp=NULL,grp.name=NULL, xlab=NULL, ylab=NULL,main="Density plot",vertical=TRUE,
   dots=FALSE, rain=FALSE,jitter=.05,
alpha=1,errors=FALSE,eyes=TRUE,adjust=1,restrict=TRUE,xlim=NULL,add=FALSE,col=NULL,pch=20,scale=NULL, ...) {

violinBy(x=x,data=data,var=var,grp=grp,grp.name=grp.name,xlab=xlab,ylab=ylab,main=main,vertical=vertical,
      dots=dots, rain=rain, jitter=jitter,alpha=alpha,errors=errors,eyes=eyes,adjust=adjust,restrict=restrict,xlim=xlim,add=add,col=col,pch=pch,scale=scale,...)
}

#switched from density = 50 to alpha =.5  to speed up the plotting 
"violinBy"  <- function(x,var=NULL,grp=NULL,data=NULL,grp.name=NULL,xlab=NULL, ylab=NULL, main="Density plot",vertical=TRUE,
    dots=FALSE,rain=FALSE, jitter=.05,alpha= 1,errors=FALSE,eyes=TRUE,adjust=1,restrict=TRUE,xlim=NULL,add=FALSE,col=NULL,pch=20,scale=NULL, ...) {
SCALE=.3  #how wide are the plots?

count.valid <- function(x) {sum(!is.na(x)) }

if(is.null(col)) {col <- c("blue","red","grey","purple","green","yellow")}

 formula <- FALSE
   if(inherits(x, "formula")) {  ps <- fparse(x)
   	formula <- TRUE
   if(is.null(data)) {data <- get(ps$y) 
      var <- NULL #we specified the entire dataframe} else { #stop("You must specify the data if you are using formula input") 
     x <- data} else {x <- data[,ps$y,drop=FALSE] }  #added the drop = FALSE 11/28/20
   grp <- data[ps$x]
   if(is.null(xlab)) xlab  <- ps$x 
   names <- ps$x
   grp.n <- ps$x
   if(is.null(ylab)) ylab <- ps$y
   }
  if(!formula) {ps <- list(x = grp, y = x)  #convert from x, grp   to y ~ x format 
     if(is.null(data))  {data <- x
       
       if(!is.null(var))  {data <- data[,var]}
        var <-  ps$y <- colnames(data)
        if(is.null(var))  var <- paste("V",1:NCOL(data))
    #    colnames(x) <- colnames(data) <- ps$y  <- var}
        colnames(data) <- ps$y  <- var}
        if(!is.null(grp)) {x <- cbind(data[,ps$y,drop=FALSE], x[,grp,drop=FALSE])} else {
     x <- data[,ps$y,drop=FALSE] }
     data <- x 
     if(is.null(ylab)) ylab = "Observed "
     if(is.null(xlab))  xlab = "" 
    }  
       if(!is.null(grp)) { 
    if(!is.data.frame(grp) && !is.list(grp) && (length(grp) < NROW(x))) grp <- x[,grp,drop=FALSE]}

 if(!is.null(var)) {if(missing(ylab) & (length(var) ==1)) {ylab <- var}
          x <- x[,var , drop = FALSE]}
 nvar <- nvarg <- NCOL(x)
 x <- char2numeric(x)
     
 col <- adjustcolor(col,alpha.f =alpha)        
#if(!is.null(grp)) {
# if(!is.data.frame(grp) && !is.list(grp) && (length(grp) < NROW(x))) grp <- x[,grp]

 if(!is.null(grp)) {stats <- describeBy(data[,ps$y],group = data[,ps$x],quant=c(0,1,.5,.25,.75),mat=TRUE,skew=FALSE) } else { stats <- describe(data[,ps$y],quant=c(0,1,.5,.25,.75),skew=FALSE)}
  meanX <- stats[,"mean"]
  minx <- stats[,"min"]
  maxx <- stats[,"max"]
  medx <- stats[,"Q0.5"]
   Q25 <- stats[,"Q0.25"]
  Q75  <- stats[,"Q0.75"]
  nvarg <- NROW(stats)
  n.col <- length(col)
  nX <- stats[,"n"]
  col <- rep(col,nvarg) #too many, but this will do

#  Qnt <-  apply(x,2,function(xx) by(xx,grp,quantile,prob=c(0,1,.5,.25,.75),na.rm=TRUE))
# meanX <- apply(x,2,function(xx) by(xx,grp,mean,na.rm=TRUE))
# nX <- apply(x,2,function(xx) by(xx,grp,count.valid))
# meanX <- matrix(unlist(meanX))
#    Qnt <- matrix(unlist(Qnt),nrow=5)
#    ngrp <- ncol(Qnt)/nvar
#    nvarg <- ncol(Qnt)
#    rangex <- matrix(c(Qnt[1,],Qnt[2,]),nrow=2,byrow=TRUE)
# } else {Qnt <- apply(x,2,quantile,prob=c(0,1,.5,.25,.75),na.rm=TRUE)
# meanX <- apply(x,2,mean,na.rm=TRUE)}
# minx <- Qnt[1,]
# maxx <- Qnt[2,]
# medx <- Qnt[3,]
# Q25 <-  Qnt[4,]
# Q75 <-  Qnt[5,]
# #rangex <- apply(x,2,range,na.rm=TRUE)
# rangex <- matrix(c(Qnt[1,],Qnt[2,]),nrow=2,byrow=TRUE)

tot.n.obs <- nrow(x)

                                                        
names <- grp.name  
if(is.null(grp.name) & !is.null(grp)) { grp.name <- paste(ps$x[1],stats[,"group1"])
names <- grp.name  } else {if(length(grp.name) != NROW(stats)) names <- rownames(stats)}

# names <- rep(grp.name, length(ps$x)) 
# #if(length(ps$x) > 1) ) ) {names <- paste(rep(names,each=ngrp),grp.name[1:ngrp],sep=" ")} else {names <- grp.name}
# if(length(ps$y) > 1) names <- rep(names,length(ps$x))
#      col <- rep(col,nvar* ngrp)
d <- list(nvar)
if(is.null(xlim)) xlim <- c(.5,nvarg+.5)

for (i in 1:nvar) {
#if(!is.null(grp)) { if(restrict) {d[[i]] <- by(x[,i], grp ,function(xx) density(xx,na.rm=TRUE,from=rangex[1,i],to=rangex[2,i]))}
if(!is.null(grp)) { if(restrict) {d[[i]] <- by(x[,i], grp ,function(xx) density(xx,na.rm=TRUE,adjust=adjust,from=min(xx,na.rm=TRUE),to=max(xx,na.rm=TRUE)))} else {
d[[i]] <- by(x[,i], grp ,function(xx) density(xx,na.rm=TRUE)) }} else {
if(restrict) {d[[i]] <- density(x[,i],na.rm=TRUE,adjust=adjust,from=minx[i],to=maxx[i])} else {
d[[i]] <- density(x[,i],na.rm=TRUE)} }
}

if(!add) {if(vertical){plot(1:nvarg,meanX,ylim=c(min(minx),max(maxx)),xlim=xlim,axes=FALSE,xlab=xlab,ylab=ylab,main=main,pch=pch,...)
  axis(1,1:nvarg,names,...)
  axis(2,...)} else {plot(meanX,1:nvarg,xlim=c(min(minx),max(maxx)),ylim=c(.5,nvarg+.5),axes=FALSE,xlab=ylab,ylab=xlab,main=main,pch=pch,...)
      axis(1,...)
     axis(2, 1:nvarg,names,...)}
  box()}
if(!is.null(grp)) d <- unlist(d,recursive=FALSE)
rev <- (length(d[[1]]$y):1) #this prevents a line down the middle


for(i in 1:nvarg) {
if(!is.null(scale)) {width <- scale*sqrt(nX[[i]]/tot.n.obs)/max(d[[i]]$y)} else {width <- SCALE/max(d[[i]]$y)}
#polygon(width*c(-d[[i]]$y,d[[i]]$y[rev])+i,c(d[[i]]$x,d[[i]]$x[rev]),density=density,col=col[i],...)
#if(vertical) {polygon(width*c(-d[[i]]$y,  d[[i]]$y[rev])+i,  c(d[[i]]$x,  d[[i]]$x[rev]),col=col[i ],...)} else {
 #        polygon(c(d[[i]]$x,  d[[i]]$x[rev]), width* c(-d[[i]]$y,  d[[i]]$y[rev])+i,col=col[i],...)}
zero <- rep(0,length(d[[i]]$x))

if(rain) {if(vertical) {polygon(width*c(d[[i]]$y,  zero)+i,c( d[[i]]$x,d[[i]]$x[rev]  ),col=col[i ],...)} else {
         polygon(c(d[[i]]$x,  d[[i]]$x[rev]), width* c(d[[i]]$y,  zero)+i,col=col[i],...)}} else {
         #not rain
   if(vertical) {polygon(width*c(-d[[i]]$y,  d[[i]]$y[rev])+i,  c(d[[i]]$x,  d[[i]]$x[rev]+i),col=col[i ],...)} else {
         polygon(c(d[[i]]$x,  d[[i]]$x[rev]), width* c(-d[[i]]$y,  d[[i]]$y[rev])+i,col=col[i],...)}}
         
#draw the quartiles and  medians        
dmd <- max(which(d[[i]]$x <= medx[i]))
d25 <- max(which(d[[i]]$x <= Q25[i]))
d75 <- max(which(d[[i]]$x <= Q75[i]))
if(vertical) {
segments(x0=width*d[[i]]$y[dmd] +i ,y0=d[[i]]$x[dmd],x1=-(rain-1)*width*d[[i]]$y[dmd]+i,y1=d[[i]]$x[dmd],lwd=2)
segments(x0=width*d[[i]]$y[d25] +i ,y0=d[[i]]$x[d25],x1=-(rain-1)*width*d[[i]]$y[d25]+i,y1=d[[i]]$x[d25])
segments(x0=width*d[[i]]$y[d75] +i ,y0=d[[i]]$x[d75],x1=-(rain-1)*width*d[[i]]$y[d75]+i,y1=d[[i]]$x[d75])
  } else {
  segments(y0=width*d[[i]]$y[dmd] +i ,x0=d[[i]]$x[dmd],y1=(rain-1)*width*d[[i]]$y[dmd]+i,x1=d[[i]]$x[dmd],lwd=2)
segments(y0=width*d[[i]]$y[d25] +i ,x0=d[[i]]$x[d25],y1=(rain-1)*width*d[[i]]$y[d25]+i,x1=d[[i]]$x[d25])
segments(y0=width*d[[i]]$y[d75] +i ,x0=d[[i]]$x[d75],y1=(rain-1)*width*d[[i]]$y[d75]+i,x1=d[[i]]$x[d75])
  }
 } 

 
 if(errors) {   #currently only works for one DV and one IV

 if(vertical) {
    error.bars.by(x[,1],group=grp[,1],add=TRUE,lines=FALSE,eyes=eyes,v.labels=names,...) 

 } else {
           message("error bars are not implemented for horizontal plots")}
 }
 
if(dots) { if(!is.null(grp)) {
if(jitter >0) {if(length(ps$x) > 1 ) {stripchart(x[,ps$y]~ grp[,ps$x[1] ] +grp[,ps$x[2]] ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch,...)} else {
        stripchart(x[,ps$y]~ grp[,ps$x[1] ] ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch,...) }} else {
  if(length(ps$x) > 1 ) {stripchart(x[,ps$y]~ grp[,ps$x[1]]+ grp[,ps$x[,2]],vertical=vertical,add=TRUE,pch=pch,...)} else {
       {stripchart(x[,ps$y]~ grp[,ps$x[1]],vertical=vertical,add=TRUE,pch=pch,...)} 
       }}
       } else {
       if(jitter >0) {stripchart(x ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch,...)} else {
        stripchart(x ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch,...) }} 
 
                   }
                   
if(rain)  { if(!is.null(grp)) {
 
 tempat <- 1:length(table(grp)) -.1
if(jitter >0) {if(length(ps$x) > 1 ) {stripchart(x[,ps$y]~ grp[,ps$x[1] ] +grp[,ps$x[2]] ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch, at=tempat,...)} else {
      
        stripchart(x[,ps$y]~ grp[,ps$x[1] ] ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch,at=tempat, ...) }} else {
  if(length(ps$x) > 1 ) {stripchart(x[,ps$y]~ grp[,ps$x[1]]+ grp[,ps$x[,2]],vertical=vertical,add=TRUE,pch=pch,...)} else {
       {stripchart(x[,ps$y]~ grp[,ps$x[1]],vertical=vertical,add=TRUE,pch=pch,...)} 
       }}
       } else {
       if(jitter >0) {stripchart(x ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch,...)} else {
        stripchart(x ,vertical=vertical,method="jitter",jitter=jitter,add=TRUE,pch=pch,...) }} 
 
                   }
       
       }
        
#created March 10, 2014 following a discussion of the advantage of showing distributional values
#modified March 22, 2014 to add the grouping variable
#basically just a violin plot.
#modified December, 2016 to allow for scaling of the widths of the plots by sample size. 
#modified  December 2020 to provide stripcharts and errorbar options

# 
# "histBy" <- function(x,grp=NULL,data=NULL,restrict=TRUE,xlim=NULL,ylab="Observed",xlab="",main="Density plot",density=20,scale=TRUE,col= c("blue","red"),...) {
# count.valid <- function(x) {sum(!is.na(x)) }
# formula <- FALSE
#    if(inherits(x, "formula")) {  ps <- fparse(x)
#    formula <- TRUE
#    if(is.null(data)) stop("You must specify the data if you are using formula input") 
#      x <- data[ps$y]
#    grp <- data[ps$x]
#    }
# if(
#(col)) {col <- c("blue","red")}
# if(restrict) {minx <- min(x,na.rm=TRUE)
#    maxx <- max(x,na.rm=TRUE)}
# x <- as.matrix(x,drop=FALSE)
# meanX <- apply(x,2,mean,na.rm=TRUE)
# nX <- apply(x,2,function(xx) by(xx,grp,count.valid))
# 
# if(!is.null(grp)) { if(restrict) {d <- by(x[,1], grp ,function(xx) density(xx,na.rm=TRUE,from=min(xx,na.rm=TRUE),to=max(xx,na.rm=TRUE)))} else {
# d <- by(x[,1], grp ,function(xx) density(xx,na.rm=TRUE)) }} else {
# if(restrict) {d <- density(x[,1],na.rm=TRUE,from=minx,to=maxx)} else {
# d <- density(x[,1],na.rm=TRUE)} }
# d <- unlist(d,recursive=FALSE)
# maxy <- max(d[[1]]$y,d[[2]]$y,na.rm=TRUE)
# plot(NA,ylim=c(0,maxy),xlim=xlim,axes=FALSE,xlab=xlab,ylab=ylab,main=main,pch=pch,...)
#   axis(1,1:nvarg,names,...)
#   axis(2,...)
#   box()
# if(!is.null(scale)) {width <- scale*sqrt(nX[[i]]/tot.n.obs)/max(d[[i]]$y)} else {width <- SCALE/max(d[[i]]$y)}
# polygon(width*c(-d[[i]]$x,d[[i]]$x[rev])+i,c(d[[i]]$y,d[[i]]$y[rev]),density=density,col=col[i],...)
# }


#11/28/17
#Improved 09/02/20

 "densityBy" <- function(x,var=NULL,grp=NULL,data=NULL,freq=FALSE,col=c("blue","red","black"),alpha=.5,adjust=1,ylim=NULL,xlim=NULL,xlab="Variable", ylab="Density",main="Density Plot",legend=NULL) {
 formula <- FALSE
   if(inherits(x, "formula")) {  ps <- fparse(x)
   formula <- TRUE
   if(is.null(data)) stop("You must specify the data if you are using formula input") 
     var <- ps$y
   grp <- ps$x
    x <- data
   
   }
   

 n.grp <- length(table(x[grp]))
 if(length(col) < n.grp) col <- rainbow(n.grp)
if(!is.null(var)) x <- x[,c(var,grp),drop=FALSE]
 x <- char2numeric(x)
 col <- adjustcolor(col,alpha.f =alpha)
 if(missing(main) && freq) main="Frequency Plot"
 if(missing(ylab) && freq) ylab <- "N * density"
 if(missing(xlab) && (length(var) ==1)) xlab <- var
 if(!is.null(grp)) { 
    if(!is.data.frame(grp) && !is.list(grp) && (length(grp) < NROW(x))) grp <- x[,grp,drop=FALSE]
   d <- by(x[,var],grp,function(xx) density(xx, adjust=adjust,na.rm=TRUE))
   quants <-by(x[,var],grp,function(xx) quantile(xx,c(.25,.5,.75),na.rm=TRUE))
   maxiy <- rep(NA,length(d))
   rangex <- matrix(NA, ncol=2,nrow= length(d))
   for(i in 1:length(d) ) {if(freq){ maxiy[i] <- max(d[[i]]$n *d[[i]]$y)} else { 
                   maxiy[i] <- max(d[[i]]$y)}
                   rangex[i,] <- range(d[[i]]$x)
                   }
   maxy <- max(maxiy)
   ranges <- range(rangex)
   if(missing(ylim)) ylim <- c(0,maxy)
   if(missing(xlim)) xlim <- ranges
   plot(d[[1]],ylim=ylim,xlim=xlim,main=main,ylab=ylab,xlab=xlab ) 
   
   for (i in 1:length(d)) {
   if(freq) {scal <- d[[i]]$n
      d[[i]]$y <- d[[i]]$y*scal} 
   polygon(d[[i]] ,col=col[i])

   #find the y value corresponding to x median
   medy <- max(which(d[[i]]$x <= quants[[i]][2]))
   segments(x0 = quants[[i]][2],y0=0,x1=quants[[i]][2],y1=d[[i]]$y[medy],  lty= i)
   }
 } else {
 d <- density(x[,var],na.rm=TRUE,adjust=adjust)
   plot(d,main=main,ylab=ylab,xlab=xlab )
    polygon(d,col=col[1])
 }
 if(!missing(legend)) {
 location <- c("topleft","topright","top","left","right")
  grp.names <- paste(names(grp),names(table(grp)))
 leg.text <- grp.names
 legend(location[legend],legend=leg.text,col=c(1:n.grp),fill=c(1:n.grp), lty=c(1:n.grp))

 }
 invisible(d)
  }