#Added mar and changes main Sept 23, 2018
#slightly revised June 12, 2021

"multi.hist" <-
function(x,nrow=NULL,ncol=NULL,density=TRUE,freq=FALSE,bcol="white",dcol=c("black","black"),dlty=c("dashed","dotted"),main=NULL,mar=c(2,1,1,1), breaks=21,global=TRUE,...) {
#if((!is.matrix(x)) & (!is.data.frame(x))) {nvar <- 1
#    x <- as.matrix(x,ncol=1) } else {
    if(is.list(x) & !is.data.frame(x)) {nvar <- length(x)
    xlabs <- main  <-names(x)
    xx <- unlist(x)
     global.min <- min(xx,na.rm=TRUE)
 global.max <- max(xx,na.rm=TRUE)} else {
    x <- as.data.frame(x)
    xlabs <-colnames(x) 
nvar <- NCOL(x)  #number of variables

 global.min <- min(x,na.rm=TRUE)
 global.max <- max(x,na.rm=TRUE)
 }
 if(global) xlim=c(global.min-.05,global.max+.05) #this makes all the plots on the same scale
   old.par <- par(no.readonly = TRUE)
if((is.null(main)) & nvar==1) main <- "Histogram, Density, and Normal Fit"
if (length(dcol)<2) dcol <- c(dcol,dcol) 
  
     if(is.null(main)) {main <- c(colnames(x)) } else {main <- rep(main,nvar)}
     nsize=ceiling(sqrt(nvar))   #size of graphic
     if(is.null(nrow) ) {nrow <- nsize} else {ncol <- nvar/nrow}
     if(is.null(ncol)) {ncol  <- ceiling(nvar/nsize )} else {nrow <- nvar/ncol}
     
    par(mfrow=c(nrow,ncol))       #set new graphic parameters
    par(mar=mar)


     for (i in 1:nvar) {
        
    	 
    	 if(is.data.frame(x) ) {xi <-x[,i]} else {xi <-x[[i]] }
    	 if(!global) {xlim = c(min(xi, na.rm=TRUE),max(xi,na.rm=TRUE))}
    	if(density) {histo.density(xi,xlab=xlabs[i],main=main[i],freq=freq,bcol,dcol=dcol,dlty=dlty,breaks=breaks,xlim=xlim,...)} else {
    				hist(xi,main=main[i],xlab=xlabs[i],freq=freq,bcol,dcol=dcol,dlty=dlty,breaks=breaks,xlim=xlim,...)}
    
    }  #draw the histograms for each variable
     on.exit(par(old.par))   #set the graphic parameters back to the original
     }
   #  }

"histo.density" <- 
function(x,main="Histogram, Density, and Normal Fit",freq=FALSE,xlab=NULL,bcol="white",dcol=c("black","black"),dlty=c("dashed","dotted"),breaks=21,...) {
h <-  hist(x,plot=FALSE,breaks=breaks)
m1 <- mean(x,na.rm=TRUE)
s1 <- sd(x,na.rm=TRUE)
d <- density(x,na.rm=TRUE)

 if(freq) {ymax <- max(h$count)} else {ymax <- max(h$density)}
 dmax <- max(d$y)
ymax <- max(ymax,dmax)
plot(h,freq=freq,ylim=c(0,ymax*1.2),main=main,xlab=xlab,col=bcol,ylab=NULL,...)
if(!freq) {lines(d,lty=dlty[1],col=dcol[1],...)
curve(dnorm(x,m1,s1),add=TRUE,lty=dlty[2],col=dcol[2],...)} else {
lines(d$x,lty=dlty[1],col=dcol[1],...)}
      
}



"histBy" <- function(x,var,group,data=NULL,density=TRUE,alpha=.5,breaks=21,col,xlab,main="Histograms by group",freq=FALSE,...) {
formula <- FALSE
   if(inherits(x, "formula")) {  ps <- fparse(x)
   formula <- TRUE
   if(is.null(data)) stop("You must specify the data if you are using formula input") 
     x <- data
   group <- ps$x
   var <- ps$y
   }
   old.par <- par(no.readonly = TRUE)
if(missing(xlab)) xlab = var
if(missing(group)) {
   if(missing(col)) col12 <- col2rgb("blue",TRUE)/255
    col <- rgb(col12[1],col12[2],col12[3],alpha)
   hist(x[,var],xlab=xlab,main=main,breaks=breaks,freq=freq,col=col,...)
  # hist(x[,var],xlab=xlab,main=main,breaks=breaks,col=col,...)
      d <- density(x[,var],na.rm=TRUE)
      if(density) lines(d)
 } else { #the normal case 
gr <- x[group]                
grp<- table(gr)

if(missing(col)) col <- rainbow(length(grp))
col12 <- col2rgb(col,TRUE)/255
col <- rgb(col12[1,],col12[2,],col12[3,],alpha)

xlim=range(x[var],na.rm=TRUE)
test <- hist(x[,var],breaks=breaks,plot=FALSE)
breaks <- test$breaks

grp <- names(grp)
d <- density(x[(gr==grp[1]),var],na.rm=TRUE)
 hist(x[(gr==grp[1]),var],xlim=xlim,col=col[1],breaks=breaks,freq=freq,xlab=xlab,main=main,...)

if(density) lines(d)
for(i in (2:length(grp))) {
# hist(x[(gr==grp[i]),var],xlim=xlim,col=col[i],freq=FALSE,breaks=breaks,add=TRUE,...)
hist(x[(gr==grp[i]),var],xlim=xlim,col=col[i],breaks=breaks,add=TRUE,freq=freq,...)
d <- density(x[(gr==grp[i]),var],na.rm=TRUE)
if(density) lines(d)
}}
   on.exit(par(old.par)) 
}


#moved to (and expanded) as plot.reliabiity
#  plot.reliability <- function(x,ncol=NULL,nrow=NULL,...) {
#    global.min <- min(x,na.rm=TRUE)
#   global.max <- max(x,na.rm=TRUE)
#   nvar <-nrow(x)
#   
#     old.par <- par(no.readonly = TRUE)
#     nsize=ceiling(sqrt(nvar))
#    if(is.null(nrow) ) {nrow <- nsize} else {ncol <- nvar/nrow}
#      if(is.null(ncol)) {ncol  <- ceiling(nvar/nsize )} else {nrow <- nvar/ncol}
#      
#     par(mfrow=c(nrow,ncol))       #set new graphic parameters
#     #par(mar=mar)
#   
#   for (i in 1:nvar){
#   
#   mp =barplot(x[i,],...)
# 
#    }  
#    on.exit(par(old.par))    
#   }