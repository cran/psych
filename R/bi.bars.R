bi.bars <- function(x,var=NULL, grp=NULL,horiz,color,label=NULL,zero=FALSE,xlab,ylab,...) {
if(missing(horiz)) horiz <- TRUE
if(missing(color)) color <- c("blue","red")

#new way is cleaner
 if(!is.null(var) & (length(var)==1)) {if(missing(ylab) & (length(var) ==1)) {ylab <- var}
          x <- x[,c(var,grp) , drop = FALSE]
          if(is.null(grp)) {stop("I am stopping. The grouping variable was not specified")}
          grp <- x[,grp,drop=FALSE]
          
          x <- as.numeric( x[,var])
          } else {grp <- var}  #the old way

if(horiz) {

  if(missing(xlab)) xlab <- "Frequency"
  if(missing(ylab)) ylab <- paste0(levels(grp))
  } else {
      if(missing(ylab)) ylab <- "Frequency"
      if(missing(xlab)) xlab <- paste0(levels(grp))}

groups <- table(grp)
max.val <- max(x,na.rm=TRUE)
min.val <- min(x,na.rm=TRUE)
#gr1 <- as.numeric(names(groups)[1])
#gr2 <- as.numeric(names(groups)[2])
gr1 <- (names(groups)[1])
gr2 <- (names(groups)[2])
g1 <- subset(x,grp==gr1)
g2 <- subset(x,grp==gr2)
t1 <-  tabulate(g1-min.val*zero,nbins=(max.val-min.val+1))
t2 <-  tabulate(g2-min.val*zero,nbins=(max.val-min.val+1)) 
#t1 <- table(g1-min.val*zero)
#t2 <- table(g2-min.val*zero)


m1 <- max(t1,t2)
m2 <- max(t1,t2)
xlim <- c(-m1,m2)*1.04
if(horiz) {
#t1 <- t1[t1 > 0]
xloc <- barplot(-t1,xlim=xlim,col=color[1],horiz=horiz,xlab=xlab,ylab=ylab,axes=FALSE,axisnames=FALSE,...)
barplot(t2,add=TRUE,col=color[2],horiz=horiz,axes=FALSE,axisnames=FALSE,...)
box()
if((max.val - min.val) < 10) {
if(is.null(label)) {axis(2,at=xloc+min.val*zero,labels=min.val:max.val,...)} else {
                    axis(2,at=xloc+min.val*zero,labels=label,las=2,...)}} else {
                    at <- axTicks(2,usr=c(min.val,max.val))
                    axis(2,at=at,labels=at + min.val*zero,las=2,...)}
atv <- axTicks(1)
axis(1,at=atv,labels=abs(atv),...)} else { #the vertical case
ylim <- c(-m1,m2)*1.04
xloc <- barplot(-t1,ylim=ylim,col=color[1],horiz=horiz,xlab=xlab,ylab=ylab,axes=FALSE,...)
barplot(t2 ,add=TRUE,col=color[2],horiz=horiz,axes=FALSE,...)
box()
atv <- axTicks(2)
axis(2,at=atv,labels=abs(atv),las=2,...)
if((max.val - min.val) < 10) {
if(is.null(label)) {axis(1,at=xloc,labels=min.val:max.val,...)} else {
     axis(1,at=xloc,labels=label,...)}
    } else {
    at <- axTicks(1,usr=c(min.val,max.val))
                    axis(1,at=at,labels=at+min.val*zero,...)
                    }}
}
#modified December 2, 2017 to be compatible with densityBy and violinBy syntax
