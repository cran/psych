bi.bars <- function(x,grp,horiz,color,label=NULL,...) {
if(missing(horiz)) horiz <- TRUE
if(missing(color)) color <- c("blue","red")
groups <- table(grp)
max.val <- max(x,na.rm=TRUE)
min.val <- min(x,na.rm=TRUE)
gr1 <- as.numeric(names(groups)[1])
gr2 <- as.numeric(names(groups)[2])
g1 <- subset(x,grp==gr1)
g2 <- subset(x,grp==gr2)
t1 <-  tabulate(g1,nbins=(max.val-min.val+1))
t2 <-  tabulate(g2,nbins=(max.val-min.val+1)) 
m1 <- max(t1,t2)
m2 <- max(t2)
xlim <- c(-m1,m2)*1.04
if(horiz) {
xloc <- barplot(-t1,xlim=xlim,col=color[1],horiz=horiz,xlab="count",axes=FALSE,...)
barplot(t2,add=TRUE,col=color[2],horiz=horiz,axes=FALSE,...)
box()
if((max.val - min.val) < 10) {
if(is.null(label)) {axis(2,at=xloc,labels=min.val:max.val,...)} else {
                    axis(2,at=xloc,labels=label,las=2,...)}} else {axis(2)}
atv <- axTicks(1)
axis(1,at=atv,labels=abs(atv),...)} else { #the vertical case
ylim <- c(-m1,m2)*1.04
xloc <- barplot(-t1,ylim=ylim,col=color[1],horiz=horiz,ylab="count",axes=FALSE,...)
barplot(t2,add=TRUE,col=color[2],horiz=horiz,axes=FALSE,...)
box()
atv <- axTicks(2)
axis(2,at=atv,labels=abs(atv),las=2,...)
if((max.val - min.val) < 10) {
if(is.null(label)) {axis(1,at=xloc,labels=min.val:max.val,...)} else {
     axis(1,at=xloc,labels=label,...)}
    } else {axis(1)}}
}