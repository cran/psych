bi.bars <- function(x,grp,horiz,color,...) {
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
xlim <- c(-m1,m2)
if(horiz) {
xloc <- barplot(-t1,xlim=xlim,col=color[1],horiz=horiz,xlab="count",...)
barplot(t2,add=TRUE,col=color[2],horiz=horiz,...)
box()
if((max.val - min.val) < 10) {
axis(2,at=xloc,labels=min.val:max.val,...)} else {axis(2)}
axis(1,...)} else {
ylim <- c(-m1,m2)
xloc <- barplot(-t1,ylim=ylim,col=color[1],horiz=horiz,ylab="count",...)
barplot(t2,add=TRUE,col=color[2],horiz=horiz,...)
box()
axis(2,...)
if((max.val - min.val) < 10) {axis(1,at=xloc,labels=min.val:max.val,...)} else {axis(1)}}
}