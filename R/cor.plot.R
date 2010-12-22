#developed April 24, 2009
#modifed November 14, 2009 to add legends
#
"cor.plot" <- 
function(r,colors=FALSE, n=10,main=NULL,zlim=c(-1,1),show.legend=TRUE,labels=NULL,...){

op <- par(no.readonly=TRUE)

if(is.null(main)) {main <- "Correlation plot" }
if(!is.matrix(r) & (!is.data.frame(r))) {if((length(class(r)) > 1) & (class(r)[1] =="psych"))  {if(class(r)[2] =="omega") {r <- r$schmid$sl
nff <- ncol(r)
r <- r[,1:(nff-2)]}  else {r <- r$loadings}
} }
r <- as.matrix(r)
if(min(dim(r)) < 2) {stop ("You need at least two dimensions to make a meaningful plot")}
if(is.null(n)) {n <- dim(r)[2]}
nf <- dim(r)[2]
nvar <- dim(r)[1]
if(is.null(labels)) {
if(is.null(rownames(r))) rownames(r) <- paste("V",1:nvar)
if(is.null(colnames(r))) colnames(r) <- paste("V",1:nf)
} else {rownames(r) <-  colnames(r) <- labels}
 max.len <- max(nchar(rownames(r)))/6
#max.len <- max( strwidth(rownames(r)))
if(is.null(zlim)) {zlim <- range(r)}
if(colors) {gr <- topo.colors(n)
 ord <- n:1
 gr <- gr[ord]} else {
gr <- grey((n:0)/n)}
ord1 <- seq(nvar,1,-1)
if(nvar != nf) {  r <- t(r) }

r <- r[,ord1]
#par.original <- par(c("mar", "las", "mfrow"))$mar
   # on.exit(par(par.original))
#mar <- par.original
mar <- c(5,5)
par(mar = c(mar[1]+max.len,mar[1]+max.len, 4, 4))
if(show.legend) {leg <- seq(from=zlim[1],to=zlim[2],by =(zlim[2] - zlim[1])/(nvar-1))
#leg <- matrix(c(rep(NA,nvar),leg),ncol=2)
#colnames(leg)  <- c("","key")
leg <- matrix(leg,ncol=1)
colnames(leg)  <- "key"
 r <- rbind(r,t(leg))
}
image(r,col=gr,axes=FALSE,main=main,zlim=zlim)
box()
at1 <- (0:(nf-1))/(nf-1)
at2 <- (0:(nvar-1)) /(nvar-1)
if(show.legend) {
   at1 <- (0:(nf))/(nf)
   abline(v=(nf-.5)/nf)
   axis(4,at=at2,labels=round(leg,2),las=2,...)}
  
if(max.len>.5) {axis(1,at=at1,labels=rownames(r),las=2,...)
              axis(2,at=at2,labels=colnames(r),las=1,...)} else {
              axis(1,at=at1,labels=rownames(r),...)
              axis(2,at=at2,labels=colnames(r),las=1,...)}
if (show.legend) {par(op)}  #return the parameters to where we started is the default case
}