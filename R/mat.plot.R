#developed April 24, 2009
#
"mat.plot" <- 
function(r,colors=FALSE, n=10,main=NULL,zlim=c(0,1)){

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
if(is.null(rownames(r))) rownames(r) <- paste("V",1:nvar)
if(is.null(colnames(r))) colnames(r) <- paste("V",1:nf)
if(is.null(zlim)) {zlim <- range(r)}
if(colors) {gr <- topo.colors(n)
 ord <- n:1
 gr <- gr[ord]} else {
gr <- grey((n:0)/n)}
ord1 <- seq(nvar,1,-1)
if(nvar != nf) {  r <- t(r) }

r <- r[,ord1]
image(r,col=gr,axes=FALSE,main=main,zlim=zlim)
box()
at1 <- (0:(nf-1))/(nf-1)
at2 <- (0:(nvar-1)) /(nvar-1)

axis(1,at=at1,labels=rownames(r))
axis(2,at=at2,labels=colnames(r))

}