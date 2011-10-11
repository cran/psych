#developed April 24, 2009
#modifed November 14, 2009 to add legends
#completely revised June 20, 2011 to do more reds for less than 0, blues for above 0
#also switched to using layout to make legend clearer
#
"cor.plot" <- 
function(r,colors=TRUE, n=51,main=NULL,zlim=c(-1,1),show.legend=TRUE,labels=NULL,n.legend=10,...){

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
if(colors) { 
    gr <- colorRampPalette(c("red","white","blue")) #added June 20
    colramp  <- gr(n)
      } else {
    colramp <- grey((n:0)/n)}
ord1 <- seq(nvar,1,-1)
if(nvar != nf) {  r <- t(r) }

r <- r[,ord1] #reorder the columns to allow image to work
MAR <- 5
par(mar = c(MAR +max.len,MAR+max.len, 4, .5))

if(show.legend) {   #set it up to do two plots
     layout(matrix(c(1,2),nrow=1),widths=c(.9,.1),heights=c(1,1))
    }

image(r,col=colramp,axes=FALSE,main=main,zlim=zlim)
box()
at1 <- (0:(nf-1))/(nf-1)
at2 <- (0:(nvar-1)) /(nvar-1)
if(max.len>.5) {axis(1,at=at1,labels=rownames(r),las=2,...)
              axis(2,at=at2,labels=colnames(r),las=1,...)} else {
              axis(1,at=at1,labels=rownames(r),...)
              axis(2,at=at2,labels=colnames(r),las=1,...)}
at1 <- (0:(nf-1))/(nf-1)

if(show.legend) {
    leg <- matrix(seq(from=zlim[1],to=zlim[2],by =(zlim[2] - zlim[1])/n),nrow=1)
#screen(2) 
   par(mar=c(MAR,0, 4,3)) 
   image(leg,col=colramp,axes=FALSE,zlim=zlim)  
   at2 <- seq(0,1,1/n.legend)
    labels =seq(zlim[1],zlim[2],(zlim[2]-zlim[1])/(length(at2)-1))
    axis(4,at=at2,labels =labels,las=2,...)
   }  
par(op)  #return the parameters to where we started 
}

