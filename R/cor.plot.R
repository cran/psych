#developed April 24, 2009
#modifed November 14, 2009 to add legends
#completely revised June 20, 2011 to do more reds for less than 0, blues for above 0
#also switched to using layout to make legend clearer
#modified May 5, 2012 to add the keep.par option to allow more control of small graphics
#Corrected feb 22, 2014 to not double plot text (reported by David Condon)
#modified May 12 to allow for selection
"cor.plot" <- 
function(r,numbers=FALSE,colors=TRUE, n=51,main=NULL,zlim=c(-1,1),show.legend=TRUE,labels=NULL,n.legend=10,keep.par=TRUE,select=NULL,pval=NULL,cuts=c(.001,.01),cex,MAR,...){
if(keep.par) op <- par(no.readonly=TRUE)
if(missing(MAR)) MAR <- 5
if(is.null(main)) {main <- "Correlation plot" }
if(!is.matrix(r) & (!is.data.frame(r))) {if((length(class(r)) > 1) & (class(r)[1] =="psych"))  {if(class(r)[2] =="omega") {r <- r$schmid$sl
nff <- ncol(r)
r <- r[,1:(nff-3)]}  else {r <- r$loadings}  #fixed 12/31/14 to match revised omega
} }
r <- as.matrix(r)
if(min(dim(r)) < 2) {stop ("You need at least two dimensions to make a meaningful plot")}


if(is.null(n)) {n <- dim(r)[2]}
nf <- dim(r)[2]
nvar <- dim(r)[1]
if(nf == nvar) r <- t(r)  # flip the matrix because grid requires it  but don't flip if a loadings matrix
if(missing(pval)) {pval <- matrix(rep(1,nvar*nf),nvar)} else {if (length(pval) != nvar*nf) {
        pr = matrix(0,nvar,nf)
        pr[lower.tri(pr)] <- pval
        pr <- pr + t(pr)
        diag(pr) <- 0
        pval <- pr}
       
        pval <- con2cat(pval,cuts=cuts)
        pval <- (length(cuts)+1-pval)/length(cuts)
        pval <- t(pval)  #flip these so that the lower off diagonal will be unadjusted
        }

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

if(nvar != nf) {  r <- t(r) }
if(!is.null(select)) {r <- r[select,select]
                      pval <- pval[select,select]
                      nvar <- length(select)
                      }
#reverse the order of the columns (if square)
 ord1 <- seq(nvar,1,-1)  
 
if(nf == nvar) {r <- r[,ord1] 
 pval <- pval[,ord1]} else {r <- r[,ord1] 
 pval <- t(pval[ord1,])}
 #reorder the columns to allow image to work
#MAR <- 5
par(mar = c(MAR +max.len,MAR+max.len, 4, .5))

if(show.legend) {   #set it up to do two plots
     layout(matrix(c(1,2),nrow=1),widths=c(.9,.1),heights=c(1,1))
    }

image(r,col=colramp,axes=FALSE,main=main,zlim=zlim)
box()
if(nf < nvar) {at1 <- (0:(nf-1))/(nf-1)} else {at1 <- (0:(nvar-1)) /(nvar-1)}
at2 <- (0:(nvar-1)) /(nvar-1)

#if(nvar != nf) {  r <- t(r) }
if(max.len>.5) {axis(2,at=at2,labels=colnames(r),las=1,...)
              axis(1,at=at1,labels=rownames(r),las=2,...)} else {
              axis(2,at=at2,labels=colnames(r),...)
              axis(1,at=at1,labels=rownames(r),las=1,...)}
#at1 <- (0:(nf-1))/(nf-1)

if(numbers) {rx <- rep(at1,ncol(r))
 ry <-rep(at2,each=nrow(r))
# rv <- round(r*100)
rv <- round(r,2)
 if(missing(cex))  cex = 10/max(nrow(r),ncol(r))
 text(rx,ry,rv,cex=pval*cex,...)}
if(show.legend) {
    leg <- matrix(seq(from=zlim[1],to=zlim[2],by =(zlim[2] - zlim[1])/n),nrow=1)
#screen(2) 
   par(mar=c(MAR,0, 4,3)) 
   image(leg,col=colramp,axes=FALSE,zlim=zlim)  
   at2 <- seq(0,1,1/n.legend)
    labels =seq(zlim[1],zlim[2],(zlim[2]-zlim[1])/(length(at2)-1))
    axis(4,at=at2,labels =labels,las=2,...)
   }  
if(keep.par) par(op)  #return the parameters to where we started 
}





