#developed April 24, 2009
#modifed November 14, 2009 to add legends
#completely revised June 20, 2011 to do more reds for less than 0, blues for above 0
#also switched to using layout to make legend clearer
#modified May 5, 2012 to add the keep.par option to allow more control of small graphics
#Corrected feb 22, 2014 to not double plot text (reported by David Condon)
#modified May 12 to allow for selection
#modified March 28, 2016 to add the upper option 
#modified Feb 4, 2017 to label plots and to find correlations by default
#modified April 15, 2017 to allow for non-symmetric matrices
#modified April 15, 2017 to allow more plotting control on the x and y rotation options
#modified November 29th, 2017 to allow for semi-transparency by adjusting alpha
#Finally changed the default to be numbers=TRUE (9/17/19)
#Added the sort option 6/13/23
#Added the ability to rotate x labels xsrt  if xlas =0
#fixed the stars option for the case of lower diagonal so that it works
#Added the n.obs paramter for people who like stars for correlation input  8/6/23


"cor.plot" <- "corPlot" <- 
function(r,numbers=TRUE,colors=TRUE, n=51,main=NULL,zlim=c(-1,1),show.legend=TRUE,
  labels=NULL,
    n.legend=10,keep.par=TRUE,select=NULL,
    pval=NULL,digits=2, trailing=TRUE,
    cuts=c(.001,.01),scale=TRUE,cex,MAR,upper=TRUE,
    diag=TRUE,symmetric=TRUE,stars=FALSE,adjust="holm",
    xaxis =1, xlas=0,ylas=2,ysrt=0, xsrt=0, gr=NULL,alpha =.75,
    min.length=NULL, sort=FALSE, n.obs=NULL, ...){
if(keep.par) op <- par(no.readonly=TRUE)
if(missing(MAR)) MAR <- 5
if(!is.matrix(r) & (!is.data.frame(r))) {if((length(class(r)) > 1) & (inherits(r, "psych")))  {
switch(class(r)[2],
   omega  = {r <- r$schmid$sl
             nff <- ncol(r)
            r <- r[,1:(nff-3)]
           if(is.null(main)) {main <- "Omega plot" }},
  cor.ci ={ pval <- 2*(1-r$ptci)
           r <- r$rho},
  fa = {r <- r$loadings
         if(is.null(main)) {main <- "Factor Loadings plot" }},
  pc = {r <- r$loadings
       if(is.null(main)) {main <- "PCA Loadings plot" }
            },
  principal = {r <- r$loadings
        if(is.null(main)) {main <- "PCA Loadings plot" }}    
  )  #end switch
  }
  } else { if(symmetric & !isCorrelation(r) & (nrow(r) !=ncol(r))) {
    cp <- corr.test(r,adjust=adjust)  #find the corelations
    r <- cp$r
    pval <- cp$p
    if(is.null(main)) {main <- "Correlation plot from data" }
    } else {scale<- FALSE 
      if(is.null(main)) {main <- "Correlation plot"}  
   } 
    if(isCorrelation(r)& !is.null(n.obs) )  pval <- corr.p(r,n.obs)$p 
     }
R <- r <- as.matrix(r)
if(!is.null(select)) r <- r[select,select]
if(sort) {ord <- iclust(r,plot=FALSE)$ord
    r <- r[ord,ord]}
if(min(dim(r)) < 2) {stop ("You need at least two dimensions to make a meaningful plot")}

if(is.null(n)) {n <- dim(r)[2]}
nf <- dim(r)[2]
nvar <- dim(r)[1]

if(!upper) r[col (r) > row(r) ] <- NA   #blank out the upper diagonal
if(!diag) r[col(r) == row(r)] <- NA     #and the diagonal
if(nf == nvar) r <- t(r)  # flip the matrix because grid requires it  but don't flip if a loadings matrix
if(missing(pval)|is.null(pval)) {pval <- matrix(rep(1,nvar*nf),nvar)} else {if (length(pval) != nvar*nf) {
        pr = matrix(0,nvar,nf)
        pr [row(pr) > col(pr)] <- pval # pr[pval[pr]] <- pval
        pr <- pr + t(pr)
        diag(pr) <- 0
        pval <- pr}
       
        if(!stars) {pval <- con2cat(pval,cuts=cuts)
        pval <- (length(cuts)+1-pval)/length(cuts)}
        pval <- t(pval)  #flip these so that the lower off diagonal will be unadjusted
        }

if(is.null(labels)) {
if(is.null(rownames(r))) rownames(r) <- paste("V",1:nvar)
if(is.null(colnames(r))) colnames(r) <- paste("V",1:nf)
} else {rownames(r) <-  colnames(r) <- labels}
if(!is.null(min.length)) {
    rownames(r) <- abbreviate(rownames(r),minlength = min.length)
    colnames(r) <- abbreviate(colnames(r),minlength = min.length)
    }
 max.len <- max(nchar(rownames(r)))/6
#max.len <- max( strwidth(rownames(r)))
if(is.null(zlim)) {zlim <- range(r)}
if(colors) { 
   if(missing(gr))  {gr <- colorRampPalette(c("red","white","blue"))}  #added June 20, 2018? 
   if(max(r,na.rm=TRUE) > 1) {#add a fudge to make all the big ones the same
   maxr <- max(r)
   n1 <- n*(zlim[2]- zlim[1])/(maxr- zlim[1]) 
   colramp <- rep(NA,n)
   n1 <- ceiling(n1)
   colramp[1:(n1+1)] <- gr(n1+1)
   colramp[(n1+1):n] <- colramp[n1+1] 
 zlim[2] <- maxr 
   } else {
    colramp  <- gr(n)}
      } else {
    colramp <- grey((n:0)/n)}
 colramp <- adjustcolor(colramp,alpha.f =alpha)   
if(nvar != nf) {  r <- t(r) }
#if(!is.null(select)) {r <- r[select,select]
#                      pval <- pval[select,select]
#                      nvar <- length(select)
#                      }
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
#if(nf < nvar) {
 			at1 <- (0:(nf-1))/(nf-1)
               at2 <- (0:(nvar-1)) /(nvar-1)
               lab1 <- rownames(r)
               lab2 <- colnames(r)
# } else {
#              at1 <- (0:(nf-1)) /(nf-1)
#              at2 <- (0:(nvar-1)) /(nvar-1)
#               lab1 <- rownames(r)
#               lab2 <- colnames(r)
#              }


#if(nvar != nf) {  r <- t(r) }
if(xaxis == 3) {line <- -.5
 tick <- FALSE} else {line <- NA
 tick <- TRUE}

if(max.len >.5) {axis(2,at=at2,labels=lab2,las=ylas,...)
          if(xsrt==0)  {  axis(xaxis,at=at1,labels=lab1,las=xlas,line=line,tick=tick,...)} else {
            axis(xaxis,at=at1,labels=FALSE,las=xlas,line=line,tick=tick,...)
          text(at1, par("usr")[3]*3, labels = lab1, srt = xsrt, xpd=TRUE,...) }} else {
              axis(2,at=at2,labels=lab2,las=ylas,...)
           if(xsrt==0)  { axis(xaxis,at=at1,labels=lab1,las=xlas,line=line,tick=tick,...)} else {
               axis(xaxis,at=at1,labels=FALSE,las=xlas,line=line,tick=tick,...)
              text(at1, par("usr")[3]*3, labels = lab1, srt = xsrt, xpd=TRUE,...) }
              }
#at1 <- (0:(nf-1))/(nf-1)

if(numbers) {rx <- rep(at1,ncol(r))
            ry <-rep(at2,each=nrow(r))
# rv <- round(r*100)
rv <- round(r,digits) 
 if(trailing) {rv <- sprintf("%.2f", rv)# format with trailing 0rv
              rv[rv=="NA"] <- ""     #drop the NAs  for lower Cors
              }

 if(stars) {#pval <- corr.p1(r,npairs,"none"
           
            symp <- symnum(pval, corr = FALSE,cutpoints = c(0,  .001,.01,.05, 1),
                symbols = c("***","**","*"," "),legend=FALSE)
                rv[(!is.na(rv))] <- paste0(rv[(!is.na(rv))],symp[!is.na(rv)])
                rv <- matrix(rv,ncol=ncol(pval))
                if(!upper) rv[ncol(rv)+diag - row(rv) < col(rv)] <- ""
                 if(missing(cex))  cex = 9/max(nrow(r),ncol(r))
                 text(rx,ry,rv,cex=cex,...)  } else {

 if(missing(cex))  cex = 9/max(nrow(r),ncol(r))
 
 
if(scale) { text(rx,ry,rv,cex=pval*cex,...) } else { text(rx,ry,rv,cex= 1.5*cex,...) }}
 }
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
invisible(R)    #added 11/26/18
}

 #used to find p values for corPlot
 #just an internal function

#notice that the matrix is flipped to do the plotting
#we need to find the lowerright
 "corr.p1" <-
function(r,n,adjust="holm") {
t <- (r*sqrt(n-2))/sqrt(1-r^2)
p <- 2*(1 - pt(abs(t),(n-2)))
p[p>1] <- 1
if (adjust !="none") {
  p[] <- p.adjust(p ,adjust)  #the case of an asymmetric matrix
 }
result <- p
return(result)
}


"reflect" <- function(m) {
  NR <- NROW(m)
  NC <- NCOL(m)
  m[NR+1 - row(m)] <- m[row(m)]
  m[NC +1 - col(m)] <- m[col(m)]}