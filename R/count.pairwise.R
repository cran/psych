#drastically simplified, March 14, 2009 from two loops to 1 matrix operation
#modified July 2, 2013 to allow not counting the diagonal
 "count.pairwise" <-
function (x, y=NULL,diagonal=TRUE) 
{ .Deprecated("pairwiseCount",msg="count.pairwise is deprecated.  Please use the pairwiseCount function.")
   if(is.null(y)) {n <- t(!is.na(x)) %*% (!is.na(x)) } else { n <- t(!is.na(x)) %*% (!is.na(y)) } 
   if(!diagonal) diag(n) <- NA
       return(n) }
    
pairwiseDescribe <- function(x,y=NULL,diagonal=FALSE,...) {
cp <- pairwiseCount(x,y=y,diagonal=diagonal)
cp <- as.vector(cp[lower.tri(cp,diag=diagonal)])
describe(cp,...)
} 
# replaces count.pairwise 
# slightly improved 
# "pairwiseCount"  <-
# function (x, y=NULL,diagonal=TRUE) 
# { x <- !is.na(x)
#    if(is.null(y)) {n <- t(x) %*% (x) } else { n <- t(x) %*% (!is.na(y)) } 
#    if(!diagonal & is.null(y)) diag(n) <- NA
#        return(n) }
       
 #improvement using crossprod reduces memory  load and speeds it up noticably.
 #For 255K cases and 953 variables, the timings are 114 seconds versus 193 before   
 #timings seem to vary as square of variables (from 100 to 800 at least)   
 "pairwiseCount"  <-
function (x, y=NULL,diagonal=TRUE) 
{ x <- !is.na(x)
   if(is.null(y)) {n <- crossprod(x)} else { n <- crossprod(x,!is.na(y)) } 
   if(!diagonal & is.null(y)) diag(n) <- NA
       return(n) }

    
 #doesn't work yet   
"pairwiseDelete" <- function(x,cut=0) {
if(isCorrelation(x)) {#we have already found correlations, get rid of the least number of NAs possible
    nvar <- ncol(x)
    nmissing <-apply(x,2, function(xx) sum(is.na(xx)))
    max.missing <- max(nmissing,na.rm=TRUE)
    select <- (nmissing == max.missing)
    newx <- x[!select,!select]
      } else { #do a count pairwise  and then get rid of the least number <=  cut
    if(ncol(x)!= nrow(x) ) {pc <- pairwiseCount(x)}
    #OK, it is a matrix of counts 
    nbad <- apply(pc,2, function(xx) sum(xx <  cut)) 
    max.bad <- max(nbad)
     select <- (nbad == max.bad)
     newx <- pc[!select,!select] }
    return(newx)
      }
 
 # created July 30, 2018 to impute correlations     
"pairwiseImpute" <- function(keys,R,fix=FALSE)  {
  cl <- match.call()
  if(!isCorrelation(R)) {message("Correlations found from raw data")
     cp <- pairwiseCount(R)
     R <- cor(R,use="pairwise")
    } else {cp <- NULL
       cpij <- NA}
  n.var <- ncol(R)
  diag(R) <- NA
 if(is.list(keys)) {select <- sub("-","",unlist(keys))
      select <- select[!duplicated(select)] } else {
      keys <- keys2list(keys)
        select <- selectFromKeyslist(colnames(R),keys)
      select <- select[!duplicated(select)]}
      R <- R[select,select,drop=FALSE]
      if(!is.null(cp)) cp <- cp[select,select,drop=FALSE]
      keynames <- colnames(keys)
      keys <- make.keys(R,keys)   
    n.keys <- dim(keys)[2]
    n.items <- dim(keys)[1]
     abskeys <- abs(keys)
     keynames <- colnames(keys)
     num.item <- diag(t(abskeys) %*% abskeys) #how many items in each scale
	n.keys <- ncol(keys)
  
  #brute force
  av.r <- matrix(NA,n.keys,n.keys)
  count.r <- matrix(NA,n.keys,n.keys)
  percent <- matrix(NA,n.keys,n.keys)
  size <- matrix(NA,n.keys,n.keys)
  for(i in 1:n.keys) {
    for( j in 1:n.keys)  {
    selectij <- R[keys[,i] >0,keys[,j] > 0]
    if(!is.null(cp)) {cpij <- cp[keys[,i] >0,keys[,j] > 0]} 
    av.r[i,j] <- mean(selectij,na.rm=TRUE)
    size[i,j] <- mean(cpij,na.rm=TRUE)
    count.r[i,j] <- sum(!is.na(selectij))
    percent[i,j] <- count.r[i,j]/(count.r[i,j] + sum(is.na(selectij)))
    }
    }
    
  colnames(av.r) <- rownames(av.r) <- colnames(keys)  
  colnames(count.r) <- rownames(count.r) <- colnames(keys)
  colnames(percent) <- rownames(percent) <- colnames(keys)
  colnames(size) <- rownames(size) <- colnames(keys) 
  if(is.null(cp)) size <- NULL 
  if(fix) {

   for(i in 1:n.keys) {
    for( j in 1:n.keys)  {
  
    temp <- which(is.na( R[keys[,i] > 0,keys[,j]>0]))
   R[keys[,i] > 0,keys[,j]>0][temp]  <- av.r[i,j]
    }
    }
    diag(R) <- 1
  result <- list(av.r =av.r,count=count.r,percent=percent,imputed=R,Call=cl) 
  } else {
  result <- list(av.r =av.r,count=count.r,percent=percent,size=size,Call=cl)  }
  class(result) <- c("psych","pairwise")
  return(result)
  } 

#modified July 20th to give more useful summaries
"pairwiseReport" <- function(x,y=NULL,cut=0,diagonal=FALSE,...) {
 cl <- match.call()
if(isCorrelation(x)) {#we have already found correlations, flag those that are missing
    report <- which(is.na(x),arr.ind=TRUE)
      } else {
       #do a count pairwise  and then get report those below <=  cut
   if(NCOL(x)!=NROW(x))  {pc <- pairwiseCount(x,y)} else {pc <- x}
   if(is.null(y)) {cp <- as.vector(pc[lower.tri(pc,diag=diagonal)])
    des <- describe(cp,...)} else {des <- NULL}
    
    report <- which((pc <= cut),arr.ind=TRUE)
    
   
      }
   if(length(report >0)) {if(is.null(y)) report <- report[report[,1] > report[,2],,drop=FALSE]
      N <- rep(0,dim(report)[1])
      for (i in 1:dim(report)[1] ){N[i] <- pc[report[i,1],report[i,2]]}
      report <- cbind(report,N)
      if(is.null(y) ) {
       df <- data.frame(rows=colnames(x)[report[,1]],cols=colnames(x)[report[,2]],N=N) } else { #although this can be pretty big
       df <- data.frame(rows =colnames(x)[report[,1]],cols =colnames(y)[report[,2]], N=N) }
       
       result <- list(description=des, rows=table(df$rows),cols=table(df$cols),cut = cut, df=df, Call=cl)} else {
      result <- list(description=des,rows=NULL,cols=NULL,cut=cut,df=NULL,Call=cl)}
      class(result) <- c("psych","pairwiseCounts")
      return(result)
      }
      
"pairwisePlot" <- function(x,y=NULL,upper=TRUE,diagonal=TRUE,labels=TRUE,show.legend=TRUE,n.legend=10,colors=FALSE,gr=NULL,minlength=6,xlas=1,ylas=2,main="Relative Frequencies",count=TRUE,...) {
    if(count){ r <- pairwiseCount(x=x,y=y,diagonal=diagonal)} else {r <- x}
     if(!upper) r[col (r) < row(r) ] <- NA   #blank out the upper diagonal
if(!diagonal) r[col(r) == row(r)] <- NA
     nvar <- NROW(r)
     nf <- NCOL(r)
     MAR <- 5
     if(is.null(colnames(x))) colnames(r) <- paste0("V",1:nf)
     if(is.null(rownames(x))) rownames(r) <- paste0("V",1:nvar)
     if(!labels) {minlength <- NULL
        max.len <- 1} else {
    max.len <- min(max(nchar(rownames(r)))/6,minlength)}
    r <- r/max(r,na.rm=TRUE)
    zlim <- c(0,1)
if(colors) { 
   if(missing(gr))  {gr <- colorRampPalette(c("red","white","blue"))}  
    colramp  <- gr(n.legend)
      } else {
    colramp <- grey((n.legend:0)/n.legend)}
 #colramp <- adjustcolor(colramp,alpha.f =alpha)   
if(nvar != nf) {  r <- t(r) }
#if(!is.null(select)) {r <- r[select,select]
#                      pval <- pval[select,select]
#                      nvar <- length(select)
#                      }
#reverse the order of the columns (if square)
 ord1 <- seq(nvar,1,-1)  
 
#if(nf == nvar) {
r <- r[,ord1] 
#}


 #reorder the columns to allow image to work
#MAR <- 5
par(mar = c(MAR +max.len,MAR+max.len, 4, .5))
 line <- NA
 tick <- TRUE


if(show.legend) {   #set it up to do two plots
     layout(matrix(c(1,2),nrow=1),widths=c(.9,.1),heights=c(1,1))
    }
    

image(r,col=colramp,axes=FALSE,main=main,zlim=zlim)
box()
if(labels) {
if(!is.null(minlength)) {
    rownames(r) <- abbreviate(rownames(r),minlength = minlength)
    colnames(r) <- abbreviate(colnames(r),minlength = minlength)
    
 max.len <- max(nchar(rownames(r)))/6}
 at1 <- (0:(nf-1))/(nf-1)
               at2 <- (0:(nvar-1)) /(nvar-1)
               lab1 <- rownames(r)
               lab2 <- colnames(r)
               
axis(2,at=at2,labels=lab2,las=ylas,...)
tick <- FALSE
 axis(1,at=at1,labels=lab1,las=xlas,line=line,tick=tick,...)
 }
 #screen 2
  leg <- matrix(seq(from=zlim[1],to=zlim[2],by =(zlim[2] - zlim[1])/n.legend),nrow=1)
  par(mar=c(MAR,0, 4,3)) 
   image(leg,col=colramp,axes=FALSE,zlim=zlim)  
   at2 <- seq(0,1,1/n.legend)
    labels =seq(zlim[1],zlim[2],(zlim[2]-zlim[1])/(length(at2)-1))
    axis(4,at=at2,labels =labels,las=2,...)  

#put them back in the logical order
ord1 <- seq(nvar,1,-1)  
r <- r[,ord1 ]
invisible(r)
}