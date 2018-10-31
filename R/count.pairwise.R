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
#replaces count.pairwise 
"pairwiseCount"  <-
function (x, y=NULL,diagonal=TRUE) 
{
   if(is.null(y)) {n <- t(!is.na(x)) %*% (!is.na(x)) } else { n <- t(!is.na(x)) %*% (!is.na(y)) } 
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

"pairwiseReport" <- function(x,cut=0) {
if(isCorrelation(x)) {#we have already found correlations, flag those that are missing
    report <- which(is.na(x),arr.ind=TRUE)
      } else { #do a count pairwise  and then get report those below <=  cut
    pc <- pairwiseCount(x)
    report <- which((pc <= cut),arr.ind=TRUE)
      }
      report <- report[report[,1] > report[,2],,drop=FALSE]
      df <- data.frame(rows=rownames(report),cols=colnames(x)[report[,2]],N= pc[report[,2]])
      return(df)
      }
