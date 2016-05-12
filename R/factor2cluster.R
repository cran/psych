"factor2cluster" <-
function (loads,cut=.0) 
{
    
     if (!is.matrix(loads) ) {l <-loads$loadings} else {l <- loads}
  
    l <- as.matrix(l)
    nrows <- dim(l)[1]
    ncols <- dim(l)[2]
   if (ncols ==1) {m1 <- matrix(rep(1,nrows),ncol=1) } else {
    m1 <- matrix(apply(t(apply(l, 1, abs)), 1, which.max), 
        ncol = 1)}
    id <- matrix(c(1:nrows, m1), ncol = 2)  #index row and column
   factor2cluster <- matrix(rep(0, ncols * nrows), ncol = ncols)
   factor2cluster[id] <- sign(l[id])*( (abs(l[id]) >cut)+0)  #only loadings > cut
  rownames(factor2cluster) <- rownames(l)
  colnames(factor2cluster) <- colnames(l)
   nitems <- colSums(abs(factor2cluster))
   for (i in ncols:1) {if (nitems[i]<1) {factor2cluster <- factor2cluster[,-i,drop=FALSE]} }#remove columns with no variables
    return(factor2cluster)
}


