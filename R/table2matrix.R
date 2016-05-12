"table2matrix" <- 
   function(x,labs = NULL) {
  n <- sum(x)
  nrows <- dim(x)[1]
  ncol <- dim(x)[2]
  rowval <- as.numeric(rownames(x))
  colval <- as.numeric(colnames(x))
  
  xm <- matrix(NaN,nrow=n,ncol=2)
  k <- 1
  for (rows in 1:nrows) {
     for (cols in 1:ncol) { 
     case <- x[rows,cols]
     if(case>0) {
     for (cases in 1:case) {
     	 xm[k,1] <- rowval[rows]
     	 xm[k,2] <- colval[cols]
     	 k <- k+1} #end cases
     	 }
       } #end cols
     } #end rows
    if(!is.null(labs)) colnames(xm)<- labs
   return(xm) 
  }
  
 
 "table2df" <- 
   function(x,count=NULL, labs = NULL) {
   if(!is.null(count)) {xm.df <- bigtable2matrix(x,count,labs) } else {
  n <- sum(x)
  nrows <- dim(x)[1]
  ncol <- dim(x)[2]
  rowval <- as.numeric(rownames(x))
  colval <- as.numeric(colnames(x))
  
  xm <- matrix(NaN,nrow=n,ncol=2)
  k <- 1
  for (rows in 1:nrows) {
     for (cols in 1:ncol) { 
     case <- x[rows,cols]
     if(case>0) {
     for (cases in 1:case) {
     	 xm[k,1] <- rowval[rows]
     	 xm[k,2] <- colval[cols]
     	 k <- k+1} #end cases
     	 }
       } #end cols
     } #end rows
    if(!is.null(labs)) colnames(xm)<- labs
    xm.df <- data.frame(xm) }
   return(xm.df) 
  }
  
 "bigtable2matrix" <-
 function(x,count,labs=NULL) {
 n <- dim(x)[1]
 ncol <- dim(x)[2]
 nrows <- sum(count)
 xm <- matrix(NaN,nrow=nrows,ncol=ncol)
 k <- 1
 for (i in 1 :n) {
   for (j in k:(k+count[i]-1)) {
   for (values in 1:ncol) {
    xm[j,values] <- x[i,values] }
    }
    k <- k+count[i]
    }
 if(!is.null(labs)) {colnames(xm) <- labs}
  return(xm) 
  }
 
 
 