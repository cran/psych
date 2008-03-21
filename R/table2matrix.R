"table2matrix" <- 
   function(x,labs = NULL) {
  n <- sum(x)
  nrows <- dim(x)[1]
  ncol <- dim(x)[2]
  rowval <- as.numeric(rownames(x))
  colval <- as.numeric(colnames(x))
  
  xm <- matrix(NA,nrow=n,ncol=2)
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
   function(x,labs = NULL) {
  n <- sum(x)
  nrows <- dim(x)[1]
  ncol <- dim(x)[2]
  rowval <- as.numeric(rownames(x))
  colval <- as.numeric(colnames(x))
  
  xm <- matrix(NA,nrow=n,ncol=2)
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
    xm.df <- data.frame(xm)
   return(xm.df) 
  }
  
   
 