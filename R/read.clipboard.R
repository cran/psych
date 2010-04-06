# a number of functions to read data from the clipboard for both Macs and PCs
"read.clipboard" <-
function(header=TRUE,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) return(read.table(file("clipboard"),header=TRUE,...))
            else return(read.table(file("clipboard"),...)) }
    else {
   if (header) {return(read.table(pipe("pbpaste"),header=TRUE,...))} else {
   return(read.table(pipe("pbpaste"),...))}}
   }

"read.clipboard.csv" <-
function(header=TRUE,sep=',',...) {  #same as read.clipboard(sep=',') 
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) read.clipboard<-read.table(file("clipboard"),header=TRUE,sep,...)
            else read.clipboard<-read.table(file("clipboard"),sep=sep,...) }
    else {
   if (header) read.clipboard<-  read.table(pipe("pbpaste"),header=TRUE,sep,...)
   else read.clipboard<- read.table(pipe("pbpaste") ,sep=sep,...)}
   }
#corrected November 8, 2008 to work with header=FALSE

"read.clipboard.tab" <-
function(header=TRUE,sep='\t',...) {  #same as read.clipboard(sep='\t') 
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) read.clipboard<-read.table(file("clipboard"),header=TRUE,sep,...)
            else read.clipboard<-read.table(file("clipboard"),sep=sep,...) }
    else {
   if (header) read.clipboard<-  read.table(pipe("pbpaste"),header=TRUE,sep,...)
   else read.clipboard<- read.table(pipe("pbpaste") ,sep=sep,...)}
   }
#corrected November 8, 2008 to work with header=FALSE


#adapted from John Fox's read.moments function
"read.clipboard.lower" <-
function( diag = TRUE,names=NULL,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {
           xij <- scan(file("clipboard")) }
    else { xij <- scan(pipe("pbpaste"))}
    m <- length(xij)
    d <- if (diag)    1  else -1
    n <- floor((sqrt(1 + 8 * m) - d)/2) #solve the quadratic for n
    if (m != n * (n + d)/2)   stop("wrong number of elements (cannot make square matrix)")
    X <- diag(n)
    X[upper.tri(X, diag = diag)] <- xij
  diagonal <- diag(X) 
   X <- t(X) + X 
    diag(X) <- diagonal
   if(is.null(names)) names <- paste("V",1:n,sep="")
    rownames(X) <- colnames(X) <- names
    return(X)
   }


"read.clipboard.upper" <-
function( diag = TRUE,names=NULL,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {
           xij <- scan(file("clipboard")) }
    else { xij <- scan(pipe("pbpaste"))}
    m <- length(xij)
    d <- if (diag)    1  else -1
    n <- floor((sqrt(1 + 8 * m) - d)/2) #solve the quadratic for n
    if (m != n * (n + d)/2)   stop("wrong number of elements (cannot make square matrix)")
    X <- diag(n)
    X[lower.tri(X, diag = diag)] <- xij
  diagonal <- diag(X) 
   X <- t(X) + X 
    diag(X) <- diagonal
   if(is.null(names)) names <- paste("V",1:n,sep="")
    rownames(X) <- colnames(X) <- names
    return(X)
   }
   


#added March, 2010 to read fixed width input
"read.clipboard.fwf" <-
function(header=FALSE,widths=rep(1,10),...) {  # 
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) read.clipboard<-read.fwf(file("clipboard"),header=TRUE,widths=widths,...)
            else read.clipboard<-read.fwf(file("clipboard"),widths=widths,...) }
    else {
   if (header) read.clipboard<- read.fwf(pipe("pbpaste"),header=TRUE,widths=widths,...)
   else read.clipboard<- read.fwf(pipe("pbpaste"),widths=widths,...)}
   }