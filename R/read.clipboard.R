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
#modified October 31, 2010 to be able to read row names as first column
#corrected September 2, 2011 to be able to read row names as first column but without the diagonal
"read.clipboard.lower" <-
function( diag = TRUE,names=FALSE,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {
          con <- file("clipboard")
           } else { con <- pipe("pbpaste" )}
          xij <- scan(con,what="char")
           close(con)
    m <- length(xij)
    d <- if (diag |names)    1  else -1
    n <- floor((sqrt(1 + 8 * m) - d)/2) 
   if(names)  {name <- xij[cumsum(1:n)]     
               xij <- xij[-cumsum(seq(1:n))]
                d <- if (diag )    1  else -1
                n <- floor((sqrt(1 + 8 * (m-n)) - d)/2) }
   xij <- as.numeric(xij)
   X <- diag(n)
    X[upper.tri(X, diag = diag)] <- xij
  diagonal <- diag(X) 
   X <- t(X) + X 
    diag(X) <- diagonal
   if(!names) name <- paste("V",1:n,sep="")
    if(!names) name <- paste("V",1:n,sep="")
    if(names && !diag) {rownames(X) <- colnames(X) <- c(name,paste("V",n,sep="")) } else  {rownames(X) <- colnames(X) <- name }
    return(X)
   }
   


"read.clipboard.upper" <-
function( diag = TRUE,names=FALSE,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
  if (!MAC ) {
          con <- file("clipboard")
           } else { con <- pipe("pbpaste" )}
          xij <- scan(con,what="char")
           close(con)
    m <- length(xij)
    d <- if (diag | names)    1  else -1
    
    n <- floor((sqrt(1 + 8 * m) - d )/2) #solve the quadratic for n
    
     if(names) { name <- xij[c(1,1+cumsum(n:2))]
               xij <- xij[-c(1,1+cumsum(n:2))]} 
   xij <- as.numeric(xij)
    X <- diag(n)
    X[lower.tri(X, diag = diag)] <- xij
  diagonal <- diag(X) 
   X <- t(X) + X 
    diag(X) <- diagonal
   if(!names) name <- paste("V",1:n,sep="")
    rownames(X) <- colnames(X) <- name
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
   
#added May, 2014 to read from https files
"read.https" <-
function(filename,header=TRUE) {
temp <- tempfile()   #create a temporary file
download.file(filename,destfile=temp,method="curl") #copy the https file to temp
result <- read.table(temp,header=header) #now, do the normal read.table command
unlink(temp) #get rid of the temporary file
return(result)}  #give us the result