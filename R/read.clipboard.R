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
   

#fixed April 30, 2016
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
    
     if(names) { name <- xij[1:n]
               xij <- xij[-c(1:n)] 
               }
            
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


#two useful helper functions
#August, 2016
#modified Jan/April 2017 to include SAS xpt 
"read.file" <- function(file=NULL,header=TRUE,use.value.labels=FALSE,to.data.frame=TRUE,sep=",",widths=NULL,f=NULL,filetype=NULL,...) {
 if(missing(f) && missing(file))  f <- file.choose()
 if(missing(f) && !missing(file)) f <- file
 suffix <- file_ext(f)
 
 if(!missing(filetype)) suffix <- filetype
 if(!missing(widths)) { result <- read.fwf(f,widths,...) 
        message("The fixed width file ", f, "has been loaded.") } else {
 switch(suffix, 
   sav = {result <- read.spss(f,use.value.labels=use.value.labels,to.data.frame=to.data.frame)
         message('Data from the SPSS sav file ', f ,' has been loaded.')},
   csv = {result <- read.table(f,header=header,sep=sep,...)
         message('Data from the .csv file ', f ,' has been loaded.')},
   txt = {result <- read.table(f,header=header,...)
          message('Data from the .txt file ', f , ' has been loaded.') },
  TXT = {result <- read.table(f,header=header,...)
          message('Data from the .TXT file ', f , ' has been loaded.') },
  text = {result <- read.table(f,header=header,...)
          message('Data from the .text file ', f , ' has been loaded.')},
   data =  {result <- read.table(f,header=header,...) 
           message('Data from the .data file ', f , ' has been loaded.')},
    dat =  {result <- read.table(f,header=header,...) 
           message('Data from the .data file ', f , ' has been loaded.')},
     DAT =  {result <- read.table(f,header=header,...) 
           message('Data from the .data file ', f , ' has been loaded.')},
    rds = {result <- readRDS(f,...)
          message('File ',f ,' has been loaded.')},
      R = {result <- dget(f,...)
          message('File ',f ,' has been loaded.')},
      r = {result <- dget(f,...)
          message('File ',f ,' has been loaded.')},
   Rds  = {result <- readRDS(f,...)
          message('File ',f ,' has been loaded.')},
    XPT = { result <- read.xport(f,...)
          message('File ',f ,' has been loaded.')},
    xpt = { result <- read.xport(f,...)
          message('File ',f ,' has been loaded.')},
    #the next options use load rather than read
    #if we return f and it has the same name as the file loaded, this wipes out the file
   Rda = {result <- f   #not helpful if the  
           load(f, .GlobalEnv)
          message('The file(s) in ',f,' have been loaded into your environment.') },
   rda  = {result <- f
          load(f, .GlobalEnv)
          message('The file(s) in ',f,' have been loaded into your environment.') },
    Rdata =   {result <- f
          load(f, .GlobalEnv)
          message('The file(s) in ',f,' have been loaded into your environment.') },
    RData =   {result <- f
          load(f, .GlobalEnv)
          message('The file(s) in ',f,' have been loaded into your environment.') },
     rdata =   {result <- f
          load(f, .GlobalEnv)
          message('The file(s) in ',f,' have been loaded into your environment.') },
        
     SYD =   {result <-  read.systat(f,to.data.frame=to.data.frame )
           message('Data from the systat SYD file ', f ,' has been loaded.')},
     syd =   {result <- read.systat(f,to.data.frame=to.data.frame )
           message('Data from the systat syd file ', f ,' has been loaded.')},
     sys =   {result <-   read.systat(f,to.data.frame=to.data.frame )
           message('Data from the systat sys file ', f ,' has been loaded.')},
    #this section handles (or complains) about jmp and SAS files.
  
    jmp  = {result <- f
          message('I am sorrry.  To read this .jmp file, it must first be saved as either a "txt" or "csv" file. If you insist on using SAS formats, try .xpt or .XPT')},
    sas7bdat = {result <- f
          message('I am sorry.  To read this .sas7bdat file, it must first be saved as either a xpt, or XPT file in SAS, or as a "txt" or "csv" file. ?read.ssd in foreign for help.')},
          
   {message ("I  am sorry. \nI can not tell from the suffix what file type is this.  Rather than try to read it, I will let you specify a better format.")
       }
   )
   }
   return (result)
   }
   
"read.file.spss" <- function(file=NULL,use.value.labels=FALSE,to.data.frame=TRUE,...) {
  if(missing(f) && missing(file))  f <- file.choose()
if(missing(f) &&!missing(file)) f <- file
   result <- read.spss(f,use.value.labels=use.value.labels,to.data.frame=to.data.frame,...)
         message('Data from the SPSS sav file ', f ,' has been loaded.')
return(result)
} 
 "read.file.csv" <- function(file=NULL,header=TRUE,f=NULL,...) {
if(missing(f) && missing(file))  f <- file.choose()
if(missing(f) &&!missing(file)) f <- file
 read.table(f,header=header,sep=",",...) }
 
 "write.file" <- function(x,file=NULL,row.names=FALSE,f=NULL,...) {
	 if(missing(f) && missing(file))  f <- file.choose(TRUE)
	 if(missing(f) &&!missing(file)) f <- file
 	 suffix <- file_ext(f)
	 switch(suffix, 
 		txt = {write.table(x,f, row.names=row.names, ...)},
 		text =   {write.table(x,f,row.names=row.names,...)},
 		csv =  {write.table(x,f,sep=",", row.names=row.names,...) },
 		R =  {dput(x,f,...) },
 		r =  {dput(x,f, ...) },
 		rda = {save(x,file=f,...)},
 		Rda ={save(x,file=f,...)},
 		Rds = {saveRDS(x,f)},
 		rds = {saveRDS(x,f)},
 		write.table(x,f,row.names=row.names) #the default for unspecified types
 		)
   }
 
  "write.file.csv" <- function(x,file=NULL,row.names=FALSE,f=NULL,...) {
if(missing(f) && missing(file))  f <- file.choose(TRUE)
	 if(missing(f) &&!missing(file)) f <- file
 write.table(x,f,sep=",",row.names=row.names,...) }
 
 