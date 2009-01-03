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
