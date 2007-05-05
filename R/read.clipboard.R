"read.clipboard" <-
function(header=TRUE,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) read.clipboard<-read.table(file("clipboard"),header=TRUE,...)
            else read.clipboard<-read.table(file("clipboard"),...) }
    else {
   if (header) read.clipboard<-  read.table(pipe("pbpaste"),header=TRUE,...)
   else read.clipboard<- read.table(pipe("pbpaste"),...)}
   }

