"read.clipboard" <-
function(header=TRUE,...) {
    MAC<-Sys.info()[1]=="Darwin"    #are we on a Mac using the Darwin system?
   if (!MAC ) {if (header) return(read.table(file("clipboard"),header=TRUE,...))
            else return(read.table(file("clipboard"),...)) }
    else {
   if (header) {return(read.table(pipe("pbpaste"),header=TRUE,...))} else {
   return(read.table(pipe("pbpaste"),...))}}
   }
