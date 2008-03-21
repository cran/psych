"winsor" <-
function(x, trim=.2,na.rm=TRUE) {  
    if(is.vector(x) ) {
   ans <-  win.mean(x,trim=trim,na.rm=na.rm) } else {
   if (is.matrix(x) | is.data.frame(x)) {ans <- apply(x,2,win.mean,trim=.2,na.rm=na.rm) } }
   return(ans)
}



"win.mean" <- 
function(x,trim=.2, na.rm=TRUE) 
  {if (na.rm) { x <-sort(x[!is.na(x)]) } else {x <- sort(x)}
    ncases <-length(x)
    if ((trim < 0) | (trim>0.5) ) 
        stop("trimming must be reasonable")
    numtrim<-trunc(trim*ncases)
   ans <-  (numtrim *x[numtrim+1]+sum(x[(numtrim+1):(ncases-numtrim)])+numtrim*x[ncases-numtrim])/ncases 
   return(ans)} 