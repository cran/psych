"winsor" <-
function(x, trim=.2,na.rm=TRUE) {  
    if(is.vector(x) ) {
   ans <-  wins(x,trim=trim,na.rm=na.rm) } else {
   if (is.matrix(x) | is.data.frame(x)) {ans <- apply(x,2,wins,trim=trim,na.rm=na.rm) } }
   return(ans)
}

"winsor.means" <-
function(x, trim=.2,na.rm=TRUE) {  
    if(is.vector(x) ) {
   ans <-  win.mean(x,trim=trim,na.rm=na.rm) } else {
   if (is.matrix(x) | is.data.frame(x)) {ans <- apply(x,2,win.mean,trim=trim,na.rm=na.rm) } }
   return(ans)
}

"win.mean" <- 
function(x,trim=.2, na.rm=TRUE) {
  if (na.rm) { x <-sort(x[!is.na(x)]) } else {x <- sort(x)}
    ncases <-length(x)
    if ((trim < 0) | (trim>0.5) ) 
        stop("trimming must be reasonable")
     if (trim < .5) {
    numtrim<-trunc(trim*ncases)
   ans <-  (numtrim *x[numtrim+1]+sum(x[(numtrim+1):(ncases-numtrim)])+numtrim*x[ncases-numtrim])/ncases } else {
    ans <- median(x)}
   return(ans)} 
   
"wins" <- 
 function(x,trim=.2, na.rm=TRUE) {
 if (na.rm) { x <-sort(x[!is.na(x)]) } else {x <- sort(x)}
    ncases <-length(x)
    if ((trim < 0) | (trim>0.5) ) 
        stop("trimming must be reasonable")
    numtrim<-trunc(trim*ncases)
    if(trim<.5) { 
   ans <-  c(rep(x[numtrim],numtrim ),x[(numtrim+1):(ncases-numtrim)],rep(x[ncases-numtrim+1],numtrim))} else {
        if((ncases %% 2) >0) {ans <- rep(median(x),ncases) } else {
          ans <- c(rep(x[numtrim],numtrim),rep(x[ncases-numtrim+1],numtrim)) }
          }
   return(ans)} 
   
 
