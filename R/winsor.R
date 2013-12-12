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

"winsor.mean" <-
function(x, trim=.2,na.rm=TRUE) {  
    if(is.vector(x) ) {
   ans <-  win.mean(x,trim=trim,na.rm=na.rm) } else {
   if (is.matrix(x) | is.data.frame(x)) {ans <- apply(x,2,win.mean,trim=trim,na.rm=na.rm) } }
   return(ans)
}

"winsor.var" <-
function(x, trim=.2,na.rm=TRUE) {  
    if(is.vector(x) ) {
   ans <-  win.var(x,trim=trim,na.rm=na.rm) } else {
   if (is.matrix(x) | is.data.frame(x)) {ans <- apply(x,2,win.var,trim=trim,na.rm=na.rm) } }
   return(ans)
}

"winsor.sd" <-
function(x, trim=.2,na.rm=TRUE) {  
    if(is.vector(x) ) {
   ans <-  sqrt(win.var(x,trim=trim,na.rm=na.rm) )} else {
   if (is.matrix(x) | is.data.frame(x)) {ans <- apply(x,2,win.var,trim=trim,na.rm=na.rm) 
    ans <- sqrt(ans) } }
   return(ans)
}

 #added winsor.var and winsor.sd and winsor.mean (to supplement winsor.means) August 28, 2009 following a suggestion by Jim Lemon
  #corrected January 15, 2009 to use the quantile function rather than sorting.
  #suggested by Michael Conklin in correspondence with Karl Healey
  #this preserves the order of the data
"wins" <- 
 function(x,trim=.2, na.rm=TRUE) {
    if ((trim < 0) | (trim>0.5) ) 
        stop("trimming must be reasonable")
      qtrim <- quantile(x,c(trim,.5, 1-trim),na.rm = na.rm)
      xbot <- qtrim[1]
      xtop <- qtrim[3]
       if(trim<.5) { 
      x[x < xbot]  <- xbot
      x[x > xtop] <- xtop} else {x[!is.na(x)] <- qtrim[2]}
     return(x) } 
    
    
"win.mean" <- 
function(x,trim=.2, na.rm=TRUE) {
    if ((trim < 0) | (trim>0.5) ) 
        stop("trimming must be reasonable")
     if (trim < .5) {
   ans <-  mean(wins(x,trim =trim,na.rm=na.rm),na.rm=na.rm)
   return(ans)} else {return(median(x,na.rm=TRUE))} 
   }
   
   
"win.var" <- 
function(x,trim=.2, na.rm=TRUE) {
    if ((trim < 0) | (trim > 0.5) )   {stop("trimming must be reasonable")}
     if (trim < .5) {
   ans <-  var(wins(x,trim =trim,na.rm=na.rm),na.rm=na.rm)
          return(ans)
    } else {return(median(x,na.rm=TRUE))
    }
}