"headtail" <- 
function(x, hlength=4,tlength=4,digits=2,ellipsis=TRUE,from=1,to=NULL) {
.Deprecated("headTail", msg = "headtail is deprecated.  Please use the headTail function")
  if(is.data.frame(x) | is.matrix(x) ) { if (is.matrix(x)) x <- data.frame(unclass(x))
  nvar <- dim(x)[2]
  if(is.null(to)) to <- nvar
  dots  <- rep("...",nvar)
   h <- data.frame(head(x[from:to],hlength))
   t <- data.frame(tail(x[from:to],tlength))
 for (i in 1:nvar) {
    if(is.numeric(h[1,i])) {h[i] <- round(h[i],digits)
                         t[i] <- round(t[i],digits)
                         } else {dots[i] <- NA}
     } 
  
  if(ellipsis) { head.tail <- rbind(h,... =  dots,t)} else {head.tail <- rbind(h,t) }
   } else {h <- head(x,hlength)
           t <- tail(x,tlength)
     if(ellipsis) {      head.tail <- rbind(h,"...       ...",t) } else { head.tail <- rbind(h,t)
               head.tail <- as.matrix(head.tail)}}  
  
   return(head.tail)}
 #revised Feb 1, 2010
 #revised August 10, 2011 to work with mixed numeric and non-numeric data
 
 #changing the name of headtail to be camelcase.
 "headTail" <- 
function(x, top=4,bottom=4,from=1,to=NULL,digits=2, hlength=4,tlength=4,ellipsis=TRUE) {
  if(is.data.frame(x) | is.matrix(x) ) { if (is.matrix(x)) x <- data.frame(unclass(x))
  nvar <- dim(x)[2]
  hlength <- top
  tlength <- bottom
   if(is.null(to)) to <- nvar
  dots  <- rep("...",nvar)
   h <- data.frame(head(x[from:to],hlength))
   t <- data.frame(tail(x[from:to],tlength))
 for (i in 1:nvar) {
    if(is.numeric(h[1,i])) {h[i] <- round(h[i],digits)
                         t[i] <- round(t[i],digits)
                         } else {dots[i] <- NA}
     } 
  
  if(ellipsis) { head.tail <- rbind(h,... =  dots,t)} else {head.tail <- rbind(h,t) }
   } else {h <- head(x,hlength)
           t <- tail(x,tlength)
     if(ellipsis) {      head.tail <- rbind(h,"...       ...",t) } else { head.tail <- rbind(h,t)
               head.tail <- as.matrix(head.tail)}}  
  
   return(head.tail)}
 #revised Feb 1, 2010
 #revised August 10, 2011 to work with mixed numeric and non-numeric data
 

 
 topBottom <- 
function (x, top=4,bottom=4,from=1,to=NULL, digits=2, hlength = 4, tlength = 4) 
{
    if (is.data.frame(x) | is.matrix(x)) {
        if (is.matrix(x)) 
            x <- data.frame(unclass(x))
        nvar <- dim(x)[2]
          hlength <- top
         tlength <- bottom  
         if(is.null(to)) to <- nvar
        ellipsis <- rep("...", nvar)
        h <- data.frame(head(x[from:to], hlength))
        t <- data.frame(tail(x[from:to], tlength))
        for (i in 1:nvar) {
            if (is.numeric(h[1, i])) {
                h[i] <- round(h[i], digits)
                t[i] <- round(t[i], digits)
            }
            else {
                ellipsis[i] <- NA
            }
        }
        head.tail <- rbind(h, t)
        head.tail <- as.matrix(head.tail)
    }
    else {
        h <- head(x, hlength)
        t <- tail(x, tlength)
        head.tail <-as.matrix( rbind(h,  t))
    }
    return(head.tail)
}
#added June, 2012

#added April 20, 2017
"quickView" <- function(x,top=8,bottom=8,from=1,to=NULL) {
  if(is.null(to)) to <- NCOL(x)
   View(x[c(1:top,(NROW(x)+1  - bottom):NROW(x)),from:to])
   }