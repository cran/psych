"head.tail" <- 
function(x, hlength=4,tlength=4,digits=2) {
  if(is.data.frame(x) | is.matrix(x) ) { if (is.matrix(x)) x <- data.frame(x)
   ellipsis <- rep("...",dim(x)[2])
   h <- data.frame(head(x,hlength))
   t <- data.frame(tail(x,tlength))
   headtail <- rbind(round(h,digits),'...' = ellipsis,round(t,digits))
   } else {h <- head(x,hlength)
           t <- tail(x,tlength)
           head.tail <- rbind(h,"...       ...",t)}  
   return(head.tail)}
