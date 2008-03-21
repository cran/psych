"headtail" <- 
function(x, hlength=6,tlength=6,digits=2) {
  if(!is.data.frame(x)) x <- data.frame(x)
   ellipsis <- rep("...",dim(x)[2])
   h <- data.frame(head(x,hlength))
   t <- data.frame(tail(x,tlength))
   headtail <- rbind(round(h,digits),'...' = ellipsis,round(t,digits))
   return(headtail)}