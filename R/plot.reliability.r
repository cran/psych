  #Created June 11-17, 2021
  #fixed 6/20/21 to avoid the problem of null cases
  
"plot.reliability" <- function(x,omega=TRUE,alpha=TRUE,split=TRUE,uni=TRUE, add=FALSE,xlim=NULL,main=NULL,...) {
   if(!is.list(x)  &&  split)  {stop("To show split halfs, you must call reliability with raw=TRUE")}
   rel.object <- x 
   if(!is.list(x)) {result.df<- x} else {result.df<- x$result.df}
  if(split) { x <- rel.object$splits
  #get rid of the empty elements 
     lx <- lapply(x,length)
     
               global.min <- min(unlist(lapply(x[lx>0],min)),result.df[,1])
               global.max <- max(unlist(lapply(x[lx >0],max)),result.df[,3])
               u.max <- max(unlist(lapply(x[lx >0],max)),result.df[,4])
              if(uni) global.max <- max(global.max,u.max)
               nvar <- length(x[lx> 0])
              labels<- names(x[lx>0])
              x <- x[lx>0]
    }  else {
       global.min <- min(result.df[,1])
       global.max <- max(result.df[,3])
       nvar <- NROW(result.df)
       labels <- rownames(result.df)}
       
  if(!add) plot.new()
  
  
    ord <- nvar:1
    labels[ord] <-labels
      linch <- max(strwidth(labels))
       lheight <- par("csi") 
     # linch <- if (!is.null(labels)) max(strwidth(labels, "inch"),label.width, na.rm = TRUE)
        ginch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.1, ginch, na.rm = TRUE) + 0.1)/lheight
   gpos <- 1:nvar
    if (!(is.null(labels) )) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + 0, ginch) + 
            0.1
        par(mai = nmai)
    }
 if(!add ) {if(missing(xlim)) xlim=c(global.min-.05,global.max+.05) #this makes all the plots on the same scale
    plot.window(xlim = xlim, ylim = c(.9,nvar+1.1))}
   
      
      linch <- max(strwidth(labels))
       lheight <- par("csi") 
     # linch <- if (!is.null(labels)) max(strwidth(labels, "inch"),label.width, na.rm = TRUE)
        ginch <- max(strwidth(labels, "inch"), na.rm = TRUE)
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
        mtext(labels, side = 2, line = goffset, at = gpos, adj = 0,  las = 2, ...)
       # if (!is.null(gdata)) {
            abline(h = gpos, lty = "dotted")
   for(i in 1:nvar) {   #show the distributions
    abline(h = i, lty = "dotted")
    if(split) { h <- density(x[[i]],na.rm=TRUE)
      hor  <- h$x
      vert <- h$y/(1 * max(h$y)) + nvar - i + 1
     points(hor,vert,type="l",...) }
   }
if(!add) axis(1)
if(!add) box()

 if(omega){ text(result.df[,1],nvar:1,expression(omega[h]))   #omega_h
            text(result.df[,3],nvar:1,expression(omega[t]))  #omega total
            }
  if(alpha) {text(result.df[,2],nvar:1,expression(alpha))
    }
  if(uni) {text(result.df[,4],nvar:1,"u")}
if(is.null(main)) {if(split) {
        if(omega & alpha & uni ){
                  main=expression(paste("Split half distributions +", ~~~~ omega[h], ~~~~ alpha,~~~~ omega[t], "   + unidim"))} else {
                   if(omega & alpha &!uni) { main=expression(paste("Split half distributions +", ~~~~ omega[h], ~~~~ alpha, ~~~~ omega[t]))} else {
                   
                  if(omega & !alpha & uni) { main=expression(paste("Split half distributions +", ~~~~ omega[h], ~~~~ omega[t], "   + unidim"))} else {
                  if(omega & !alpha) { main=expression(paste("Split half distributions +", ~~~~ omega[h], ~~~~ omega[t]))} else {
                  if(!omega & alpha) {main=expression(paste("Split half distributions +",  ~~~~ alpha,))}
                  }         
                  }
                   }
                  }
                
                  } else {
                  if(omega & alpha & uni) {
                  main=expression(paste(" ", ~~~~ omega[h], ~~~~ alpha,~~~~ omega[t], "   + unidim" ))} else {
                   if(omega & !alpha & uni) { main=expression(paste( ~~~~ omega[h], ~~~~ omega[t], "   + unidim"))} else {
                  if(!omega & alpha) {main=expression(paste(  ~~~~ alpha,))}}
                  }
                  } }
         title(main) 
   }
 