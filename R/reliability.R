reliability <- function(keys,items,nfactors=2,split=TRUE,raw=FALSE,plot=FALSE,hist=FALSE) {
 cl <- match.call()
  result <- list()
  splits <- list()
  if(hist) raw <-TRUE 
  if(raw) split <- TRUE
  n.scales <- length(keys)
  
  for (scales in 1:n.scales) {
  scale.key <- keys[[scales]]
 select <- selectFromKeys(scale.key)

   om <- omegah(items[,select], nfactors=nfactors,plot=plot,two.ok=TRUE)
  if(split){ split.half <- suppressWarnings(splitHalf(items[,select],raw=raw,brute=TRUE))
          result[[scales]] <- list(omega_h = om$omega_h, alpha = om$alpha, omega.tot = om$omega.tot,maxrb=split.half$maxrb,minrb=split.half$minrb,
          mean.r=split.half$av.r, med.r <- split.half$med.r )
          if(raw) splits[[scales]] <- split.half$raw} else {      
   result[[scales]] <- list(omega_h = om$omega_h, alpha = om$alpha, omega.tot = om$omega.tot) }
   
  }
  names(result) <- names(keys)
  if(split) {ncol <- 7} else {ncol <- 3}
   result.df <- matrix(unlist(result), ncol=ncol,byrow=TRUE)
  if(split) {   colnames(result.df) <- c("omega_h", "alpha", "omega.tot","max.split","min.split","mean.r", "med.r") } else {
  colnames(result.df) <- c("omega_h", "alpha", "omega.tot")}
  rownames(result.df) <- names(keys)
  if(raw) {
  names(splits) <-names(keys)
 # splits.mat <- matrix(unlist(splits),ncol=length(keys))
 # colnames(splits.mat) <- names(keys)
   class(result.df) <- c("psych","reliability", "matrix")
  result <- list(result.df = result.df,splits= splits, Call = cl)
  if(hist) {multi.hist(splits)}
  class(result) <-  c("psych","reliability")
  return(result) } else {
  
   class(result.df) <- c("psych","reliability", "matrix")
  return(result.df)
  }
  }
  
  #Created June 11-17, 2021
  
 plot.reliability <- function(x,omega=TRUE,alpha=TRUE,split=TRUE, add=FALSE,xlim=NULL,main=NULL,...) {
   if(!is.list(x)  &&  split)  {stop("To show split halfs, you must call reliability with raw=TRUE")}
   rel.object <- x 
   if(!is.list(x)) {result.df<- x} else {result.df<- x$result.df}
  if(split) { x <- rel.object$splits
               global.min <- min(unlist(lapply(x,min)),result.df[,1])
               global.max <- max(unlist(lapply(x,max)),rel.object$result.df[,3])
               nvar <- length(x)
              labels<- names(x)
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
        goffset <- (max(linch + 0.2, ginch, na.rm = TRUE) + 0.1)/lheight
   gpos <- 1:nvar
    if (!(is.null(labels) )) {
        nmai <- par("mai")
        nmai[2L] <- nmai[4L] + max(linch + 0, ginch) + 
            0.1
        par(mai = nmai)
    }
 if(!add ) {if(missing(xlim)) xlim=c(global.min-.05,global.max+.05) #this makes all the plots on the same scale
    plot.window(xlim = xlim, ylim = c(1,nvar+1))}
   
      
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
if(is.null(main)) {if(split) {
        if(omega & alpha){
                  main=expression(paste("Split half distributions +", ~~~~ omega[h], ~~~~ alpha,~~~~ omega[t]))} else {
                  if(omega & !alpha) { main=expression(paste("Split half distributions +", ~~~~ omega[h], ~~~~ omega[t]))} else {
                  if(!omega & alpha) {main=expression(paste("Split half distributions +",  ~~~~ alpha,))}
                  }         
                  
                   }
                  } else {
                  if(omega & alpha) {
                  main=expression(paste(" ", ~~~~ omega[h], ~~~~ alpha,~~~~ omega[t]))} else {
                   if(omega & !alpha) { main=expression(paste( ~~~~ omega[h], ~~~~ omega[t]))} else {
                  if(!omega & alpha) {main=expression(paste(  ~~~~ alpha,))}}
                  }
                  } }
         title(main) 
##title(main)
   }
  

   