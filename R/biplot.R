#rewritten Sept 18,2013 to allow much more controll over the points in a two dimensional biplot
#the three or more dimensional case remains the same as before
#The code for the two dimensional case is adapted (heavily) from the stats:::biplot function

"biplot.psych" <-
function(x, labels=NULL,cex=c(.75,1),main="Biplot from fa",hist.col="cyan",xlim.s=c(-3,3),ylim.s=c(-3,3),xlim.f=c(-1,1),ylim.f=c(-1,1),maxpoints=100,adjust=1.2,col,pos, arrow.len = 0.1,pch=16,...) {
if(is.null(x$scores)) stop("Biplot requires factor/component scores:")
op <- par()
old.par <-par(no.readonly = TRUE)
on.exit(par(old.par)) 

#fa.poly nests the fa and scores within a list
if(is.list(x$scores)) x$scores <- x$scores$scores  #the case of fa.poly output
if(is.list(x$fa)) x$loadings <- x$fa$loadings  #once again, the case of fa.poly

colnames(x$scores) <- colnames(x$loadings)

#if(missing(pos)) pos <- NULL
#if(is.null(labels)) {if(nrow(x$scores) > maxpoints ) {labels = rep(".",dim(x$scores)[1] )} else {labels = rep("o",dim(x$scores)[1] )}}

n.dims <- dim(x$loadings)[2]
#this is taken directly from biplot
if (missing(col)) {
        col <- par("col")
        if (!is.numeric(col)) 
            col <- match(col, palette(), nomatch = 1L)
        col <- c(col, col + 1L)
    }
    else if (length(col) == 1L) 
        col <- c(col, col)
 
#here is where we add some plotting controls that are missing from stats:::biplot       
if (n.dims == 2) {
 # biplot(x$scores[,1:2],x$loadings[,1:2],xlabs=labels, cex=cex,main=main,xlim=xlim.s,ylim=ylim.s,...) 
  op <- par(pty = "s")
  if (!is.null(main)) 
        op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0))) #give room for the title

#There are three different conditions
#1 just show the (unlabeled) points
#2 don't show the points, but put in labels at their locations
#3 show the points with labels near them according to pos
 choice <- "one"
 if(!is.null(labels)) {
     if(missing(pos))  {choice <- "two" } else {choice <- "three"}
     }
 
switch(choice,
 "one" = {plot(x$scores, xlim = xlim.s, ylim = ylim.s,  col=col[1],cex=cex[1],main=main,pch=pch, ...) } ,
 "two" = {plot(x$scores,typ='n',  xlim = xlim.s, ylim = ylim.s,  col=col[1],cex=cex[1],main=main,pch=pch, ...)   
       text(x$scores,labels=labels,col=col[1],pos=NULL,cex=cex[1])
       },
 "three" = {plot(x$scores,  xlim = xlim.s, ylim = ylim.s,  col=col[1],cex=cex[1],main=main,pch=pch, ...)   
       text(x$scores,labels=labels,col=col[1],pos=pos,cex=cex[1])
       })

   par(new = TRUE)  #this sets it up so that we can plot on top of a plot 
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    plot(x$loadings, axes = FALSE, type = "n", xlim = xlim.f, ylim = ylim.f,
      xlab = "", ylab = "", col = col[1L], ...)
    
    text(x$loadings, labels = rownames(x$loadings), cex = cex[2L], col = col[2L], ...)
     arrows(0, 0, x$loadings[, 1L] * 0.8, x$loadings[, 2L] * 0.8, col = col[2L], 
            length = arrow.len)
    axis(3, col = col[2L], ...)
    axis(4, col = col[2L], ...)
    box(col = col[1L])
   
    
  
  } else {  #the case of 3 or more factors -- we do the equivalent of a pairs plot
	op1 <- par(mfrow=c(n.dims,n.dims), mar=c(2,3,3,2))
	if(nrow(x$scores) > maxpoints) {labels <- rep(".",nrow(x$scores))} else {labels <- rep("o",nrow(x$scores))}
	for (i in 1:n.dims) {
   		for (j in 1:n.dims){ 
  		   if(i==j) {h <- hist(x$scores[,i],freq=FALSE, main=colnames(x$loadings)[i],xlab="",ylab="",col=hist.col)
  		        breaks <- h$breaks; nB <- length(breaks)
   		 tryd <- try( d <- density(x$scores[,i],na.rm=TRUE,bw="nrd",adjust=adjust),silent=TRUE)
     	if(class(tryd) != "try-error") {
    	 lines(d)}
  } else {
        biplot(x$scores[,c(j,i)],x$loadings[,c(j,i)],xlabs=labels,xlab="",ylab="",cex=cex,xlim=xlim.s,ylim=ylim.s,pch=pch,...)}
           }  }
   }
}


