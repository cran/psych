#rewritten Sept 18,2013 to allow much more controll over the points in a two dimensional biplot
#the three or more dimensional case remains the same as before
#The code for the two dimensional case is adapted (heavily) from the stats:::biplot function
#corrected April 9, 2015 to allow multiple biplots in the same window 
#Seriously modified June 25, 2016 to allow better control for labels, as well as just generally cleaner code

"biplot.psych" <-
function(x, labels=NULL,cex=c(.75,1),main="Biplot from fa",hist.col="cyan",xlim.s=c(-3,3),ylim.s=c(-3,3),xlim.f=c(-1,1),ylim.f=c(-1,1),maxpoints=100,adjust=1.2,col,pos, arrow.len = 0.1,pch=16,choose=NULL,cuts=1,cutl=.0,group=NULL,smoother = FALSE,vars=TRUE,...) {
if(is.null(x$scores)) stop("I am sorry, but Biplot requires factor/component scores. \nYou need to  run fa/pca from the raw data")
op <- par()
old.par <- par(no.readonly = TRUE)
on.exit(par(old.par)) 
MAR <- par("mar")
MFROW <- par("mfrow")
title <- main
main <- NULL
 
#fa.poly nests the fa and scores within a list
if(is.list(x$scores)) x$scores <- x$scores$scores  #the case of fa.poly output
if(is.list(x$fa)) x$loadings <- x$fa$loadings  #once again, the case of fa.poly
if(!is.null(choose)) {  #plot just a pair 
     x$scores <- x$scores[,choose,drop=FALSE]
     x$loadings <- x$loadings[,choose,drop=FALSE]
     }
colnames(x$scores) <- colnames(x$loadings)
if((missing(group)) || (is.null(group))) group <- rep(1,nrow(x$scores))

#if(missing(pos)) pos <- NULL
#if(is.null(labels)) {if(nrow(x$scores) > maxpoints ) {labels = rep(".",dim(x$scores)[1] )} else {labels = rep("o",dim(x$scores)[1] )}}

n.dims <- dim(x$loadings)[2]
#this is taken directly from biplot
# if (missing(col)) {
#         col <- par("col")
#         if (!is.numeric(col)) 
#             col <- match(col, palette(), nomatch = 1L)
#         col <- c(col, col + 1L)
#     }
#     else if (length(col) == 1L) 
#         col <- c(col, col)
if(missing(col)) {col <-  c("black","red","blue","#FF0000FF", "#00FF00FF", "#00FFFFFF", "#0000FFFF" ,"#FF00FFFF")}  #rainbow + black, red
 ch.col <- col
 
#here is where we add some plotting controls that are missing from stats:::biplot  
     
if (n.dims == 2) { #just do a one panel graph

  op <- par(pty = "s")
  
  #we no longer resize the margins, but rather adjust where the title goes
  #if (!is.null(main)) op1 <- c(op, par("mar" = MAR + c(0, 0, 1, 0))) #give room for the title  -- 
 # if (!is.null(main)) par("mar" = MAR + c(0, 0, 1, 0))  #give room for the title  -
  #plotone does the work

plotone(x$scores,x$loading,labels=labels,main=main,xlim.s=xlim.s,ylim.s=ylim.s,xlim.f=xlim.f,ylim.f=ylim.f,maxpoints=maxpoints,adjust=adjust,col=col,pos=pos, arrow.len = arrow.len,pch=pch,choose=choose,cuts=cuts,cutl=cutl,group=group,ch.col=ch.col,smoother=smoother,vars=vars,... )
  par(new = TRUE)  #this sets it up so that we can plot on top of a plot 
    dev.hold()
   
    on.exit(dev.flush(), add = FALSE)
  } else {  #the case of 3 or more factors -- we do the equivalent of a pairs plot
	op1 <- par(mfrow=c(n.dims,n.dims), mar=c(2,3,3,2))
	if(nrow(x$scores) > maxpoints) {labels <- rep(".",nrow(x$scores))} else {labels <- rep("o",nrow(x$scores))}
	for (i in 1:n.dims) {
   		for (j in 1:n.dims){ 
  		   if(i==j) {h <- hist(x$scores[,i],freq=FALSE, main=colnames(x$loadings)[i],xlab="",ylab="",col=hist.col)
  		        breaks <- h$breaks; nB <- length(breaks)
   		 tryd <- try( d <- density(x$scores[,i],na.rm=TRUE,bw="nrd",adjust=adjust),silent=TRUE)
     	if(!inherits(tryd,"try-error")) {
    	 lines(d)}
  } else {
   
     #  biplot(x$scores[,c(j,i)],x$loadings[,c(j,i)],xlabs=labels,xlab="",ylab="",cex=cex,xlim=xlim.s,ylim=ylim.s,pch=pch,...)}

      plotone(x$scores[,c(j,i)],x$loadings[,c(j,i)],main=NULL,xlim.s=xlim.s,ylim.s=ylim.s,xlim.f=xlim.f,ylim.f=ylim.f,maxpoints=maxpoints,adjust=adjust,col=col,pos=pos, arrow.len = arrow.len,pch=pch,choose=choose,cuts=cuts,cutl=cutl,group=group,ch.col=ch.col,smoother=smoother,vars=vars,... )}  #work on this way of doing it
           }  }
   }
   #We do not want to reset the margins back to their prior values, because then we lose the ability to add lines to the figure
  # par(old.par)
  # par("mar"=MAR) #this puts the margins back to what they were when we started. Important for multiple biplots
 #  par("mfrow"=MFROW) 
title(title,line=2)
}   #End of biplot.psych

 plotone <- function( scores,loadings,labels=NULL,cex=c(.75,1),main=main,hist.col="cyan",xlim.s=c(-3,3),ylim.s=c(-3,3),xlim.f=c(-1,1),ylim.f=c(-1,1),maxpoints=100,adjust=1.2,col,pos, arrow.len = 0.1,pch=16,choose=NULL,cuts=1,cutl=.0,group=NULL,ch.col=c("black","blue"),smoother=FALSE, vars=TRUE,...)  { 
#There are three different conditions
#1 just show the (unlabeled) points
#2 don't show the points, but put in labels at their locations
#3 show the points with labels near them according to pos

 choice <- "one"
 if(!is.null(labels)) {
     if(missing(pos))  {choice <- "two" } else {choice <- "three"}
     }
if(smoother) choice="smoother"

switch(choice,   #we plot the scores here for scores > abs(cut) on x and y
 "one" = {plot(scores, xlim = xlim.s, ylim = ylim.s, cex=cex[1],main=main,pch=pch[group],bg=ch.col[group],col=col[group],...) } ,
 "two" = {plot(scores,typ='n',  xlim = xlim.s, ylim = ylim.s,cex=cex[1],main=main,pch=pch[group],bg=ch.col[group],col=col[group],...) 
       labels[sqrt((abs(scores[,1])^2  +   abs(scores[,2])^2 ) ) <  cuts] <-   NA
       text(scores,labels=labels,col=ch.col[group],pos=NULL,cex=cex[1])
       },
 "three" = {plot(scores,  xlim = xlim.s, ylim = ylim.s,cex=cex[1],main=main,pch=pch[group],bg=ch.col[group],col=col[group],...) 
        labels[sqrt((abs(scores[,1])^2 + abs(scores[,2])^2)) <  cuts] <-   NA  
       text(scores,labels=labels,pos=pos,cex=cex[1],col=ch.col[group])},
  "smoother" = {smoothScatter(scores, nrpoints=0)
    }
  )
       

     par(new = TRUE)  #this sets it up so that we can plot on top of a plot 
    dev.hold()
   
    on.exit(dev.flush(), add = FALSE)
    plot(loadings, axes = FALSE, type = "n", xlim = xlim.f, ylim = ylim.f,
      xlab = "", ylab = "", col = col[1L], ...)
    labels <- rownames(loadings)
    labels[sqrt(loadings[,1]^2  + loadings[,2]^2) <  cutl] <-   NA
    text(loadings, labels = labels, cex = cex[2L], col = col[2L], ...)
    if(vars) {arrows(0, 0, loadings[, 1L] * 0.8, loadings[, 2L] * 0.8, col = col[2L], 
            length = arrow.len)} else {
            arrows(0, 0, scores[, 1L] * xlim.f/xlim.s, scores[, 2L] * ylim.f/ylim.s, col = col[2L], 
            length = arrow.len)}
    axis(3, col = col[2L], ...)
    axis(4, col = col[2L], ...)
    box(col = col[1L])
  
}    
  

