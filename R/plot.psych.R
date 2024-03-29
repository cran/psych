"plot.psych" <-
function(x,labels=NULL,...)
  {

   result <- NULL
  
     vss <- iclust <- omega <- fa <-  irt.fa <- irt.poly <-  principal <- parallel <- set.cor <- residuals <- reliability <-    FALSE 
  


   if(length(class(x)) > 1)  {
   #This next section was added December 7, 2019 to change from class(x)[2] to inherits(x, ...)
   		names <- cs(irt.fa,irt.poly,vss,iclust,fa, principal,omega,parallel,set.cor,residuals,reliability)
        value <- inherits(x,names,which=TRUE)   # value <- class(x)[2]
         if(any(value > 1) ) { value <- names[which(value > 0)]} else {value <- "None"}
    
     } else {stop ("I am sorry, this is not an object I know how to plot")}
 
switch(value,

   vss = {
   	    
  		n=dim(x)
  		symb=c(49,50,51,52)              #plotting sym
		plot(x$cfit.1,ylim=c(0,1),type="b",ylab="Very Simple Structure Fit",xlab="Number of Factors",pch=49)
    	
		x$cfit.3<- x$vss.stats$cfit.3
		x$cfit.4<- x$vss.stats$cfit.4
	
  title <- x$title
		title(main=title)
		x$cfit.2[1]<-NA
		x$cfit.3[1]<-NA
		x$cfit.3[2]<-NA
		x$cfit.4[1]<-NA
		x$cfit.4[2]<-NA
		x$cfit.4[3]<-NA
		lines(x$cfit.2)
		points(x$cfit.2,pch=50)
		lines(x$cfit.3)
		points(x$cfit.3,pch=symb[3])
		lines(x$cfit.4)
		points(x$cfit.4,pch=symb[4])
 },	
 
iclust = {
 op <- par(no.readonly = TRUE) # the whole list of settable par's.
  cut <- 0
# if(iclust) {
  load <- x$loadings
          cat("Use ICLUST.diagram to see the  hierarchical structure\n") 
#          } else {load <- x$schmid$orthog 
#         cat("Use omega.diagram to see the  hierarchical structure\n") }
  nc <- dim(load)[2]
  nvar <- dim(load)[1]
ch.col=c("black","blue","red","gray","black","blue","red","gray")

	cluster <- rep(nc+1,nvar)
	cluster <- apply( abs(load)  ,1,which.max)
	cluster[(apply(abs(load),1,max) < cut)] <- nc+1
	
if (nc > 2 ) {
 pairs(load,pch = cluster+19,col=ch.col[cluster],bg=ch.col[cluster]) }
 else {
 plot(load,pch = cluster+20,col=ch.col[cluster],bg=ch.col[cluster],...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc <3) text(load,labels,pos=1)
 par(op) },
 
 omega = {
 op <- par(no.readonly = TRUE) # the whole list of settable par's.
  cut <- 0
 if(iclust) { load <- x$loadings
          cat("Use ICLUST.diagram to see the  hierarchical structure\n") } else {load <- x$schmid$orthog 
         cat("Use omega.diagram to see the  hierarchical structure\n") }
  nc <- dim(load)[2]
  nvar <- dim(load)[1]
ch.col=c("black","blue","red","gray","black","blue","red","gray")

	cluster <- rep(nc+1,nvar)
	cluster <- apply( abs(load)  ,1,which.max)
	cluster[(apply(abs(load),1,max) < cut)] <- nc+1
	
if (nc > 2 ) {
 pairs(load,pch = cluster+19,cex=1.5,col=ch.col[cluster],bg=ch.col[cluster]) }
 else {
 plot(load,pch = cluster+20,col=ch.col[cluster],bg=ch.col[cluster],...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc <3) text(load,labels,pos=1)
 par(op) },
 
 
set.cor ={ plot(x$cancor2,typ="b",ylab="Squared Canonical Correlation",xlab="Canonical variate",main="Scree of canonical correlations" ,ylim=c(0,1),...)},
lmCor ={ plot(x$cancor2,typ="b",ylab="Squared Canonical Correlation",xlab="Canonical variate",main="Scree of canonical correlations" ,ylim=c(0,1),...)},
irt.fa =  {result <- plot.irt(x,labels=labels,...)},

irt.poly = { result <-  plot.poly(x,labels=labels,...)},

fa =  {fa.plot(x,labels=labels,...)},

principal = {fa.plot(x,labels=labels,...)},

parallel = {plot_fa.parallel(x,...)},

reliability = {plot.reliability(x,...)},

residuals = {plot.residuals(x,...)},

none = {stop ("I am sorry, this is not an object I know how to plot")} )


if(!is.null(result))  {class(result) <- c("psych","polyinfo")
   invisible(result)}
}