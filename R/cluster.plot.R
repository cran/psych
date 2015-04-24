#revised Sept 16, 2013 to give control over the positon (pos) and size (cex) of the labels
 "cluster.plot" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title="Cluster plot",pch=18,pos,...) {
 if (!is.matrix(ic.results) ) {if (!is.null(class(ic.results)) )   {
  if(class(ic.results)[1] == "kmeans") { load <- t(ic.results$centers) }  else {
      load <-ic.results$loadings} }} else {load <- ic.results}
nc <- dim(load)[2]
nvar <- dim(load)[1]
if(missing(pos)) pos <- rep(1,nvar)   #this allows more control over plotting
ch.col=c("black","blue","red","gray","black","blue","red","gray")
if (is.null(cluster)) {
cluster <- rep(nc+1,nvar)
cluster <- apply( abs(load)  ,1,which.max)
cluster[(apply(abs(load),1,max) < cut)] <- nc+1
}
if (nc > 2 ) {
 pairs(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) }
 else {
 plot(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc <3) text(load,labels,pos=pos,...)
 }
 
 
  "factor.plot" <-  
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title,jiggle=FALSE,amount=.02,pch=18,pos,...) { #deprecated
fa.plot(ic.results,cluster=cluster,cut=cut,labels=labels,title=title,jiggle=jiggle,amount=amount,pch=pch,pos=pos,...)
}


 "fa.plot" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title,jiggle=FALSE,amount=.02,pch=18,pos,...) {
if(missing(title) ) { title="Plot"
              if (length(class(ic.results)) >1 )  {if (class(ic.results)[2] == "fa") {title = "Factor Analysis"} else {
                     if (class(ic.results)[2] == "principal") {title = "Principal Component Analysis"} 
                     } }
                     }
 if (!is.matrix(ic.results)) {
        if (!is.null(class(ic.results))) {
       if(class(ic.results)[1] == "kmeans") { load <- t(ic.results$centers) }  else {
      load <-ic.results$loadings} }} else {load <- ic.results}
      
      
   if(is.null(colnames(load))) colnames(load) <- paste("F",1:ncol(load),sep="")
nc <- dim(load)[2]
nvar <- dim(load)[1]
if(missing(pos)) pos <- rep(1,nvar)   #this allows more control over plotting
ch.col=c("black","blue","red","gray","black","blue","red","gray")
if (is.null(cluster)) {
cluster <- rep(nc+1,nvar)
cluster <- apply( abs(load)  ,1,which.max)
cluster[(apply(abs(load),1,max) < cut)] <- nc+1
}
if(jiggle) load <- jitter(load,amount=amount)
if (nc > 2 ) {
 pairs(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) }
 else {
 plot(load,pch = cluster+pch,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc < 3) text(load,labels,pos=pos,...)
 }
 
 "plot.factor" <-  
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title="Factor plot",...) { #deprecated
 if (!is.matrix(ic.results) ) {if (!is.null(class(ic.results)) )   {
  if(class(ic.results)[1] == "kmeans") { load <- t(ic.results$centers) }  else {
      load <-ic.results$loadings} }} else {load <- ic.results}
   if(is.null(colnames(load))) colnames(load) <- paste("F",1:ncol(load),sep="")
nc <- dim(load)[2]
nvar <- dim(load)[1]
ch.col=c("black","blue","red","gray","black","blue","red","gray")
if (is.null(cluster)) {
cluster <- rep(nc+1,nvar)
cluster <- apply( abs(load)  ,1,which.max)
cluster[(apply(abs(load),1,max) < cut)] <- nc+1
}
if (nc > 2 ) {
 pairs(load,pch = cluster+19,col=ch.col[cluster],bg=ch.col[cluster],main=title) }
 else {
 plot(load,pch = cluster+20,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc <3) text(load,labels,pos=1,...)
 }
 
  "plot.cluster" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title="Cluster plot",...) { #deprecated
 if (!is.matrix(ic.results) ) {if (!is.null(class(ic.results)) )   {
  if(class(ic.results)[1] == "kmeans") { load <- t(ic.results$centers) }  else {
      load <-ic.results$loadings} }} else {load <- ic.results}
nc <- dim(load)[2]
nvar <- dim(load)[1]
ch.col=c("black","blue","red","gray","black","blue","red","gray")
if (is.null(cluster)) {
cluster <- rep(nc+1,nvar)
cluster <- apply( abs(load)  ,1,which.max)
cluster[(apply(abs(load),1,max) < cut)] <- nc+1
}
if (nc > 2 ) {
 pairs(load,pch = cluster+19,cex=1.5,col=ch.col[cluster],bg=ch.col[cluster],main=title) }
 else {
 plot(load,pch = cluster+20,cex=1.5,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc <3) text(load,labels,pos=1)
 }