 "cluster.plot" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title="Cluster plot",...) {
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
 
 
  "factor.plot" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title,jiggle=FALSE,amount=.02,...) {
if(missing(title) ) {if (class(ic.results)[2] == "fa") title = "Factor Analysis"
                     if (class(ic.results)[2] == "principal") title = "Principal Component Analysis"
                     }
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
if(jiggle) load <- jitter(load,amount=amount)
if (nc > 2 ) {
 pairs(load,pch = cluster+19,cex=1.5,col=ch.col[cluster],bg=ch.col[cluster],main=title) }
 else {
 plot(load,pch = cluster+20,cex=1.5,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc <3) text(load,labels,pos=1)
 }
 
 "plot.factor" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title="Factor plot",...) {
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
 pairs(load,pch = cluster+19,cex=1.5,col=ch.col[cluster],bg=ch.col[cluster],main=title) }
 else {
 plot(load,pch = cluster+20,cex=1.5,col=ch.col[cluster],bg=ch.col[cluster],main=title,...) 
 abline(h=0)
 abline(v=0)}
 if(is.null(labels)) labels <- paste(1:nvar)
 if(nc <3) text(load,labels,pos=1)
 }
 
  "plot.cluster" <- 
function(ic.results,cluster=NULL,cut = 0.0,labels=NULL, title="Cluster plot",...) {
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