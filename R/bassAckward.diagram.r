bassAckward.diagram <- function(x,digits=2,cut = .3,labels=NULL,marg=c(1.5,.5,1.0,.5),
main="BassAckward",items=TRUE,sort=TRUE,lr=TRUE,curves=FALSE,organize=TRUE,...) {
 old.par<- par(mar=marg)  #give the window some narrower margins
    on.exit(par(old.par))  #set them back

if(organize) x <- ba.organize(x)    
nf = length(x$bass.ack) #this counts how many results are there
if(!items) nf <- nf-1
if(sort){ x$bass.ack[[nf]] <- fa.sort(x$bass.ack[[nf]])
       x$labels[[nf]] <- rownames(x$bass.ack[[nf]]) }
if(lr) {ylim <- c(0,NROW(x$bass.ack[[nf]]))
xlim <- c(-1,(nf-2)) } else {xlim <- c(0,NROW(x$bass.ack[[nf]]))
ylim <- c(-1,(nf-2))}
lower <- list()
upper <- list()
if(is.null(labels)) labels <- x$labels
plot(0,type="n",xlim=xlim,ylim=ylim,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)
#first draw the bottom row
nvar <- NROW(x$bass.ack[[nf]])
max.var <- nvar
rname <- labels[[nf]] 
 for(j in 1:nvar) {
   if(lr) {lower [[j] ] <- dia.rect(-1,nvar-j +1, rname[j],...) } else {lower [[j] ] <- dia.rect(j,-1, rname[j],...)}
 } 

#now draw the next row and then repeat until the top 

for(j in (nf):2) {
if((j < nf)  & organize) x <- ba.organize(x,j)
nvar <- NCOL(x$bass.ack[[j]]) 
scale <- max.var/(nvar+1)
for(i in 1:nvar) {
  cname <- labels[[j-1]]
  
  if(lr) {upper[[i]] <-  dia.rect(nf-j,(nvar-i + 1) *scale, labels= cname[i],...)} else {  upper[[i]] <-  dia.rect(i*scale,nf-j, labels= cname[i],...) }
    
    }
  
  #connect them

for(i in 1:nvar) {#do it for every top factor
if(length(x$Phi)>0) {Phi <- x$Phi[[j-1]]} else {Phi <- NULL}
nfact  <- NROW(x$bass.ack[[j]])




 if(!is.null(Phi) && (ncol(Phi) >1) && curves) {

   if(i < nvar) {for(k in ((i+1):(nvar))) {
     if(abs(Phi[i,k]) > cut) {
		 if(lr){dia.curve(from=upper[[i]]$right,to=upper[[k]]$right,labels=round(Phi[i,k],digits),scale = .2 , ...) } else {dia.curve(from=upper[[i]]$top,to=upper[[k]]$top,labels=round(Phi[i,k],digits),scale = .2 , ...)}
              }
              }}
              }

for(k in 1:nfact) {
if(abs(x$bass.ack[[j]][k,i]) >  cut ) {
   value <- x$bass.ack[[j]][k,i]
   
   if(lr) {dia.arrow(upper[[i]]$left,lower[[k]]$right,adj=((i-k) %% 3)   ,labels = round(value,digits),
   col=(sign(value <0) +1),lty=(sign(value<0)+1),...)
   												
   } else {
   dia.arrow(upper[[i]]$bottom,lower[[k]]$top,adj=((i-k) %% 3)   ,labels = round(value,digits),
   col=(sign(value <0) +1),lty=(sign(value<0)+1),...)}
 
   
} 
}

}
lower <- upper
}
invisible(x)
}

#organize the lowest two levels to get somewhat cleaner structures

ba.organize <- function(x,level=NULL){
if(is.null(level)) {nf = length(x$bass.ack) #this counts how many results are there
level0 <- fa.sort(x$bass.ack[[nf]])
x$labels[[nf]] <- rownames(level0) 
fa <- x$fa$loadings[[nf-1] ]   #added as fa$loadings to match change if bassAckward
fa <- fa[x$labels[[nf]],]
x$fa[[nf-1] ] <- fa

level1 <- fa.sort(x$bass.ack[[nf-1]]) 
ord1 <- rownames(level1)
level0 <- level0[,ord1]
colnames(level0) <- paste0("F",1:NCOL(level0))
x$bass.ack[[nf]] <- level0
x$bass.ack[[nf-1]] <- level1 } else {nf <- level  #just organize the factors, not the items



}
return(x)
}
