#modified 6/620  to increase drastically increase speed by taking the drawing out of a loop
#fixed a problem of plotting in wrong order   8/20/20
#perhaps finally fixed the variables to appear in the right order for both lr and vertical mode  09/12/20
bassAckward.diagram <- function(x,digits=2,cut = .3,labels=NULL,marg=c(1.5,.5,1.0,.5),
main="BassAckward",items=TRUE,sort=TRUE,lr=TRUE,curves=FALSE,organize=TRUE,...) {
 old.par<- par(mar=marg)  #give the window some narrower margins
    on.exit(par(old.par))  #set them back

if(organize) x <- ba.organize(x)    
nf = length(x$bass.ack) #this counts how many results are there
if(!items) nf <- nf-1
if(sort){ x$bass.ack[[nf]] <- fa.sort(x$bass.ack[[nf]])
       x$labels[[nf]] <- rownames(x$bass.ack[[nf]]) }
if(lr) {ylim <- c(0,NROW(x$bass.ack[[nf]]))   #this is the number of variables 
xlim <- c(-1,(nf-2)) } else {xlim <- c(0,NROW(x$bass.ack[[nf]]))
ylim <- c(-1,(nf-2))}
lower <- list()
upper <- list()
if(is.null(labels)) labels <- x$labels
labels[[nf]] <- x$labels[[nf]]  #this puts in the bottom row/left hand side (the items) 
plot(0,type="n",xlim=xlim,ylim=ylim,frame.plot=FALSE,axes=FALSE,ylab="",xlab="",main=main)
#first draw the bottom row
nvar <- NROW(x$bass.ack[[nf]])
max.var <- nvar
rname <- labels[[nf]] 
if(lr) {
	all.rects.x <- rep(-1,nvar)
	#all.rects.y <- seq(1:nvar)  #this was not doing it before  (seq(nvar:1) was wronge)
	all.rects.y <- seq(nvar,1, -1 )  #this was not doing it before  (seq(nvar:1) was wronge)
	all.rects.rname <- rname[1:nvar] 
	} else {
 	all.rects.y <- rep(-1,nvar)
	all.rects.x <- seq(1:nvar) 
	#all.rects.rname <- rname[seq(nvar,1,-1)] }
	all.rects.rname <- rname[seq(1, nvar, 1)] }



#first define the various locations but don't draw them
 for(j in 1:nvar) {
   if(lr) {lower [[j]] <- dia.rect(-1, nvar - j + 1, rname[ nvar - j + 1],draw=FALSE,...) } else {lower [[j]] <- dia.rect(j,-1, rname[j],draw=FALSE,...)}
  #  if(lr) {lower [[j]] <- dia.rect(-1,j , rname[j],draw=FALSE,...) } else {lower [[j]] <- dia.rect(j,-1, rname[j],draw=FALSE,...)}  #12/07/20
 
 
 } 
dia.rect(all.rects.x, all.rects.y,all.rects.rname)  #now draw them 

#now draw the next row and then repeat until the top 

for(j in (nf):2) {
	if((j < nf)  & organize) x <- ba.organize(x,j)
	nvar <- NCOL(x$bass.ack[[j]]) 
	scale <- max.var/(nvar+1)

#store the locations
if(lr) {
	all.rects.x <- rep(nf-j,nvar)
	all.rects.y <- seq(nvar,1,-1) * scale
	all.rects.rname <- labels[[j-1]]
	} else {
 	all.rects.y <- rep(nf-j,nvar)
	all.rects.x <- seq(1:nvar) *scale
	all.rects.rname <- labels[[j-1]] }

for(i in 1:nvar) {  #which is actually the number of lower level variables or factors
  cname <- labels[[j-1]]
  
#  if(lr) {upper[[i]] <-  dia.rect(nf-j,(nvar-i + 1) *scale, labels= cname[i],draw=FALSE,...)} else {  upper[[i]] <-  dia.rect(i*scale,nf-j, labels= cname[i],draw=FALSE,...) }
    if(lr) {upper[[i]] <-  dia.ellipse(nf-j,(nvar-i + 1) *scale, labels= cname[i],draw=FALSE,e.size=1,...)} else {  upper[[i]] <-  dia.ellipse(i*scale,nf-j, labels= cname[i],draw=FALSE,e.size=1,...) }
    
    }
#dia.rect(all.rects.x,all.rects.y,all.rects.rname)
 dia.multi.ellipse(all.rects.x,all.rects.y,all.rects.rname) 
  #connect them and then put  in the correlation values

text.values <- list()     #save the text values from the arrows
ki <- 1   #set the counter to 1
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
if(abs(x$bass.ack[[j]][k,i]) >  cut ) { #just draw the large loadings
   value <- x$bass.ack[[j]][k,i]
   
   if(lr) {#text.values[[ki]] <- dia.arrow(upper[[nvar-i +1]]$left,lower[[nfact-k+1]]$right,adj=((i-k) %% 3)   ,labels = round(value,digits),
            #                   col=(sign(value <0) +1),lty=(sign(value<0)+1),draw=FALSE,...)
            text.values[[ki]] <- dia.arrow(upper[[i ]]$left,lower[[k]]$right,adj=((i-k) %% 3)   ,labels = round(value,digits),
                               col=(sign(value <0) +1),lty=(sign(value<0)+1),draw=FALSE,...)
   												
   } else {
     text.values[[ki]] <-  dia.arrow(upper[[i]]$bottom,lower[[k]]$top,adj=((i-k) %% 3)   ,labels = round(value,digits),
                              col=(sign(value <0) +1),lty=(sign(value<0)+1),draw=FALSE,...)}
 
  #text.values <- list(x0,y0,xr,yr,  length = (both+0) * .1*scale, angle = 30, code = 1, xl,yl,xe,ye, length2 = 0.1*scale, angle = 30, code2 = 2)
   ki <- ki +1 
} 
   
}


}

tv <- matrix(unlist(text.values),byrow=TRUE,ncol=21)
text(tv[,1],tv[,2],tv[,3]) # ,tv[,5])    #don't use the adj parameter
      arrows(x0=tv[,6],y0=tv[,7],x1=tv[,8],y1=tv[,9],length=tv[1,10],angle=tv[1,11],code=1,col=tv[,20],lty=tv[,21])
        arrows(x0=tv[,13],y0=tv[,14],x1=tv[,15],y1=tv[,16],length=tv[1,17],angle=tv[1,18],code=2,col=tv[,20],lty=tv[,21])


lower <- upper

}
invisible(x)
}

#organize the lowest two levels to get somewhat cleaner structures

ba.organize <- function(x,level=NULL){
if(is.null(level)) {nf = length(x$bass.ack) #this counts how many results are there
level0 <- fa.sort(x$bass.ack[[nf]])
x$labels[[nf]] <- rownames(level0) 
fa <- x$fa$loadings[[nf-1] ]   #added as fa$loadings to match change in bassAckward
fa <- fa[x$labels[[nf]],]
x$fa[[nf-1] ] <- fa

level1 <- fa.sort(x$bass.ack[[nf-1]]) 
ord1 <- rownames(level1)
level0 <- level0[,ord1,drop=FALSE]
colnames(level0) <- paste0("F",1:NCOL(level0))
x$bass.ack[[nf]] <- level0
x$bass.ack[[nf-1]] <- level1 } else {nf <- level  #just organize the factors, not the items



}
return(x)
}
