"bassAckward" <- function(r,nfactors=1,fm="minres",rotate="oblimin",scores="tenBerge",adjust=TRUE,plot=TRUE,cut=.3,use="pairwise", cor="cor",weight=NULL,correct=.5,...) {
 cl <- match.call()
#find r if data matrix
    if (!isCorrelation(r)) {  matrix.input <- FALSE  #return the correlation matrix in this case
                       n.obs <- dim(r)[1]
                       cnames <- colnames(r)
   
    # if given a rectangular matrix, then find the correlation or covariance 
    #multiple ways of find correlations or covariances
    #added the weights option to tet, poly, tetrachoric, and polychoric  June 27, 2018
    switch(cor, 
       cor = {r <- cor(r,use=use)},
       cov = {r <- cov(r,use=use) 
              covar <- TRUE},
       wtd = { r <- cor.wt(r,w=weight)$r},
       tet = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
       poly = {r <- polychoric(r,correct=correct,weight=weight)$rho},
       tetrachoric = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
       polychoric = {r <- polychoric(r,correct=correct,weight=weight)$rho},
       mixed = {r <- mixed.cor(r,use=use,correct=correct)$rho},
       Yuleb = {r <- YuleCor(r,,bonett=TRUE)$rho},
       YuleQ = {r <- YuleCor(r,1)$rho},
       YuleY = {r <- YuleCor(r,.5)$rho } 
       )
      colnames(r) <- rownames(r) <- cnames 
  }

r.n <- list()
fa <- list()
Phi <- list()
num.fac <- length(nfactors)
if (num.fac == 1L) { num.fac <- nfactors 
        nfactors <- 1:num.fac}
  
if(fm =="pca") {
#do the conventional pca bass-akwards approach with components
pc <- pca(r,nfactors[1],rotate=rotate)
#pc <- pca(r)
pc.weights <- pc$weights
colnames(pc.weights) <- paste0("C",nfactors[1])
for(nf in 1:num.fac) {
fn <- pca(r,nfactors[nf],rotate=rotate)
colnames(fn$loadings) <-    paste0("C",1:ncol(fn$loadings))
fa[[nf]] <- fn$loadings
Phi[[nf]] <- fn$Phi
pcn.weights <- fn$weights
colnames(pcn.weights) <- paste0("C",1:ncol(pcn.weights))
colnames(pc$weights) <- paste0("C",1:ncol(pc$weights))
r.n[[nf]] <-  t(pcn.weights) %*% r %*% pc$weights  
pc <- fn
}
} else { #factor analyze
f <- fa(r,nfactors[1],fm=fm,rotate=rotate,scores=scores)
for(nf in 1:num.fac) {
fn <- fa(r,nfactors[nf],fm=fm,rotate=rotate,scores=scores)
colnames(fn$loadings) <-  paste0("F",1:ncol(fn$loadings))
fa[[nf]] <- fn$loadings
Phi[[nf]] <- fn$Phi
fn.weights <- fn$weights
colnames(fn.weights) <- paste0("F",1:ncol(fn$weights))
colnames(f$weights) <- paste0("F",1:ncol(f$weights))

rf  <-  t(fn.weights) %*% r %*% f$weights  #need to adjust by  variances if not using tenBerge
rs1 <- diag(t(f$weights) %*% r %*% f$weights )
rs2 <- diag(t(fn$weights) %*% r %*% fn$weights )
if(adjust) rf <- (diag(1/sqrt(rs2)) %*% rf %*% diag(1/sqrt(rs1)))
colnames(rf) <- paste0("F",1:ncol(rf))
rownames(rf) <- paste0("F",1:nrow(rf))
r.n[[nf]] <- rf
f <- fn
}
}

#Now summarize the results 
sumlist <-  sumnames <- labels <- list()

for(f in 1:nf) {
sumlist[[f]] <- apply(r.n[[f]],2,function(x) {which(max(abs(x))==abs(x))})
sumnames[[f]] <- rownames(r.n[[f]])[sumlist[[f]]]
labels[[f]] <- rownames(r.n[[f]])
}
labels[[nf+1]] <- rownames(fn$loadings)
r.n[[nf+1]] <- fn$loadings

result <- list(Call=cl,fm=fm,bass.ack= r.n,Phi=Phi,r = r,summary=sumlist,sumnames=sumnames,labels =labels,fa=fa)
class(result) <- c("psych","bassAck")
if(plot) bassAckward.diagram(result,cut=cut,...)
return(result)
}

print.psych.back<- function(x,digits=2 ,short=TRUE) {

   cat("\nCall: ")
   print(x$Call)
nf <- length(x$bass.ack)-1
for (f in 1:nf) {
cat("\n",f, 
x$sumnames[[f]])}

if(!short) {
for (f in 1:nf) {
cat("\nFactor correlations\n ")
print(round(x$bass.ack[[f]],digits=digits))}
}
}

summary.psych.back <- function(x,digits=2) {
cat("\nCall: ")
   print(x$Call)
nf <- length(x$bass.ack)-1
for (f in 1:nf) {
cat("\nFactor correlations\n ")
print(round(x$bass.ack[[f]],digits=digits))
}
}


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
Phi <- x$Phi[[j-1]]
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
fa <- x$fa[[nf-1] ]
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
