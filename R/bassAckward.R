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
       mixed = {r <- mixedCor(r,use=use,correct=correct)$rho},
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
colnames(pc.weights) <- paste0("C",1:nfactors[1])
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
#two cases:
#normal case is do a regular factor analysis
#but, if rotate = "schmid" we do a schmid leiman rotation 
if(rotate != "schmid") {
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
f <- fn}

}  else {  #do schmid leiman extractions

f <- schmid(r,nfactors[1],fm=fm)    #the top level

f$weights <- solve(r,f$sl[,1:(nfactors[1]+1)]) 
for(nf in 1:num.fac) {
fn <- schmid(r,nfactors[nf],fm=fm)
fn$loadings <- fn$sl[,1:(nfactors[nf] + 1)]
colnames(fn$loadings) <- c("g", paste0("F*",1:(ncol(fn$loadings)-1)))
fa[[nf]] <- fn$loadings
fn$weights <- solve(r,fn$loadings)
colnames(fn$weights) <- c("g",paste0("F*",1:(ncol(fn$weights)-1)))
colnames(f$weights) <- c("g", paste0("F*",1:(ncol(f$weights)-1)))
Phi <- NULL
rf  <-  t(fn$weights) %*% r %*% f$weights  #need to adjust by  variances if not using tenBerge
rs1 <- diag(t(f$weights) %*% r %*% f$weights )
rs2 <- diag(t(fn$weights) %*% r %*% fn$weights )
if(adjust) rf <- (diag(1/sqrt(rs2)) %*% rf %*% diag(1/sqrt(rs1)))
colnames(rf) <- c("g",paste0("F*",1:(ncol(rf)-1)))
rownames(rf) <- c("g", paste0("F*",1:(nrow(rf)-1)))
r.n[[nf]] <- rf
f <- fn
}
}
}

#Now summarize the results 
sumlist <-  sumnames <- labels <- list()
fa.loading.phi <-list()
for(f in 1:nf) {
sumlist[[f]] <- apply(r.n[[f]],2,function(x) {which(max(abs(x))==abs(x))})
sumnames[[f]] <- rownames(r.n[[f]])[sumlist[[f]]]
labels[[f]] <- rownames(r.n[[f]])
if(length(Phi)>0) { #added this check March 3, 2020 
fa.loading.phi [[f]] <-list(loadings = fa[[f]],Phi=Phi[[f]])} else {
fa.loading.phi [[f]] <-list(loadings = fa[[f]],Phi=NA) }

class(fa.loading.phi[[f]]) <- cs(psych,fa)
}
labels[[nf+1]] <- rownames(fn$loadings)
r.n[[nf+1]] <- fn$loadings


result <- list(Call=cl,fm=fm,bass.ack= r.n,Phi=Phi,r = r,summary=sumlist,sumnames=sumnames,labels =labels,fa=fa.loading.phi)
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


