#find the correlation between two sets of factors extracted differently

"faCor" <- function(r,nfactors=c(1,1),fm=c("minres","minres"),rotate=c("oblimin","oblimin"),scores=c("tenBerge","tenBerge"), adjust=c(TRUE,TRUE),use="pairwise", cor="cor",weight=NULL,correct=.5,Target=list(NULL,NULL)) {
 cl <- match.call()
#find r if data matrix
    if (!isCorrelation(r)) {  matrix.input <- FALSE  #return the correlation matrix in this case
                       n.obs <- dim(r)[1]  
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
       
  }
 #do the factor 2 different ways 
 if(fm[1]!="pca") {if(is.null(Target[[1]])) {f1 <- fa(r,nfactors=nfactors[1],fm=fm[1],rotate=rotate[1],scores=scores[1])} else {f1 <- fa(r,nfactors=nfactors[1],fm=fm[1],rotate=rotate[1],scores=scores[1],Target=Target[[1]]) }
 
 } else {f1 <- pca(r,nfactors=nfactors[1],rotate=rotate[1])}
 
 if(fm[2]!="pca") {if(is.null(Target[[2]])) {f2 <- fa(r,nfactors=nfactors[2],fm=fm[2],rotate=rotate[2],scores=scores[2])} else {f2 <- fa(r,nfactors=nfactors[2],fm=fm[2],rotate=rotate[2],scores=scores[2],Target=Target[[2]]) }
 } else {f2 <- pca(r,nfactors=nfactors[2],rotate=rotate[2])}
 
 #Find the interfactor correlations

colnames(f1$weights) <- paste0("F",1:ncol(f1$weights))
colnames(f2$weights) <- paste0("F",1:ncol(f2$weights))
rf <-  t(f1$weights) %*% r %*% f2$weights  #adjust by factor variances
rs1 <- diag(t(f1$weights) %*% r %*% f1$weights )
rs2 <- diag(t(f2$weights) %*% r %*% f2$weights )
if(adjust[1]) rf <- diag(1/sqrt(rs1)) %*% rf
if(adjust[2]) rf <- rf %*% diag(1/sqrt(rs2))
rownames(rf) <- colnames(f1$loadings)
colnames(rf) <- colnames(f2$loadings)
fc <- factor.congruence(f1,f2)

result <-list(Call=cl,r=rf,congruence=fc, f1=f1,f2=f2,rs1=rs1,rs2=rs2)
class(result) <- c("psych","faCor")
return(result)
}
 
  
  
