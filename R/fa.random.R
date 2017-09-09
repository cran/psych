"fa.random" <- function(data,nfactors=1,fix=TRUE,n.obs = NA,n.iter=1,rotate="oblimin",scores="regression", residuals=FALSE,SMC=TRUE,covar=FALSE,missing=FALSE,impute="median", min.err = .001,max.iter=50,symmetric=TRUE,warnings=TRUE,fm="minres",alpha=.1, p =.05,oblique.scores=FALSE,np.obs=NULL,use="pairwise",cor="cor",weight=NULL,...) {
subject <- rowMeans(data,na.rm=TRUE)
r <- cor(data,subject,use="pairwise")
colnames(r) <-"within"
data <- data - subject+ fix * rnorm(NROW(data),0,.03)
 f <- fa(r=data,nfactors=nfactors,n.obs=n.obs,rotate=rotate,scores=scores,residuals=residuals,SMC = SMC,covar=covar,missing=missing,impute=impute,min.err=min.err,max.iter=max.iter,symmetric=symmetric,warnings=warnings,fm=fm,alpha=alpha,oblique.scores=oblique.scores,np.obs=np.obs,use=use,cor=cor, weight=weight,...=...) #call fa with the appropriate parameters
 f$subject <- subject
 f$within.r <- r
 return(f)
 }
 

