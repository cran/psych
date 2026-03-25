#convert a CFA analysis output to an IRT output
#March 1, 2026
#finally got scores to work March 22, 2026

irt.CFA <- irtCFA <- function(model=NULL,x=NULL,all=TRUE,n.obs=NA, orthog=FALSE, weight=NULL, correct=0, method="regression", 
	missing=FALSE,impute="none" ,Grice=FALSE,
	plot=TRUE,n.iter=1,
	scores =FALSE) {
	cl <- match.call()
	if(is.matrix(model) ) {x <- model
	          model<- NULL}
	rho <- polychoric(x,correct=correct)
	   n.obs <- rho$n.obs
	   tau <- rho$tau
	   
	cf <- CFA(model=model,x=rho$rho,n.obs=n.obs,all=all,orthog=orthog,Grice=Grice, n.iter=n.iter)
	irt <- fa2irt(cf$cfa,rho,n.obs=n.obs)     #this plots as well
	if(scores) scores <- score.irt.cfa(irt$irt,items = x, X = cf$cfa$model) 
	result <- list(irt=irt,fa=cf,rho=rho,tau=tau,scores=scores,Call=cl)
	
	if(inherits(rho[2],"poly" )) {class(result) <-  c("psych","irt.poly") } else {class(result) <- c("psych","irt.fa")}
class(result) <- c("psych","irt.cfa")
	 return(result)}
	 
	 
score.irt.cfa <- function(stats,items,X, bounds=c(-4,4),mod="logisitic" ) {
nfact <- ncol(X)
scores <- list()
fa.stats <- list()
class(stats) <- c("psych","irt.cfa")
 scores <- scoreIrt(stats,items,keys=X,cut=1,bounds=bounds,mod=mod)
return(scores)
}
 	 

